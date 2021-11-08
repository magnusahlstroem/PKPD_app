library(tidyverse)
library(rlang)
library(zoo)

drug_concentration_app <- function(drug_name = "",
                               bacterium = "",
                               VD_kg = 0.19, 
                               t.5 = 1.5,
                               weight = 70, 
                               infusion_time = 0.5, 
                               dose = 1500, 
                               dosing_interval = 8, 
                               n_intervals = 3,
                               bioavailability = 100,
                               meassure_at = "SS",
                               protein_binding = 0.4,
                               resolution = 100,
                               pk_pd_index = "Time>MIC",
                               MIC= 4) {
  #if(bacterium != "") {
  #  MIC_def <- MIC
  #  MIC <- read_csv2(paste("H:/decay_curves/Hjaelpedata/mic_fordelinger/", bacterium, ".csv", sep = "")) %>%
  #    filter(Drugname == drug_name) %>%
  #    pull(ECOFF)
  #  if(is.na(MIC)) {
  #    warning("no e-coff defined for this bug/drug combination, check spelling or make sure to indicate a correct mic. \nOr else the default value of 4 mg/mL is used")
  #    MIC <- MIC_def
  #  }
  #  rm(MIC_def)
  #}
  
  #if(drug_name != "") {
  #  named_args <- lapply(ls(), function(x) get(x))
  #  names(named_args) <- ls()[(ls() != "named_args")]
  #  drug_args <- c(drug_name = drug_name, bacterium = bacterium, get(drug_name))
  #  drug_args <- drug_args[order(names(drug_args))]
  #  formal_args <- fn_fmls()
  #  formal_args <- formal_args[order(names(formal_args))]
  #  named_args <- mapply(function(named_args, formal_args) if(named_args != formal_args) named_args else "", formal_args = formal_args, named_args = named_args, SIMPLIFY = F)
  #  arg_val <- mapply(function(named_args, drug_args) if(named_args != drug_args & named_args != "") named_args else drug_args, drug_args = drug_args, named_args = named_args, SIMPLIFY = F)
  #  #def_val <- def_val[order(names(def_val))]
  #} else {
  #  arg_val <- lapply(ls(), function(x) get(x))
  #  names(arg_val) <- ls()[(ls() != "arg_val")]
  #}
  arg_val <- list(drug_name = drug_name,  bacterium = bacterium, VD_kg = VD_kg, t.5 = t.5,  weight = weight,
                  infusion_time = infusion_time, dose = dose, dosing_interval = dosing_interval, 
                  n_intervals = n_intervals, bioavailability = bioavailability, meassure_at = meassure_at,
                  protein_binding = protein_binding, resolution = resolution, pk_pd_index = pk_pd_index,
                  MIC = MIC)
  
  #calc_ekstra <- min(which(arg_val$dosing_interval*0:6 / (arg_val$t.5*5) > 1)) - 1
  calc_ekstra <- ceiling(arg_val$t.5*5 / arg_val$dosing_interval)
  if(meassure_at == "SD") {
    calc_ekstra <- 0
    arg_val$n_intervals = 1
  }
  if(meassure_at == "SP") {
    calc_ekstra <- 0
  }
  
  expA <- log(1/2)/arg_val$t.5
  VD <- arg_val$VD_kg * arg_val$weight
  time_step <- 1/arg_val$resolution
  time <- seq(0, arg_val$dosing_interval*(arg_val$n_intervals + calc_ekstra), time_step)
  
  
  concentration <- vector()
  plus_conc <- vector()
  minus_conc <- vector()
  for(i in seq_along(time)) {
    if(time[i] == 0) {
      #plus_conc <- dose/(infusion_time * time_step * VD)
      minus_conc <- 0
      concentration <- 0
    } 
    if(time[i] != 0 & time[i] %% arg_val$dosing_interval < arg_val$infusion_time) {
      plus_conc <- (arg_val$dose/(arg_val$infusion_time * arg_val$resolution))/(VD)
      minus_conc <- concentration[i-1] - (concentration[i-1] * exp(expA)^time_step)
      concentration[i] <- concentration[i-1] + plus_conc - minus_conc
    }
    if(time[i] != 0 & (time[i] %% arg_val$dosing_interval >= arg_val$infusion_time | time[i] %% arg_val$dosing_interval == 0)) {
      plus_conc <- 0
      minus_conc <- concentration[i-1] - (concentration[i-1] * exp(expA)^time_step)
      concentration[i] <- concentration[i-1] + plus_conc - minus_conc
    }
  }
  #data.frame(time = time - (calc_ekstra * arg_val$dosing_interval), 
  #           concentration = concentration * (1 - arg_val$protein_binding), 
  #           MIC = arg_val$MIC)


  
  time_conc <- data.frame(time = time - (calc_ekstra * arg_val$dosing_interval), 
                          concentration = concentration * (1 - arg_val$protein_binding) * (arg_val$bioavailability/100), 
                          MIC = arg_val$MIC) %>%
    filter(time >= 0)
  #calc_ekstra
  #time_conc
  
  
  j <- time_conc %>%
    mutate(above_mic = concentration > arg_val$MIC,
           tstart_below = ifelse(!above_mic & is.na(lag(above_mic)), 1, 0),
           tstart_below = ifelse(!above_mic & lag(above_mic) & !is.na(lag(above_mic)), 2, tstart_below),
           tstart_below = ifelse(above_mic & !lead(above_mic) & !is.na(lead(above_mic)), 2, tstart_below),
           
           tstop_below = ifelse(!above_mic & (lead(above_mic) | is.na(lead(above_mic))), 2, 0),
           tstop_below = ifelse(!above_mic & is.na(lead(above_mic)), 1, tstop_below),
           tstop_below = ifelse(above_mic & !lag(above_mic), 2, tstop_below),
           
           tstart_above = ifelse(above_mic & !lag(above_mic), 2, 0),
           tstart_above = ifelse(!above_mic & lead(above_mic) & !is.na(lead(above_mic)), 2, tstart_above),
           
           tstop_above = ifelse(!above_mic & lag(above_mic)  & !is.na(lag(above_mic)), 2, 0),
           tstop_above = ifelse(above_mic & !lead(above_mic) & !is.na(lead(above_mic)), 2, tstop_above),
           tstop_above = ifelse(above_mic & is.na(lead(above_mic)), 1, tstop_above))
  
  
  
  
  if(sum(j$above_mic) != 0) {
    if(sum(j$above_mic) == nrow(j)) {
      j$tstart_above[1] <- 1
      j$tstop_below[1] <- 0
      vars <- c("tstart_above", "tstop_above")
      lty_above <- quote(c(rep(2, length(start_stop[[1]]))))
    } else {
      vars <- c("tstart_below", "tstop_below", "tstart_above", "tstop_above")
      lty_above <- quote(c(rep(2, length(start_stop[[1]])), rep(1, length(start_stop[[3]]))))
    }
    
    #vars <- c("tstart_below", "tstop_below", "tstart_above", "tstop_above")
    start_stop <- lapply(vars, function(x) {
      select(j, time:MIC, x) %>%
        filter(.data[[x]] > 0) %>%
        #arrange(desc(tstart_below)) %>%
        mutate(grouper = if(.data[[x]][1] == 1) ceiling((seq(1,n()) + .5)/2) else ceiling((seq(1,n()))/2)) %>%
        group_by(grouper)  %>%
        summarise(a = if(n() == 2) (concentration[2] - concentration[1])/(time[2] - time[1]) else 1,
                  b = concentration[1] - time[1] * a,
                  time = if(n() == 2) (MIC - b)/a else time) %>%
        filter(!duplicated(time)) %>%
        pull(time)
    })
    
    l <- length(start_stop)
    even <- (1:l)[which((1:l) %% 2 == 0)]
    odd <- (1:l)[which((1:l) %% 2 == 1)]
    
    
    #above_mic_frame <- data.frame(tstart = c(start_stop[[1]], start_stop[[3]]), 
    #                              tstop = c(start_stop[[2]], start_stop[[4]]), 
    #                              lty = c(rep(2, length(start_stop[[1]])), rep(1, length(start_stop[[3]]))))

    above_mic_frame <- data.frame(tstart = do.call(c, start_stop[odd]), 
                                  tstop = do.call(c, start_stop[even]) , 
                                  lty = eval(lty_above))
    
    time_above_mic <- group_by(above_mic_frame, lty) %>%
      summarise(time_above_mic = sum(tstop - tstart)) %>%
      ungroup() %>%
      summarise(time_above_mic = time_above_mic[1]/sum(time_above_mic)) %>%
      pull(time_above_mic)
    
    if(arg_val$pk_pd_index == "24-AUC") {
      auc_poly <- auc_poly(time_conc, arg_val$MIC, arg_val$dosing_interval, arg_val$n_interval)
    } else {
      auc_poly <- list()
    }
    
    
  } else {
    above_mic_frame <- data.frame(tstart = 0, 
                                  tstop = arg_val$n_intervals * arg_val$dosing_interval, 
                                  lty = 2)
    time_above_mic = 0
    auc_poly <- list()
  }
  
  peak_conc <- list(time = arg_val$infusion_time, 
                    concentration = filter(time_conc, time ==  arg_val$infusion_time) %>% pull(concentration), 
                    Cmax_MIC = (filter(time_conc, time ==  arg_val$infusion_time) %>% pull(concentration))/arg_val$MIC)
  
  
  out <- list(drug_name = drug_name, bacterium = bacterium, time_conc = time_conc, above_mic_frame = above_mic_frame, time_above_mic = time_above_mic, auc_poly = auc_poly, arg_val = arg_val, peak_conc = peak_conc)
  class(out) <- append(class(out), "drug_conc")
  out
}



