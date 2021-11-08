library(shiny)
library(readr)
library(tidyverse)
library(DT)

#setwd("H:/decay_curves/pkpd_app")
Bakt <- list.files("Hjaelpedata/mic_distributions")
Antibiotika <- read_csv2("Hjaelpedata/other/pkpd.csv") %>%
  mutate_at("drug_name", function(x) gsub("_", "-", x))
#mapply(function(name, list) assign(name, list, envir = globalenv()), name = Antibiotika$drug_name, list = apply(Antibiotika, 1, function(x) as.list(x[-1])), SIMPLIFY = F)
setwd("functions")
f <- list.files()
lapply(f, source)
setwd("..")

ui <- shinyUI(fluidPage(
  
  titlePanel("Bakterier og antibiotika følsomhed"),
  
  htmlOutput("output_table"),
  plotOutput("conc_plot"),
  
  fluidRow(
    column(3,
           
           selectizeInput(
             'Bakterie', 'Bakterie', choices = sub(".csv$", "", Bakt),
             selected = "Escherichia coli",
             #multiple = T,
             options = list(
               placeholder = 'Søg efter bakterie',
               onInitialize = I('function() { this.setValue(""); }')
             )
           ),
           selectizeInput(
             'Antibiotikum', 'Antibiotikum', choices = pull(Antibiotika, drug_name),
             selected = "Cefuroxime",
             #multiple = T,
             options = list(
               placeholder = 'Søg efter antibiotikum',
               onInitialize = I('function() { this.setValue(""); }')
             )
           ),
           selectInput(
             'MIC', 'MIC (mg/L)', 
             c("0,002", "0,004", "0,008", 
               "0,016", "0,03", "0,06", 
               "0,125", "0,25", "0,5", "1", "2", "4",
               "8", "16", "32", "64", "128", "256", "512"), 
             selected = "4"
           )
    ),
    column(3,
           selectInput(
             'pk_pd_index', 'PK/PD Index', choices = c("24-AUC", "Time>MIC", "Cmax-MIC ratio"), selected = "Time>MIC"
           ),
           numericInput(
             'protein_binding', 'Proteinbinding (%)', value = 40, min = 0, max = 100, step = 1
           ),
           numericInput(
             't.5', 'Plasma halveringstid (t)', value = 0.5, min = 0.1, max = 250, step = 0.1
           ),
           numericInput(
             'VD_kg', 'Fordelingsvolumen (L/kg)', value = 0.5, min = 0, max = 10, step = 0.01
           ),
           numericInput(
             'bioavailability', 'Biotilgængelighed (%)', value = 100, min = 0, max = 100, step = 1
           )
    ),
    column(3,
           numericInput(
             'infusion_time', 'Infusiontid (t)', value = 0.5, min = 0.1, max = 48, step = 0.1
           ),
           numericInput(
             'dose', 'Stofmængde pr. administration (mg)', value = 1500, min = 0, max = 50000, step = 1
           ),
           numericInput(
             'n_intervals', 'Antal administrationer', value = 3, min = 1, max = 100, step = 1
           ),
           numericInput(
             'dosing_interval', 'Tid mellem administrationer (t)', value = 8, min = 0.1, max = 48, step = 0.1
           )
    ),
    column(3,
           numericInput(
             'weight', 'Vægt (kg)', value = 70, min = 0, max = 500, step = 1
           ),
           selectInput(
             'formulation', 'Administrationmåde', choices = c("iv", "po")
           ),
           selectInput(
             'meassure_at', 'Single dose or steady state', choices = c("SS", "SD", "SP")
           )
    )
  )
)
)

server <- shinyServer(function(input, output) {
  
  ab <- reactive({
    filter(Antibiotika, drug_name == input$Antibiotikum)  
  })
  
  bakterie <- reactive({
    if(input$Bakterie != "") {
      read_csv2(paste("Hjaelpedata/mic_distributions/", input$Bakterie, ".csv", sep = ""),
                col_types = "cnnnnnnnnnnnnnnnnnnnnncccc") %>%
        mutate_at(c("ECOFF", "Sensitive", "Resistant"), function(x) sub(",0*$", "", x)) %>%
        mutate_at( c("ECOFF", "Sensitive", "Resistant"), function(x) sub("0*$", "", x)) %>%
        mutate_at("Drugname", function(x) gsub("_", "-", x))
    } else {
      read_csv2("Hjaelpedata/mic_distributions/Escherichia coli.csv",
                col_types = "cnnnnnnnnnnnnnnnnnnnnncccc") %>%
        mutate_at(c("ECOFF", "Sensitive", "Resistant"), function(x) sub(",0*$", "", x)) %>%
        mutate_at( c("ECOFF", "Sensitive", "Resistant"), function(x) sub("0*$", "", x)) %>%
        mutate_at("Drugname", function(x) gsub("_", "-", x))
    }
  })
  
  for_mic <- reactive({
    if(input$Bakterie != "" & input$Antibiotikum != "") {
      j <- bakterie() %>% filter(Drugname == input$Antibiotikum) %>%
        select(starts_with("MIC")) %>%
        as.numeric()
      
      if(!sum(is.na(j)) > 0) {
        j2 <- cumsum(j)
        frac_j2 <- j2/j2[length(j2)]
        ind <- min(which(frac_j2 > .90))
        MIC_txt <- names(select(bakterie(), starts_with("MIC")))[ind+1]
        dist_mic <- #as.numeric(
          sub("MIC_", "", MIC_txt)
        #)
      } else {
        dist_mic <- NA
      }
      
      bind_cols(bakterie() %>% 
                  filter(Drugname == input$Antibiotikum) %>%
                  select(Sensitive, Resistant, ECOFF), dist_mic = sub("\\.", ",", as.character(dist_mic)))
    } else {
      bind_cols(Sensitive = 4, Resistant = 4, ECOFF = 4, dist_mic = 4)  
    }
  })
  
  ##Rettes noget omkring MIC breakpoints
  #MIC <- reactive({
  #  if(input$Bakterie != "") {
  #    read_csv2(paste("hjaelpedata/mic_distributions/", input$Bakterie, ".csv", sep = "")) %>%
  #      filter(Drugname == input$Antibiotikum) %>%
  #      pull(ECOFF)
  #  } else {
  #    4
  #  }
  #})
  
  output$conc_plot <- renderPlot({
    if(!is.na(input$VD_kg) & !is.na(input$t.5) & !is.na(input$weight) & !is.na(input$infusion_time) &
       !is.na(input$dose) & !is.na(input$dosing_interval) & !is.na(input$n_intervals) & 
       !is.na(input$bioavailability) & !is.na(input$meassure_at) & !is.na(input$protein_binding) &
       !is.na(input$pk_pd_index) & !is.na(input$MIC)) {
      plot(output_object()) #output_object()
    } else {
      plot(1,1)
    }
    
  })
  
  output$output_table <- renderText({
    if(!is.na(input$VD_kg) & !is.na(input$t.5) & !is.na(input$weight) & !is.na(input$infusion_time) &
       !is.na(input$dose) & !is.na(input$dosing_interval) & !is.na(input$n_intervals) & 
       !is.na(input$bioavailability) & !is.na(input$meassure_at) & !is.na(input$protein_binding) &
       !is.na(input$pk_pd_index) & !is.na(input$MIC)) {
      txt_bact <- input$Bakterie
      txt_ab <- input$Antibiotikum
      txt_pkpd_index_choice <- switch((ab() %>% pull(pk_pd_index) == input$pk_pd_index) + 1,
                                      " baseret på dit valg som tilsidesætter standardindstillingen",
                                      " baseret på standardindstilling")
      txt_pkpd_index <- switch(input$pk_pd_index,
                               "Time>MIC" = "Tid over MIC",
                               "24-AUC" = "AUC/MIC")
      
      txt_pkpd_index_value <- if(input$pk_pd_index == "Time>MIC") { 
        format(output_object()$time_above_mic*100, digits = 3) 
      } else {
        format(output_object()$auc_poly$AUC, digits = 3)
      }
      
      txt_mic_choice <- switch(min(which(!is.na(as.numeric(sub(",", ".", for_mic()))) & for_mic() != "0,001")), 
                               "det kliniske brydepunkt for <font color='#F0027F'><b>susceptible, standard dosing regimen</b></font>", 
                               "det kliniske brydepunkt for <font color='#F0027F'><b>susceptible, increased exposure</b></font>", 
                               "<font color='#F0027F'><b>det epidemiologiske cut-off (ECOFF)</b></font>", 
                               "<font color='#F0027F'><b>MIC fordelingen med et 90% cut-off baseret på EUCAST data</b></font>")
      txt_mic <- input$MIC
      
      paste("Du har valgt bakterien <font color='#F0027F'><b>", em(input$Bakterie), "</b></font>, og antibiotikummet <font color='#F0027F'><b>", txt_ab, 
            "</b></font><br/>Det valgte PKPD-index er <font color='#F0027F'><b>", txt_pkpd_index, ",", txt_pkpd_index_choice, 
            "</b></font><br/>Den valgte minimale inhibitoriske concentration (MIC) er <font color='#F0027F'><b>", txt_mic, " &microg/mL</b></font> baseret på ", txt_mic_choice, 
            " <br/>På baggrund af de ovenstående parametre samt de øvrige parametre der er valgt nedenfor, ",
            "har vi estimeret en <font color='#F0027F'><b>", txt_pkpd_index, "</b></font> på: <font color='#F0027F'><b>", txt_pkpd_index_value, "%</b></font>", sep = "") #output_object()
    } else {
      "Vælg parametre nedenfor. Start eksempelvis med at vælge en bakterie og et antibiotikum."
    }
    
    #for_mic()$Sensitive, for_mic()$Resistant, for_mic()$ECOFF, for_mic()$dist_mic,
    #      min(which(!is.na(as.numeric(sub(",", ".", for_mic()))))))
    #paste(class(input$Antibiotikum),
    #      class(input$Bakterie),for_mic()
    #      class(input$VD_kg), 
    #      class(input$t.5),
    #      class(input$weight), 
    #      class(input$infusion_time), 
    #      class(input$dose), 
    #      class(input$dosing_interval), 
    #      class(input$n_intervals),
    #      class(input$bioavailability),
    #      class(input$protein_binding),
    #      100,
    #      class(input$pk_pd_index),
    #      class(input$MIC), sep = "_")
  })
  
  #output$bakterie <- renderTable({
  #  ab()
  #})
  
  observeEvent(input$Bakterie, {
    updateSelectizeInput(inputId = "Antibiotikum", 
                         choices = intersect(pull(Antibiotika, drug_name), pull(bakterie(), Drugname)))
  })
  
  #observeEvent(input$Antibiotikum, {
  #  updateNumericInput(inputId = "MIC", 
  #                     value = if(input$Bakterie != "") {
  #                       bakterie() %>%
  #                         filter(Drugname == input$Antibiotikum) %>%
  #                         pull(Breakpoint)
  #                     } else {
  #                       "4"
  #                     })
  #})
  observeEvent(input$Antibiotikum, {
    updateTextInput(inputId = "MIC", 
                       value = switch(min(which(!is.na(as.numeric(sub(",", ".", for_mic()))) & for_mic() != "0,001")), 
                                      for_mic()$Sensitive, 
                                      for_mic()$Resistant,
                                      for_mic()$ECOFF, 
                                      for_mic()$dist_mic)
                       )
    
  })
  
  #observeEvent(input$Bakterium, {
  #  updateNumericInput(inputId = "MIC", 
  #                     value = if(input$Antibiotikum != "") {
  #                       bakterie() %>%
  #                         filter(Drugname == input$Antibiotikum) %>%
  #                         pull(Breakpoint)
  #                     } else {
  #                       "4"
  #                     })
  #})
  observeEvent(input$Bakterium, {
    updateTextInput(inputId = "MIC", 
                       value = switch(min(which(!is.na(as.numeric(sub(",", ".", for_mic()))) & for_mic() != "0,001")), 
                                      for_mic()$Sensitive, 
                                      for_mic()$Resistant, 
                                      for_mic()$ECOFF, 
                                      for_mic()$dist_mic)
                       )
  })
  
  observeEvent(input$Antibiotikum, {
    updateSelectInput(inputId = "pk_pd_index", 
                      selected = filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(pk_pd_index))
  })
  observeEvent(input$Antibiotikum, {
    updateNumericInput(inputId = "dose", 
                       value = filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(dose))
  })
  observeEvent(input$Antibiotikum, {
    updateNumericInput(inputId = "n_intervals", 
                       value = filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(n_intervals))
  })
  observeEvent(input$Antibiotikum, {
    updateNumericInput(inputId = "dosing_interval", 
                       value = filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(dosing_interval))
  })
  observeEvent(input$Antibiotikum, {
    updateNumericInput(inputId = "protein_binding", 
                       value = filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(protein_binding))
  })
  observeEvent(input$Antibiotikum, {
    updateNumericInput(inputId = "t.5", 
                       value = filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(t.5))
  })
  observeEvent(input$Antibiotikum, {
    updateNumericInput(inputId = "VD_kg", 
                       value = filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(VD_kg))
  })
  observeEvent(input$Antibiotikum, {
    updateNumericInput(inputId = "formulation", 
                       value = filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(formulation))
  })
  observeEvent(input$Antibiotikum, {
    updateNumericInput(inputId = "meassure_at", 
                       value = filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(meassure_at))
  })
  observeEvent(input$formulation, {
    updateNumericInput(inputId = "bioavailability", 
                       value = if(input$formulation == "po") {
                         filter(Antibiotika, drug_name == input$Antibiotikum) %>% pull(bioavailability) 
                       } else {
                           100
                         }
                       )
  })
  
  output_object <- reactive({
    drug_concentration_app(drug_name = input$Antibiotikum,
                           bacterium = input$Bakterie,
                           VD_kg = input$VD_kg, 
                           t.5 = input$t.5,
                           weight = input$weight, 
                           infusion_time = input$infusion_time, 
                           dose = input$dose, 
                           dosing_interval = input$dosing_interval, 
                           n_intervals = input$n_intervals,
                           bioavailability = input$bioavailability,
                           meassure_at = input$meassure_at,
                           protein_binding = input$protein_binding,
                           resolution = 100,
                           pk_pd_index = input$pk_pd_index,
                           MIC = as.numeric(sub(",", ".", input$MIC)))
  }) 
  
}) 


shinyApp(ui = ui, server = server)