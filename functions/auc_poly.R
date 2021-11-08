#auc_poly <- function(above_mic_frame, concentration_frame, MIC, dosing_interval, n_interval) {
#    above_mic_frame <- filter(above_mic_frame, lty == 1)
#   
#   out <- mapply(function(tstart, tstop) {
#        first_whole <- ceiling((tstart*100))/100
#        last_whole <- floor((tstop*100))/100
#        
#        ind_first_conc <- which(concentration_frame$time == first_whole)
#        ind_last_conc <- which(concentration_frame$time == last_whole)
#        
#        first_conc <- concentration_frame$concentration[ind_first_conc]
#        last_conc <- concentration_frame$concentration[ind_last_conc]
#        
#        x_poly <- c(tstart, seq(first_whole, last_whole, 0.01), tstop, tstart)
#        y_poly <- c(MIC, concentration_frame$concentration[ind_first_conc:ind_last_conc], MIC, MIC)
#        id <- order(x_poly)
#        
#        AUC <- sum(diff(x_poly[id])*rollmean(y_poly[id] - MIC,2))/(dosing_interval * n_interval)
#        
#        list(x_poly = x_poly, y_apoly = y_poly, AUC = AUC)
#    },
#    tstart = above_mic_frame$tstart,
#    tstop = above_mic_frame$tstop,
#    SIMPLIFY = F
#    )
#    out
#}


auc_poly <- function(concentration_frame, MIC, dosing_interval, n_interval) {
    #above_mic_frame <- filter(above_mic_frame, lty == 1)
    first_whole <- 0
    last_whole <- dosing_interval*n_interval
    
    ind_first_conc <- which(concentration_frame$time == first_whole)
    ind_last_conc <- which(concentration_frame$time == last_whole)
    
    first_conc <- concentration_frame$concentration[ind_first_conc]
    last_conc <- concentration_frame$concentration[ind_last_conc]
    
    x_poly <- c(seq(0,dosing_interval*n_interval,0.01), dosing_interval*n_interval, 0)
    y_poly <- c(concentration_frame$concentration, 0, 0)
    id <- order(x_poly)
    
    AUC <- (sum(diff(x_poly[id])*rollmean(y_poly[id],2))/MIC)*(24/(dosing_interval*n_interval))
    
    list(x_poly = x_poly, y_poly = y_poly, AUC = AUC)
}