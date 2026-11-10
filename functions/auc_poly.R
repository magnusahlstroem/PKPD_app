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
    stopifnot(is.data.frame(concentration_frame), MIC > 0,
              dosing_interval > 0, n_interval > 0)

    period_end <- dosing_interval * n_interval
    in_period <- concentration_frame$time >= 0 & concentration_frame$time <= period_end
    x <- concentration_frame$time[in_period]
    y <- concentration_frame$concentration[in_period]

    if (length(x) < 2 || any(!is.finite(x)) || any(!is.finite(y))) {
        stop("concentration_frame must contain at least two finite observations in the dosing period")
    }

    # Use the simulated time grid rather than assuming a fixed resolution of 0.01 h.
    auc <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    auc_24 <- (auc / MIC) * (24 / period_end)

    list(x_poly = c(x, rev(x)),
         y_poly = c(y, rep(0, length(y))),
         AUC = auc_24)
}
