t.5_egfr <- function(eGFR, Hep_CL, t.5, t.5_esrd, liverfunction, VD) {
    if(missing(eGFR) | missing(Hep_CL) | missing(t.5) | missing(t.5_esrd) | missing(VD)) {
        t.5_est <- 0.5 
    } else  if(!is.na(Hep_CL) & !is.na(t.5_esrd) | (is.na(Hep_CL) & !is.na(t.5_esrd))) {
        CL_total_egfr80 <- 0.693 * VD / t.5
        CL_other <- 0.693 * VD / t.5_esrd
        CL_part_other <- CL_other/CL_total_egfr80
        CL_part_urin <- (CL_total_egfr80 - CL_other)/CL_total_egfr80
        other_CL <- CL_other * CL_total_egfr80
        urin_CL_egfr80 <- CL_part_urin * CL_total_egfr80
        t.5_est <- 0.693 * VD /(CL_other * liverfunction + urin_CL_egfr80 * eGFR/80)
    } else if(is.na(t.5_esrd)) {
        #if(Hep_CL > 1 | Hep_CL < 0 | is.na(Hep_CL)) {
        #    warning("Hep_CL must be numeric between 0 and 1") 
        #} else {
            CL_total_egfr80 <- 0.693 * VD / t.5
            CL_other <- Hep_CL * CL_total_egfr80
            CL_part_other <- CL_other/CL_total_egfr80
            CL_part_urin <- 1-Hep_CL
            other_CL <- CL_other * CL_total_egfr80
            urin_CL_egfr80 <- CL_part_urin * CL_total_egfr80
            t.5_est <- 0.693 * VD /(CL_other * liverfunction + urin_CL_egfr80 * eGFR/80)
        #}
    } else if(is.na(Hep_CL) & is.na(t.5_esrd)) {
        t.5_est <- t.5
    } else {
        t.5_est <- t.5
    }
    if(is.na(t.5_est)) 1.5 else round(t.5_est, 2)
}
