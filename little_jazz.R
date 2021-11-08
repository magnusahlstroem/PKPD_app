Vd <- 0.19
t.5_normal <- 1.5
t.5_esrd <- 17
weight <- 70
Cl_total_normal <- 0.693 * Vd * weight/t.5_normal
Cl_liver <- 0.693 * Vd * weight/t.5_esrd
Cl_total <- Cl_liver + Cl_renal

t.5 <- 0.693*Vd/Cl