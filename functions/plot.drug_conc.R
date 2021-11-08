plot.drug_conc <- function(drug_conc, 
                           incl_mic_line = T,
                           incl_info_table = F,...) {
    mytitle = drug_conc$drug_name
    mysubtitle = paste("vs. ", drug_conc$bacterium, sep = "")
    ylim_max <- if(drug_conc$arg_val$MIC > max(drug_conc$time_conc$concentration)) drug_conc$arg_val$MIC * 1.05 else max(drug_conc$time_conc$concentration) * 1.05 
    
    if(drug_conc$arg_val$pk_pd_index == "Time>MIC") {
        legend_text <- paste("MIC = ", drug_conc$arg_val$MIC, " mg/L",
              "\nTime > MIC = ",  format(drug_conc$time_above_mic * 100, digits = 0), "%",
              "\nHalf-life = ", drug_conc$arg_val$t.5, "h",
              "\nDosing interval = ", drug_conc$arg_val$dosing_interval,
              "\nDrug pr. dose = ",  drug_conc$arg_val$dose, "mg",
              "\nInfussion time = ", drug_conc$arg_val$infusion_time, "h",
              "\nVD (L/kg) = ", drug_conc$arg_val$VD_kg,
              "\nProtein binding = ", drug_conc$arg_val$protein_binding * 100, "%",
              "\nWeight = ", drug_conc$arg_val$weight, "kg",
              sep = "")
    } 
    if(drug_conc$arg_val$pk_pd_index == "24-AUC") {
        legend_text <- paste("MIC = ", drug_conc$arg_val$MIC, " mg/L",
                             "\nfAUC/MIC = ",  format(drug_conc$auc_poly$AUC, digits = 2),
                             "\nHalf-life = ", drug_conc$arg_val$t.5, "h",
                             "\nDosing interval = ", drug_conc$arg_val$dosing_interval,
                             "\nDrug pr. dose = ",  drug_conc$arg_val$dose, "mg",
                             "\nInfussion time = ", drug_conc$arg_val$infusion_time, "h",
                             "\nVD (L/kg) = ", drug_conc$arg_val$VD_kg,
                             "\nProtein binding = ", drug_conc$arg_val$protein_binding * 100, "%",
                             "\nWeight = ", drug_conc$arg_val$weight, "kg",
                             sep = "")
    }
    if(drug_conc$arg_val$pk_pd_index == "Cmax-MIC ratio") {
        legend_text <- paste("MIC = ", drug_conc$arg_val$MIC, " mg/L",
                             "\nCmax/MIC ratio = ",  format(drug_conc$peak_conc$Cmax_MIC, digits = 2),
                             "\nHalf-life = ", drug_conc$arg_val$t.5, "h",
                             "\nDosing interval = ", drug_conc$arg_val$dosing_interval,
                             "\nDrug pr. dose = ",  drug_conc$arg_val$dose, "mg",
                             "\nInfussion time = ", drug_conc$arg_val$infusion_time, "h",
                             "\nVD (L/kg) = ", drug_conc$arg_val$VD_kg,
                             "\nProtein binding = ", drug_conc$arg_val$protein_binding * 100, "%",
                             "\nWeight = ", drug_conc$arg_val$weight, "kg",
                             sep = "")
    }
        
    #max(drug_conc$time_conc$concentration)
#}
    plot(drug_conc$time_conc$concentration ~ drug_conc$time_conc$time, 
         type = "l",
         ylab = "Fri stof koncentration (mg/L)",
         xlab = "Tid (t)",
         #xlim = c(0,drug_conc$arg_val$dosing_interval*drug_conc$arg_val$n_intervals * 1.28),
         ylim = c(0,ylim_max))
    mtext(side=3, line=2, at=-0.07, adj=0, cex=1.3, mytitle)
    mtext(side=3, line=1, at=-0.07, adj=0, cex=1, mysubtitle)
    if(drug_conc$arg_val$pk_pd_index == "Cmax-MIC ratio") {
        lines(x = rep(drug_conc$arg_val$infusion_time, 2), 
              y = c(0, drug_conc$peak_conc$concentration))
        lines(x = c(-0.5, 0.5) + drug_conc$arg_val$infusion_time, y = rep(drug_conc$peak_conc$concentration, times = 2))
    }
    if(incl_mic_line) {
        if(drug_conc$arg_val$pk_pd_index == "24-AUC") {
            polygon(drug_conc$auc_poly$x_poly, drug_conc$auc_poly$y_poly, 
                    #density = 50,
                    col = "#7FC97F")
            #, density = 15)
        }
        ##FARVER
        mapply(function(x1, x2, y, lty) lines(x = c(x1, x2), y = c(y, y), lty = lty, col = "#666666"), 
               x1 = drug_conc$above_mic_frame$tstart,
               x2 = drug_conc$above_mic_frame$tstop, 
               y = drug_conc$arg_val$MIC, 
               lty = drug_conc$above_mic_frame$lty)
        
            
    }
    if(incl_info_table) {
        legend(drug_conc$arg_val$dosing_interval*drug_conc$arg_val$n_intervals * 0.98, 
               max(drug_conc$time_conc$concentration)*1.05, 
               legend = legend_text, 
               cex = 0.8,
               #yjust = 0.5,
               adj = c(0, 0.05)
               , box.col = "white"
               )
    }
}