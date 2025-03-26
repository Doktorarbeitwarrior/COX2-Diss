library(forestplot)
library(dplyr)
# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <- structure(list(mean  = c(NA, NA,  4.67, 4.33, 2.9,  NA, NA), 
                                      lower = c(NA, NA, 2.37, 1.866, 1.63,  NA, NA),
                                      upper = c(NA, NA,  9.206, 9.954, 5.16,  NA, NA)),
                                 .Names = c("mean", "lower", "upper"), 
                                 row.names = c(NA, -10L), 
                                 class = "data.frame")

tabletext <- cbind(c("", "RiskModel/Endpoint",  "OS/OS AY MALE", "OS/DSS AY MALE","OS/PFI AY MALE",NA, NA),
                   c("", "pValue","8,51E-06","5,5E-04", "2,9E-04",NA, NA),
                   c("", "HR", "4,67","4,33","2,9", NA, NA))

cochrane_from_rmeta %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(rep(TRUE, 2), rep(FALSE, 10), TRUE),
             clip = c(0, 9.99), 
             xlog = FALSE,
             boxsize = 0.2,
             txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks = gpar(cex=1)),
             xlab = "HazardRatio",
             graphwidth = unit(70,"mm"),
             col = fpColors(box = "black",
                            line = "darkblue",
                            summary = "royalblue"))