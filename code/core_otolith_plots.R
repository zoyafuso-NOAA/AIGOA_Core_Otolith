###############################################################################
## Project:       
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov), modfied from 
##                Lewis Barnett (lewis.barnett@noaa.gov)
##
## Notes:         
###############################################################################
rm(list = ls())

###############################################################################
####   Import packages
############################################################################### 
library(readxl)

###############################################################################
####   Load bootstrap results and collection spreadsheet
###############################################################################
load("results/bootstrap_results.RData")
collection <- readxl::read_xlsx(path = "data/AI_core_collections.xlsx", 
                                sheet = "2022 collection") 
collection <- subset(x = collection, subset = total_requested > 0)

###############################################################################
####   Plot boxplot time series of total bootstrapped collected otoliths 
####       for each species
###############################################################################
{
  ## Open Device
  pdf(file = "results/species_bootstrap.pdf",
     width = 8.5, height = 11,
     family = "serif", version = "1.7")
  
  ## Plot layout
  par(mfrow = c(6, 3), mar = c(2, 3, 3, 1), oma = c(5, 5, 3, 5))
  
  for (irow in 1:nrow(collection) ) { ## Loop over species -- start
    
    ## What species code, common name, and target?
    species_name <- collection$common_name[irow]
    target_requested <- !is.na(collection$total_requested[irow])
    target_n <- collection$total_requested[irow]

    ## Calculate quantiles of the total otoliths collected
    temp_dist <- apply(X = bootstrap_oto_main[, irow, ], 
                       MARGIN = 1, 
                       FUN = quantile, 
                       probs = c(0.025, 0.25, 0.5, 0.75, 0.975) )
    
    if(is.na(target_n)) target_n <- max(temp_dist["50%", ])
    
    ## Base plot
    plot(temp_dist["50%", ], type = "n", axes = F, ann = F,
         ylim = range(c(c(1.25, 0.75) * target_n, 
                        temp_dist[c("2.5%", "97.5%"), ])) )
    
    ## Plot boxplot time series
    polygon(x = c(1:5, 5:1), 
            y = c(temp_dist["2.5%", ], rev(temp_dist["97.5%", ])),
            col = "cornflowerblue", border = F)
    
    polygon(x = c(1:5, 5:1), 
            y = c(temp_dist["25%", ], rev(temp_dist["75%", ])),
            col = "darkgreen", border = F)
    
    lines(x = 1:5, y = temp_dist["50%", ], lwd = 2)
    points(x = 1:5, y = temp_dist["50%", ], pch = 16, cex = 1.5)
    
    ## Add target total collection
    if(target_requested) abline(h = target_n, lty = 'dashed')
    
    ## Add observed total collection
    # points(x = 1:5, y = total_ages[species_name, ], 
    #        pch = "*", col = "red", cex = 2)
    
    ## Axis Labels
    box()
    axis(side = 1, at = 1:5, labels = c(2010, 2012, 2014, 2016, 2018))
    axis(side = 2, las = 1)
    mtext(side = 3, text = species_name, line = 0.5)
  }  ## Loop over species -- close
   
  ## Add Legend
  plot(1, type = "n", axes = F, ann = F, xlim = c(0, 10), ylim = c(0, 100))
  rect(xleft = 4, xright = 6, ybottom = 5, ytop = 95, 
       col = "cornflowerblue", border = "cornflowerblue")
  rect(xleft = 4, xright = 6, ybottom = 25, ytop = 75, 
       col = "darkgreen", border = "darkgreen")
  points(x = 5, y = 50, pch = 16)
  text(x = 6, y = c(5, 25, 50, 75, 95), 
       labels = paste0(c(2.5, 25, 50, 75, 97.5), "%"), pos = 4)
  mtext(side = 3, text = "percentile")
  
  ## Add main panel axis labels
  mtext(side = 1, text = "Year", font = 2, outer = T)
  mtext(side = 2, text = "Total Otoliths Collected", font = 2, outer = T)
  
  ############################################################################
  ####   Plot boxplot time series of total per-haul otoliths
  ####   
  ####   For each bootstrapped sample, you have a calculation of the quantiles
  ####      (2.5, 25, 50, 75, 97.5%) of total number of otoliths collected 
  ####      across the 420 stations. For each quantile, there is a distribution
  ####      of 1000 values. Calculate the medians of those distributions
  ############################################################################
  total_work <- apply(X = bootstrap_total_work,
                      MARGIN = c(1, 2),
                      FUN = median)
  
  ###############################################################################
  ####   
  ############################################################################### ## Base plot
  plot(1, type = "n", axes = F, ann = F,
       xlim = c(1, 5),
       ylim = range(total_work) )
  
  ## Plot boxplot time series
  polygon(x = c(1:5, 5:1), 
          y = c(total_work[, "2.5%"], rev(total_work[, "97.5%"])),
          col = "cornflowerblue", border = F)
  
  polygon(x = c(1:5, 5:1), 
          y = c(total_work[, "25%"], rev(total_work[, "75%"])),
          col = "darkgreen", border = F)
  
  lines(x = 1:5, y = total_work[, "50%"], lwd = 2)
  points(x = 1:5, y = total_work[, "50%"], pch = 16, cex = 1.5)
  
  ## Axis Labels
  box()
  axis(side = 1, at = 1:5, labels = c(2010, 2012, 2014, 2016, 2018))
  axis(side = 2, las = 1)
  mtext(side = 3, text = "Haul-Level total otoliths", line = 0.5)
  
  dev.off()
}


