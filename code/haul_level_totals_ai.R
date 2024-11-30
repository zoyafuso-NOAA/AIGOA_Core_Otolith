##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Optimize otolith collection rules
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For the AI/GOA BTS, calculate optimal otolith collection rules
##                    across species. 
##                Constraints are: total otoliths requested, distribution 
##                    and abundance of species, median haul-level collection 
##                    around 20 and avoid haul-level collections that are > 30
##                    (~ 10% of hauls are okay). 
##                Bootstrap hauls to add some variation to the analysis.
##
## Notes:         If downloading data for the first time, use the sumfish
##                   package. In the future, this may be replaced by the
##                   AFSC.GAP.DBE R package
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(tidyr)
# library(xlsx)
# library(gapindex)
# library(reshape)
# chl <- gapindex::get_connected()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Retrive previous RACE data from the sumfish package
##      Only do once and save locally.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
ai_collection <- 
  xlsx::read.xlsx(file = "data/AI/2024/AI_core_collections_2024.xlsx",
                  sheetIndex = "Sheet1")
names(x = ai_collection) <- gsub(x = names(x = ai_collection), 
                                 pattern = "\\.", 
                                 replacement = " ")

ai_collection_long <- 
  reshape::melt(data = subset(x = ai_collection,
                              select = c("SPECIES_CODE", "criteria_value",
                                         "Southern Bering Sea", 
                                         "Eastern Aleutians", 
                                         "Central Aleutians",
                                         "Western Aleutians")), 
                id.vars = c('SPECIES_CODE', 'criteria_value'), 
                measure.vars = c("Southern Bering Sea", 
                                 "Eastern Aleutians", 
                                 "Central Aleutians",
                                 "Western Aleutians"), 
                variable_name = "AREA_NAME")

## copied from G:\ALEUTIAN\AI 2024\Station Allocation\ai_2024_station_allocation.xlsx
ai_allocation <- 
  xlsx::read.xlsx(file = "data/AI/2024/ai_2024_station_allocation.xlsx",
                  sheetIndex = "Sheet1", colIndex = 1:2)

ai_allocation <- table(ai_allocation$STRATUM[
  ai_allocation$STATION_TYPE %in% c("prescribed", "new_stn")
])

ai_data <- gapindex::get_data(year_set = c(2018, 2022),
                              survey_set = "AI",
                              spp_codes = ai_collection$SPECIES_CODE, 
                              sql_channel = chl, 
                              pull_lengths = F)
ai_haul <- ai_data$haul
ai_haul$YEAR <- floor(x = ai_haul$CRUISE / 100)

ai_areas <- merge(x = subset(x = ai_data$stratum_groups,
                             subset = AREA_ID %in% c(299, 799, 3499, 5699),
                             select = c("AREA_ID", "STRATUM")),
                  y = subset(x = ai_data$subarea, 
                             select = c(AREA_ID, AREA_NAME)), 
                  by = "AREA_ID")

ai_catch <- merge(x = ai_data$catch[, c("HAULJOIN", "SPECIES_CODE", "NUMBER_FISH")],
                  y = ai_data$haul[, c("HAULJOIN", "CRUISE", "STRATUM")],
                  by = "HAULJOIN")
ai_catch$YEAR <- floor(x = ai_catch$CRUISE / 100)
ai_catch <- merge(x = ai_catch,
                  y = subset(x = ai_areas, select = -AREA_ID),
                  by = "STRATUM")

ai_catch <- merge(x = ai_catch,
                  y = ai_collection_long, 
                  by = c("SPECIES_CODE", "AREA_NAME"))

ai_catch$collect <- 
  ifelse(test = ai_catch$NUMBER_FISH > ai_catch$criteria_value, 
         yes = apply(X = ai_catch[, c("NUMBER_FISH", "value")],
                     MARGIN = 1, 
                     FUN = min),
         no = 0)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set constants for analysis
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
spp_codes <- paste(ai_collection$species_code)
n_spp <- length(x = spp_codes)

years <- sort(x = unique(x = ai_data$cruise$YEAR))
n_years <- length(x = years)

n_boots <- 100
n_haul <- sum(ai_allocation)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create result object
##   bootstrap_opt_collect will hold the optimal collection rule for each
##      species and bootstrap iteration. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
haul_total_stats <- array(dim = c(n_boots, n_years, 6), 
                          dimnames = list(NULL, paste(years), 
                                          c("Min.", "1st Qu.", "Median", 
                                            "Mean", "3rd Qu.", "Max.") ))

gt_30 <- array(dim = c(n_boots, n_years), dimnames = list(NULL, paste(years)))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Start the bootstrap
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (iter in 1:n_boots) { ## Loop over iterations -- start
  
  ## Set seed
  set.seed(iter * 23432)
  
  ## For each year, draw a vector of bootstrap indices, making sure to 
  ## sample with replacement according to the stratum allocation.
  boot_idx <- c()
  for (iyear in years) {
    for (istratum in names(x = ai_allocation) ) {
      isample <- ai_allocation[istratum]
      boot_idx <- c(boot_idx, 
                    sample(x = which(ai_haul$YEAR == iyear &
                                       ai_haul$STRATUM == istratum), 
                           size = isample, 
                           replace = TRUE))
    }
  }
  
  bootstrap_hauls <- ai_haul[boot_idx, c("YEAR", "HAULJOIN")]
  bootstrap_hauls$haulno <- 1:nrow(x = bootstrap_hauls)
  
  x <- 
    merge(x = bootstrap_hauls,
          y = ai_catch,
          all.x = TRUE,
          by = c("YEAR","HAULJOIN"))
  
  y <- aggregate(collect ~ YEAR + haulno,
                 FUN = sum,
                 data = x)
  
  haul_total_stats[iter, , ] <-
  do.call(what = rbind,
          args = lapply(X = split(x = y, f = y$YEAR), 
                        FUN = function(df) summary(df$collect)))
  
  gt_30[iter, ] <- 
    do.call(what = rbind,
            args = lapply(X = split(x = y, f = y$YEAR),
                          FUN = function(df)
                            sum(df$collect > 30) / nrow(x = df)))
  
  if(iter%%10 == 0) print(paste0("Finished with ", iter, " of ", n_boots))
} ## Loop over iterations -- end

boxplot(gt_30)
boxplot(haul_total_stats[, , "Median"])
