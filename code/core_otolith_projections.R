##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Optimize otolith collection rules
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For the AI BTS, calculate optimal otolith collection rules
##                    across species. 
##                Constraints are: total otoliths requested, distribution 
##                    and abundance of species, median haul-level collection 
##                    around 20 and avoid haul-level collections that are > 30
##                    (~ 10% of hauls are okay). 
##                Bootstrap hauls to add some variation to the analysis.
##
## Notes:         If downloading data for the first time, use the sumfish
##                   package (for Z. Oyafuso: sumfish was installed on rVersion
##                   3.6.1).
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(readxl)
library(tidyr)
library(rgdal)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Retrive AI RACE data (2010-2018) from the sumfish package
##      Only do once and save locally, takes a long time to pull data
##
##   If the AI RACE are saved locally, load it into the environment.
##   Also, load the AI allocation for 420 stations. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# library(sumfish)
# 
# sumfish::setUser()
# AI = sumfish::getRacebase(year = c(2010, 2018), survey = "AI")
# saveRDS(object = AI,
#         file = paste0("C:/Users/zack.oyafuso/Desktop/",
#                       "AIGOA_Core_Otolith/AI_RACEBASE.RDS"))

AI <- readRDS("data/AI_RACEBASE.RDS")
specimen <- AI$specimen
catch <- AI$catch

allocation <- table(read.csv("data/AIallocation420.csv")$STRATUM)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import core requests for a given year
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
collection <- readxl::read_xlsx("data/AI_core_collections.xlsx", 
                                sheet = "2022 collection")
pollock_collection <- readxl::read_xlsx("data/AI_core_collections.xlsx", 
                                        sheet = "2022 pollock collection")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    Calculate total otolith collection across species and year  
##    Note: Age data is located in the specimen slot in AI
##      with SPECIMEn_haul_TYPE == 1. These objects aren't really used, 
##      just for reference about historical collections. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
specimen <- subset(x = AI$specimen,
                   subset = SPECIMEN_SAMPLE_TYPE == 1 & 
                     SPECIES_CODE %in% collection$species_code)

total_ages <- table(specimen$SPECIES_CODE, specimen$CRUISE)     

spp_age_names <- AI$species$COMMON_NAME[match(rownames(total_ages), 
                                              AI$species$SPECIES_CODE)]
rownames(total_ages) <- spp_age_names
total_ages

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set constants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
spp_codes <- collection$species_code
n_spp <- nrow(collection)

years <- colnames(total_ages)
n_years <- length(years)

n_boots <- 10
n_haul <- 420

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Result objects
##   bootstrap_oto_main will hold the otoliths collected for a given
##      haul, year, species, and bootstrap iteration.
##   bootstrap_opt_collect will hold the optimal collection rule for each
##      species and bootstrap iteration. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
bootstrap_oto_main <- 
  array(dim = c(n_haul, n_years, n_spp, n_boots),
        dimnames = list(NULL, years, collection$species_code, NULL))

bootstrap_opt_collect <- matrix(nrow = n_boots, ncol = n_spp,
                                dimnames = list(NULL, spp_codes))

bootstrap_oto_reg <- array(dim = c(4, n_years, n_boots),
                           dimnames = list(c("CAI", "EAI", "SBS", "WAI"),
                                           years, NULL))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create a wide dataset that holds catch data at the haul level
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## Subset the catches for iyear for just the interested species
sub_catch <- subset(x = AI$catch, 
                    subset = SPECIES_CODE %in% collection$species_code,
                    select = c("CRUISE", "HAULJOIN", 
                               "SPECIES_CODE", "NUMBER_FISH"))

## Widen sub_catch to fill in zeros
catch_wide <- tidyr::spread(data = sub_catch, 
                            value = NUMBER_FISH, 
                            key = SPECIES_CODE, 
                            fill = 0)

## Get lat/lon locations of each haul in iyear
## Append to lat/lon info to catch_wide. This will be used to figure which
##   managment area the haul is in. 
## Remove hauls w/o lat/lon info.
haul_locs <- subset(x = AI$haul, 
                    subset = HAULJOIN %in% catch_wide$HAULJOIN,
                    select = c(HAULJOIN, START_LONGITUDE, 
                               START_LATITUDE, STRATUM))

catch_wide[, c("LON", "LAT", "STRATUM")] <- 
  haul_locs[match(catch_wide$HAULJOIN, haul_locs$HAULJOIN), 
            c("START_LONGITUDE", "START_LATITUDE", "STRATUM")]
catch_wide <- catch_wide[!is.na(catch_wide$LON), ]

## Based on the longitude of the station, assign management area
catch_wide$LON_E <- ifelse(test = catch_wide$LON < 0, 
                           yes = catch_wide$LON + 360,  
                           no = catch_wide$LON)
catch_wide$REG <- as.character(cut(x = catch_wide$LON_E, 
                                   breaks = c(170, 177, 183, 190, 195), 
                                   labels = c("WAI", "CAI", "EAI", "SBS")))

region_order <- catch_wide$REG[1:n_haul]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   modifications
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
collection[c(1, 7, 10), "criteria_value"] <- c(4, 3, 3)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Start the bootstrap
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
for (iter in 1:n_boots) { ## Loop over iterations -- start
  
  ## Set seed
  set.seed(iter * 23432)
  
  ## For each year, draw a vector of bootstrap indices, stratified by stratum
  boot_idx <- c()
  
  for (iyear in years) {
    for (istratum in sort(unique(catch_wide$STRATUM)) ) {
      isample <- allocation[istratum]
      boot_idx <- c(boot_idx, 
                    sample(x = which(catch_wide$CRUISE == iyear &
                                       catch_wide$STRATUM == istratum), 
                           size = isample, 
                           replace = TRUE))
    }
  }
  
  bootstrap_hauls <- catch_wide[boot_idx, ]
  
  ## For each species, optimize the haul-level otolith collection that would
  ## minimize the summed absolute difference between the total otolith request
  ## and the totals collected for each year
  for (ispp in spp_codes[]) {
    
    if (ispp %in% spp_codes[-16]) {
      
      ## Optimizer function
      fn_ <- function(collect, 
                      bootstrap_hauls_,
                      crit_val, 
                      target_val,
                      return_what = c("obj", 
                                      "total_by_haul",
                                      "total_by_region",
                                      "total_by_cruise")[1]) {
        
        ## Result output to put ototliths collected
        bootstrap_otos <- vector(length = nrow(bootstrap_hauls_))
        
        for (ihaul in 1:nrow(bootstrap_hauls_)) { ## Loop over hauls -- start  
          
          ## How many ispp_code was caught in ihaul?
          icatch <- bootstrap_hauls_[ihaul, paste0(ispp)]
          
          ## Should an otolith subsample take place?
          collect_otos <- icatch >= crit_val
          
          ## If collect_otos is TRUE, collect otoliths based on the threshold
          bootstrap_otos[ihaul] <- ifelse(test = collect_otos == TRUE,
                                          no = 0,
                                          yes = min(icatch, collect))
        } ## Loop over hauls -- end 
        
        ## Calculate sum of the absolute difference between the total cruise 
        ## collection and the target
        total_by_haul <- 
          matrix(bootstrap_otos, 
                 ncol = length(unique(bootstrap_hauls_$CRUISE)),
                 dimnames = list(NULL, 
                                 unique(bootstrap_hauls_$CRUISE)))
        
        total_by_region <- tapply(X = bootstrap_otos, 
                                  INDEX = list(bootstrap_hauls_$CRUISE,
                                               bootstrap_hauls_$REG), 
                                  FUN = sum)
        
        total_by_cruise <- tapply(X = bootstrap_otos, 
                                  INDEX = bootstrap_hauls_$CRUISE, 
                                  FUN = sum)
        
        obj <- sum(abs(total_by_cruise - target_val)^2)
        
        return(get(return_what))
      }
      
      optim_collect <- 
        optim(par = c(5), 
              fn = fn_,
              bootstrap_hauls_ = bootstrap_hauls,
              crit_val = collection$criteria_value[collection$species_code == ispp],
              target_val = collection$total_requested[collection$species_code == ispp],
              return_what = c("obj", 
                              "total_by_haul",
                              "total_by_cruise")[1], 
              method = "Brent", lower = 0, upper = 15)
      
      bootstrap_opt_collect[iter, paste0(ispp)] <- optim_collect$par
    } 
    
    if (ispp == "21740")  ## Special collection rules for pollock
      bootstrap_opt_collect[iter, paste0(ispp)]  <- 5
    
  } ## Loop over species -- end
  
  if(iter%%2 == 0) print(paste0("Finished with ", iter, " of ", n_boots))
} ## Loop over iterations -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Recalculate bootstrap species totals using the median collection rule, 
##   round up to the next integer
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
for (iter in 1:n_boots) { ## Loop over iterations -- start
  
  ## Set seed
  set.seed(iter * 23432)
  
  ## For each year, draw a vector of bootstrap indices, stratified by stratum
  boot_idx <- c()
  
  for (iyear in years) {
    for (istratum in sort(unique(catch_wide$STRATUM)) ) {
      isample <- allocation[istratum]
      boot_idx <- c(boot_idx, 
                    sample(x = which(catch_wide$CRUISE == iyear &
                                       catch_wide$STRATUM == istratum), 
                           size = isample, 
                           replace = TRUE))
    }
  }
  
  bootstrap_hauls <- catch_wide[boot_idx, ]
  
  for (ispp in paste0(spp_codes)[]) {
    if (ispp %in% spp_codes[-16]) {
      ## Optimizer function
      rounded_collect <- ceiling(median(bootstrap_opt_collect[, paste0(ispp)]))
      
      bootstrap_oto_main[, , paste0(ispp), iter] <- 
        fn_(collect = rounded_collect,
            bootstrap_hauls_ = bootstrap_hauls,
            crit_val = collection$criteria_value[collection$species_code == ispp],
            target_val = collection$total_requested[collection$species_code == ispp],
            return_what = "total_by_haul")
      
    } 
    
    if (ispp == "21740") { ## Special collection rules for pollock
      bootstrap_otos <- vector(length = nrow(bootstrap_hauls)) 
      
      for (ihaul in 1:nrow(bootstrap_hauls)) {
        iarea <- bootstrap_hauls$REG[ihaul]
        
        icatch <- bootstrap_hauls[ihaul, paste(ispp)]
        icollect <- cut(x = icatch,
                        breaks = c(-1, subset(x = pollock_collection,
                                              area == iarea)$max_threshold),
                        labels = subset(x = pollock_collection,
                                        area == iarea)$collect_n)
        bootstrap_otos[ihaul] <- as.integer(paste(icollect))
      }
      
      bootstrap_oto_main[, , paste0(ispp), iter] <- 
        matrix(bootstrap_otos, 
               ncol = length(unique(bootstrap_hauls$CRUISE)),
               dimnames = list(NULL, 
                               unique(bootstrap_hauls$CRUISE) ))
    }
  }
  
  bootstrap_haul_totals <- 
    apply(X = bootstrap_oto_main[, , , iter],
          MARGIN = c(1, 2),
          FUN = function(x) sum(x))
  
  region_order <- bootstrap_hauls$REG[1:n_haul]
  
  bootstrap_oto_reg[, , iter] <- 
    apply(X = bootstrap_haul_totals, 
          MARGIN = 2, 
          FUN = function(x) tapply(X = x, 
                                   INDEX = region_order, 
                                   FUN = function(y) sum(y > 30) ))
  
}

over_30_TF <- apply(X = bootstrap_oto_main,
                    MARGIN = c(1, 2, 4),
                    FUN = function(x) sum(x) > 30 ) 
over_30 <- apply(X = over_30_TF, 
                 MARGIN = 3:2, 
                 FUN = function(x) round(x = sum(x) / n_haul, digits = 2))
cruise_totals <- apply(X = bootstrap_oto_main, 
                       MARGIN = c(2:4), 
                       FUN = sum)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plots
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(4, 4), mar = c(2, 4, 2, 3))
for (ispp in paste0(spp_codes)[]) {
  boxplot(t(cruise_totals[, ispp, ]),  
          main = collection$common_name[collection$species_code == ispp], 
          las = 1)
  abline(h = collection$total_requested[collection$species_code == ispp])
}


for (ispp in paste0(spp_codes)[]) {
  boxplot(bootstrap_opt_collect[, ispp],
          las = 1,
          main = collection$common_name[collection$species_code == ispp])
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Update collection rules 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
collection$collection_rule_final <- 
  c(2, 8, 8, 5, 10, 10, 3, 5, 2, 3, 5, 10, 10, 5, 10, 5)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create table for deck poster...
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

###############################################################################
####   Save Results
###############################################################################
# save(list = c("bootstrap_oto_main", "bootstrap_total_work", 
#               "bootstrap_over_thirty", "total_ages"), 
#      file = "results/bootstrap_results.RData")
