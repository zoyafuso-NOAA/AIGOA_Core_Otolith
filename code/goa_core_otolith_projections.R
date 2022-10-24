##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Optimize otolith collection rules
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For the GOA BTS, calculate optimal otolith collection rules
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
##   Retrive GOA RACE data (2017-2021) from the sumfish package
##      Only do once and save locally, takes a long time to pull data
##
##   If the GOA RACE are saved locally, load it into the environment.
##   Also, load the AI allocation for 420 stations. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# library(sumfish)
# 
# sumfish::setUser()
# GOA = sumfish::getRacebase(year = c(2017, 2021), survey = "GOA")
# saveRDS(object = GOA, file = "data/GOA_RACEBASE.RDS")

GOA <- readRDS("data/GOA_RACEBASE.RDS")
specimen <- GOA$specimen
catch <- GOA$catch

allocation <- table(read.csv("data/GOAallocation2021.csv")$stratum)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import core requests for a given year
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
collection <- read.csv("data/GOA_core_collections.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    Calculate total otolith collection across species and year  
##    Note: Age data is located in the specimen slot in AI
##      with SPECIMEn_haul_TYPE == 1. These objects aren't really used, 
##      just for reference about historical collections. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
specimen <- subset(x = GOA$specimen,
                   subset = SPECIMEN_SAMPLE_TYPE == 1)

total_ages <- table(specimen$SPECIES_CODE, specimen$CRUISE)

spp_age_names <- GOA$species$COMMON_NAME[match(rownames(total_ages),
                                               GOA$species$SPECIES_CODE)]
rownames(total_ages) <- spp_age_names
total_ages

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set constants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
spp_codes <- collection$species_code
n_spp <- length(spp_codes)

years <- colnames(total_ages)[-1]
n_years <- length(years)

n_boots <- 30
n_haul <- 540

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

bootstrap_oto_cruise <- 
  array(dim = c(n_years, n_spp, n_boots),
        dimnames = list(years, collection$species_code, NULL))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create a wide dataset that holds catch data at the haul level
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## Subset the catches for iyear for just the interested species
sub_catch <- subset(x = GOA$catch, 
                    subset = SPECIES_CODE %in% spp_codes & CRUISE %in% years,
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
haul_locs <- subset(x = GOA$haul, 
                    subset = HAULJOIN %in% catch_wide$HAULJOIN,
                    select = c(HAULJOIN, START_LONGITUDE, 
                               START_LATITUDE, STRATUM))

catch_wide[, c("LON", "LAT", "STRATUM")] <- 
  haul_locs[match(catch_wide$HAULJOIN, haul_locs$HAULJOIN), 
            c("START_LONGITUDE", "START_LATITUDE", "STRATUM")]
catch_wide <- catch_wide[!is.na(catch_wide$LON), ]

catch_wide$DEPTH <- GOA$haul$BOTTOM_DEPTH[match(catch_wide$HAULJOIN, GOA$haul$HAULJOIN)]

## Based on the longitude of the station, assign management area
# catch_wide$LON_E <- ifelse(test = catch_wide$LON < 0, 
#                            yes = catch_wide$LON + 360,  
#                            no = catch_wide$LON)
# catch_wide$REG <- as.character(cut(x = catch_wide$LON_E, 
#                                    breaks = c(170, 177, 183, 190, 195), 
#                                    labels = c("WAI", "CAI", "EAI", "SBS")))

# region_order <- catch_wide$REG[1:n_haul]

lengths <- subset(GOA$raw_length, SPECIES_CODE %in% collection$species_code[collection$sample_type == "stratified"] )

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   modifications
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

collection <- read.csv("data/GOA_core_collections.csv")

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
  for (ispp in spp_codes[5]) {
    
    ## Optimizer function
    fn_ <- function(collect,
                    bootstrap_hauls_,
                    crit_val,
                    target_val,
                    max_depth,
                    return_what = c("obj",
                                    "total_by_cruise")[1]) {
      
      ## Result output to put ototliths collected
      bootstrap_otos <- vector(length = nrow(bootstrap_hauls_))
      
      for (ihaul in 1:nrow(bootstrap_hauls_)) { ## Loop over hauls -- start  
        
        ## How many ispp_code was caught in ihaul?
        icatch <- bootstrap_hauls_[ihaul, paste0(ispp)]
        
        ## Should an otolith subsample take place?
        collect_otos <- 
          icatch >= crit_val & bootstrap_hauls_[ihaul, "DEPTH"] < max_depth
        
        ## If collect_otos is TRUE, collect otoliths based on the threshold
        bootstrap_otos[ihaul] <- ifelse(test = collect_otos == TRUE,
                                        no = 0,
                                        yes = min(icatch, collect))
      } ## Loop over hauls -- end 
      
      ## Calculate sum of the absolute difference between the total cruise 
      ## collection and the target
      total_by_cruise <- tapply(X = bootstrap_otos, 
                                INDEX = bootstrap_hauls_$CRUISE, 
                                FUN = sum)
      
      obj <- sum(abs(total_by_cruise - target_val)^2)
      
      return(get(return_what))
    }
    
    if (collection$sample_type[collection$species_code == ispp] == "random") {
      optim_collect <-
        optim(par = collection$init_collect[collection$species_code == ispp],
              fn = fn_,
              bootstrap_hauls_ = bootstrap_hauls,
              crit_val = collection$min_val[collection$species_code == ispp],
              target_val = collection$target[collection$species_code == ispp],
              max_depth = collection$max_depth[collection$species_code == ispp],
              return_what = c("obj",
                              "total_by_cruise")[1],
              method = "Brent", lower = 0, upper = 30)
      
      bootstrap_opt_collect[iter, paste0(ispp)] <- optim_collect$par
      
      bootstrap_oto_cruise[ , paste0(ispp), iter] <- 
        fn_(collect = optim_collect$par,
            bootstrap_hauls_ = bootstrap_hauls,
            crit_val = collection$min_val[collection$species_code == ispp],
            target_val = collection$target[collection$species_code == ispp],
            max_depth = collection$max_depth[collection$species_code == ispp],
            return_what = c("obj",
                            "total_by_cruise")[2])
    }
    
    
    if (collection$sample_type[collection$species_code == ispp] == "stratified") {
      
      bootstrap_opt_collect[iter, paste0(ispp)] <- 
        collection$init_collect[collection$species_code == ispp]
      
      fn_ <- function(collect,
                      bootstrap_hauls_) {
        
        ## Result output to put ototliths collected
        bootstrap_otos <- vector(length = nrow(bootstrap_hauls_))
        
        for (ihaul in 1:nrow(bootstrap_hauls_)) { ## Loop over hauls -- start  
          
          ## How many ispp_code was caught in ihaul?
          icatch <- bootstrap_hauls_[ihaul, paste0(ispp)]
          
          if (icatch == 0) bootstrap_otos[ihaul] <- 0
          if (icatch > 0) {
            table_sex_length <- 
              with(subset(x = lengths, 
                          subset = CRUISE == bootstrap_hauls_$CRUISE[ihaul] &
                            HAULJOIN == bootstrap_hauls_$HAULJOIN[ihaul] & 
                            SPECIES_CODE == ispp),
                   table(SEX, LENGTH) )
            
            collect_otos <- apply(X = table_sex_length, 
                                  MARGIN = 1:2, 
                                  FUN = function(x) ifelse(test = x > 0, 
                                                           yes = min(collect, x), 
                                                           no = 0))
            
            ## If collect_otos is TRUE, collect otoliths based on the threshold
            bootstrap_otos[ihaul] <- sum(collect_otos)
          }
        } ## Loop over hauls -- end 
        
        total_by_cruise <- tapply(X = bootstrap_otos, 
                                  INDEX = bootstrap_hauls_$CRUISE, 
                                  FUN = sum)
        return(total_by_cruise)
      }
      
      bootstrap_oto_cruise[ , paste0(ispp), iter] <- 
        fn_(collect = collection$init_collect[collection$species_code == ispp],
            bootstrap_hauls_ = bootstrap_hauls)
      
    }
    
    
    
    
  } ## Loop over species -- end
  
  if(iter%%10 == 0) print(paste0("Finished with ", iter, " of ", n_boots))
} ## Loop over iterations -- end

par(mar = c(2, 6, 1, 1), mfrow = c(1, 1))
boxplot(bootstrap_opt_collect, horizontal = T, las = 1, ylim = c(0, 30))

par(mar = c(2, 3, 3, 1), mfrow = c(3, 4))
for (ispp in spp_codes) {
  boxplot( t(bootstrap_oto_cruise[ , paste0(ispp), ]), main = ispp)
  abline(h = collection$target[collection$species_code == ispp]) 
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calcculate which species are the most frequent
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
names(sort(apply(X = catch_wide[, paste(spp_codes)], MARGIN = 2, FUN = function(x) sum(x > 0) / length(x)), decreasing = TRUE))

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
            crit_val = collection$min_val[collection$species_code == ispp],
            target_val = collection$target[collection$species_code == ispp],
            max_depth = collection$max_depth[collection$species_code == ispp],
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
