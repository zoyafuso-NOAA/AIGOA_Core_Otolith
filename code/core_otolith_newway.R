###############################################################################
## Project:       Synthesize Alaska Bottom Trawl Arctic otter trawl survey data 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov), modfied from 
##                Lewis Barnett (lewis.barnett@noaa.gov)
##
## Notes:         If downloading data for the first time, use the sumfish
##                   package (for Z. Oyafuso: sumfish was installed on rVersion
##                   3.6.1)
###############################################################################
rm(list = ls())

###############################################################################
####   Import packages
###############################################################################  
library(readxl)
library(tidyr)

##################################################
####   Retrive AI RACE data (2010-2018) from the sumfish package
####      Only do once and save locally, takes a long time to pull data
####
####   If the AI RACE are saved locally, load it into the environment
##################################################
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

###############################################################################
####   Import core requests for a given year
###############################################################################  
collection <- readxl::read_xlsx("data/AI_core_collections.xlsx", 
                                sheet = "2022 collection")
collection <- subset(x = collection, subset = common_name != "walleye pollock")
pollock_collection <- readxl::read_xlsx("data/AI_core_collections.xlsx", 
                                        sheet = "2022 pollock collection")

##################################################
####   Constants
##################################################  
years <- sort(unique(catch$CRUISE))
n_years <- length(years)

n_boot <- 420
n_iters <- 1000

species <- collection$common_name
species_codes <- collection$species_code
n_spp <- nrow(collection)

##################################################
####    Calculate total otolith collection across species and year  
####    note: Age data is located in the specimen slot in AI
####      with SPECIMEN_SAMPLE_TYPE == 1
##################################################
specimen <- subset(x = AI$specimen,
                   subset = SPECIMEN_SAMPLE_TYPE == 1 & 
                     SPECIES_CODE %in% species_codes)

total_ages <- table(specimen$SPECIES_CODE, specimen$CRUISE)     

spp_age_names <- AI$species$COMMON_NAME[match(rownames(total_ages), 
                                              AI$species$SPECIES_CODE)]
rownames(total_ages) <- spp_age_names

##################################################
####   Widen the catch dataset and append management area to each haul
##################################################  
## Subset the catches for just the interested species
sub_catch <- subset(x = AI$catch, 
                    subset = SPECIES_CODE %in% species_codes,
                    select = c("HAULJOIN","SPECIES_CODE", 
                               "NUMBER_FISH", "CRUISE"))

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

##################################################
####   Empty Result Objects
##################################################  
boot_spp <- array(dim = c(n_iters, n_spp, n_years), 
                  dimnames = list(NULL, species, years))
boot_haul <- array(dim = c(n_iters, n_boot, n_years), 
                   dimnames = list(NULL, NULL, years))
boot_collect <- array(dim = c(n_iters, n_spp),
                      dimnames = list(NULL, species))

##################################################
####   Conduct bootstrap
##################################################  
for (iter in 1:n_iters) {
  
  ## Set seed
  set.seed(iter * 23432)
  
  ## Draw the bootstrap sample, stratified by stratum and year
  bootstrap_hauls <- list()
  
  for (iyear in 1:n_years) {
    boot_idx <- c()
    for (istratum in sort(unique(catch_wide$STRATUM)) ) {
      isample <- allocation[istratum]
      boot_idx <- c(boot_idx, 
                    sample(x = which(catch_wide$STRATUM == istratum
                                     & catch_wide$CRUISE == years[iyear]), 
                           size = isample, 
                           replace = TRUE))
    }
    bootstrap_hauls[[paste0(years[iyear])]] <- catch_wide[boot_idx, ]
  }
  
  ## Optim function 
  fn_ <- function(collect, 
                  bootstrap_hauls_,
                  collection_input, 
                  return_what = c("obj", "total_by_spp", "total_by_haul")[1]) {
    
    target <- collection_input$total_requested
    threshold <- collection_input$criteria_value
    spp_codes <- collection_input$species_code
    names(collect) <- names(target) <- names(threshold) <- spp_codes
    
    collection_by_haul <- array(data = 0,
                                dim = c(nrow(bootstrap_hauls_[[1]]), 
                                        length(collect), 
                                        length(bootstrap_hauls_)),
                                dimnames = list(NULL, 
                                                spp_codes, 
                                                names(bootstrap_hauls_)))
    
    for (ispp in paste0(spp_codes)) {
      for (iyear in names(bootstrap_hauls_)) {
        
        collect_idx <- rep(TRUE, dim(collection_by_haul)[1])
        if (ispp == "10262") {
          region = bootstrap_hauls_[[iyear]]$REG 
          collect_idx[region %in% c("WAI", "CAI")] <- FALSE
        }
        
        collection_by_haul[collect_idx, ispp, iyear] <-
          sapply(X = bootstrap_hauls_[[iyear]][collect_idx, ispp], 
                 FUN = function(x) ifelse(test = x >= threshold[ispp], 
                                          yes = min(x, collect[ispp]), 
                                          no = 0))
      }
    }
    
    collection_by_year <- apply(X = collection_by_haul, 
                                MARGIN = 2:3, 
                                FUN = sum)
    
    abs_diff <- sum(abs(sweep(x = collection_by_year, 
                              MARGIN = 1, 
                              STATS = target,
                              FUN = "-")))
    
    if(return_what == "obj") return(abs_diff)
    if(return_what == "total_by_spp") return(collection_by_year)
    if(return_what == "total_by_haul") return(apply(X = collection_by_haul, 
                                                    MARGIN = c(1, 3), 
                                                    FUN = sum))
  }
  
  ## Run Optimization
  test = optim(par = rep(1, n_spp), 
               fn = fn_, method = "BFGS", 
               bootstrap_hauls_ = bootstrap_hauls,
               collection_input = collection,
               return_what = "obj")
  
  ## Save optimal results
  boot_spp[iter, , ] <- fn_(collect = test$par,
                            bootstrap_hauls_ = bootstrap_hauls,
                            collection_input = collection,
                            return_what = "total_by_spp")
  
  boot_haul[iter, , ] <- fn_(collect = test$par,
                             bootstrap_hauls_ = bootstrap_hauls,
                             collection_input = collection,
                             return_what = "total_by_haul")
  
  boot_collect[iter, ] <- test$par
  
  ## Print results
  print(paste0("Iteration: ", iter))
  print(cbind(data.frame(name = collection$common_name, 
                         collect = round(test$par, 1)),
              boot_spp[iter, ,]) )
  
  ## Periodially save results
  if (iter%%10 == 0) save(list =  c("years", "n_years", "n_boot", "n_iters", 
                                    "species", "species_codes", "n_spp",
                                    "total_ages", 
                                    "boot_spp", "boot_haul", "boot_collect"), 
                          file = "results/bootstrap_results.RData")
  
  
}
