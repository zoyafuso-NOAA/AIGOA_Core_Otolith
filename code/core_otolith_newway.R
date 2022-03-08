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
library(rgdal)

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
pollock_collection <- readxl::read_xlsx("data/AI_core_collections.xlsx", 
                                        sheet = "2022 pollock collection")

##################################################
####    Calculate total otolith collection across species and year  
####    note: Age data is located in the specimen slot in AI
####      with SPECIMEN_SAMPLE_TYPE == 1
##################################################
specimen <- subset(x = AI$specimen,
                   subset = SPECIMEN_SAMPLE_TYPE == 1 & 
                     SPECIES_CODE %in% collection$species_code)

total_ages <- table(specimen$SPECIES_CODE, specimen$CRUISE)     

spp_age_names <- AI$species$COMMON_NAME[match(rownames(total_ages), 
                                              AI$species$SPECIES_CODE)]
rownames(total_ages) <- spp_age_names

iyear = 201801

## Subset the catches for iyear for just the interested species
sub_catch <- subset(x = AI$catch, 
                    subset = SPECIES_CODE %in% collection$species_code & 
                      CRUISE == iyear,
                    select = c("HAULJOIN","SPECIES_CODE", "NUMBER_FISH"))

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
                    select = c(HAULJOIN, START_LONGITUDE, START_LATITUDE, STRATUM))

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

bootstrap_collections <- matrix(nrow = 100, ncol = 14)

for (iter in 1:100) {
  ## Set seed
  set.seed(iter * 23432)
  
  ## Draw the bootstrap sample, stratified by stratum
  boot_idx <- c()
  
  for (istratum in sort(unique(catch_wide$STRATUM)) ) {
    isample <- allocation[istratum]
    boot_idx <- c(boot_idx, sample(x = which(catch_wide$STRATUM == istratum), 
                                   size = isample, 
                                   replace = TRUE))
  }
  
  bootstrap_hauls <- catch_wide[boot_idx, ]
  
  catch <- bootstrap_hauls[, paste(collection$species_code[1:14])]
  # target <- collection$total_requested[1:14]
  # threshold <- collection$criteria_value[1:14]

  
  fn_ <- function(collect, target, threshold, 
                  return_what = c("obj", "total_by_spp")[1]) {
    
    collection_by_haul <- matrix(data = 0,
                                 nrow = nrow(catch), 
                                 ncol = length(collect))
    
    abs_diff <- 0
    collection <- rep(0, length(collect))
    
    for (ispp in 1:length(collect)){
      collect_idx <- catch[, ispp] >= threshold[ispp]
      
      if(!all(collect_idx == F)) {
        collection_by_haul[collect_idx, ispp] <- 
          sapply(X = catch[collect_idx, ispp], 
                 FUN = function(x) min(x, 
                                       collect[ispp]) )
        collection[ispp] <- sum(collection_by_haul[collect_idx, ispp])
        
        abs_diff <- abs_diff + abs(collection[ispp] - target[ispp])
      }
    }
    
    if(return_what == "obj") return(abs_diff)
    if(return_what == "total_by_spp") return(collection_by_haul)
  }
  
  threshold_ <- collection$criteria_value[1:14]
  
  test = optim(par = rep(1, 14), 
               fn = fn_, method = "BFGS",
               target = collection$total_requested[1:14],
               threshold = threshold_,
               return_what = "obj")
  
  best_collect <- ceiling(test$par)
  best_collect[best_collect > 15] <- 15
  
  temp_collection <- fn_(collect = best_collect,
                         target = collection$total_requested[1:14],
                         threshold = threshold_, 
                         return_what = "total_by_spp")
  hist(rowSums(temp_collection))
  sum(rowSums(temp_collection) > 30)
  
  colSums(temp_collection)
  
  bootstrap_collections[iter, ] <- test$par
  
  print(paste0("Iteration: ", iter))
  print(test$par)
}

for (ispp in 1:14) {
  hist(bootstrap_collections[, ispp], freq = F, 
       las = 1, 
       xlim = c(0, 1.05 * max(bootstrap_collections[, ispp], 
                       collection$SBS[ispp])),
       main = collection$common_name[ispp], 
       xlab = "Collection Rule (random per haul)")
  abline(v = collection$SBS[ispp], lwd = 3, col = "red")
}

ceiling(apply(X = bootstrap_collections, MARGIN = 2, FUN = median))
