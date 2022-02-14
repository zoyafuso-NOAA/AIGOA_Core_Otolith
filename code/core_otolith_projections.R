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

###############################################################################
####   Bootstrap 420 hauls for each year
############################################################################### 
years <- colnames(total_ages)
n_boots <- 100
n_sample <- 420

bootstrap_oto_main <- 
  array(dim = c(length(years),
                nrow(collection),
                n_boots),
        dimnames = list(years, 
                        collection$species_code, 
                        NULL))

bootstrap_total_work <- 
  array(dim = c(length(years), 5, n_boots),
        dimnames = list(years,
                        paste0(c(2.5, 25, 50, 75, 97.5), "%"),
                        NULL))

for (iyear in years) { ## Loop over years -- start
  
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
                      select = c(HAULJOIN, START_LONGITUDE, START_LATITUDE))
  
  catch_wide[, c("LON", "LAT")] <- 
    haul_locs[match(catch_wide$HAULJOIN, haul_locs$HAULJOIN), 
              c("START_LONGITUDE", "START_LATITUDE")]
  catch_wide <- catch_wide[!is.na(catch_wide$LON), ]
  
  ## Based on the longitude of the station, assign management area
  catch_wide$LON_E <- ifelse(test = catch_wide$LON < 0, 
                             yes = catch_wide$LON + 360,  
                             no = catch_wide$LON)
  catch_wide$REG <- as.character(cut(x = catch_wide$LON_E, 
                                     breaks = c(170, 177, 183, 190, 195), 
                                     labels = c("WAI", "CAI", "EAI", "SBS")))
  
  ## Start the bootstrap
  for (iter in 1:n_boots) { ## Loop over iterations -- start
    
    ## Set seed
    set.seed(iter * 23432)
    
    ## Draw the bootstrap sample
    bootstrap_sample <- sample(x = 1:nrow(catch_wide),
                               size = n_sample, 
                               replace = TRUE)
    bootstrap_hauls <- catch_wide[bootstrap_sample, ]
    
    ## Result output to put ototliths collected
    bootstrap_otos <- matrix(nrow = n_sample, 
                             ncol = nrow(collection),
                             dimnames = list(NULL, collection$species_code))
    
    for (ihaul in 1:nrow(bootstrap_hauls)) { ## Loop over hauls -- start   
      for (ispp in 1:nrow(collection)) { ## Loop over species -- start
        
        ## what management area is ihaul in and what is the species code
        iarea <- bootstrap_hauls$REG[ihaul]
        ispp_code <- paste0(collection$species_code[ispp])
        
        ## How many ispp_code was caught in ihaul?
        icatch <- bootstrap_hauls[ihaul, ispp_code]
        
        ## Should an otolith subsample take place?
        collect_otos <- icatch >= collection$criteria_value[ispp]
        
        ## If collect_otos is TRUE, collect otoliths based on the
        bootstrap_otos[ihaul, ispp_code] <-
          ifelse(test = collect_otos == TRUE,
                 no = 0,
                 yes = min(icatch, as.integer(collection[ispp, iarea]) ))
      } ## Loop over species -- end
      
      ## Special collection rules for pollock
      ispp_code <- "21740"
      iarea <- bootstrap_hauls$REG[ihaul]
      icatch <- bootstrap_hauls[ihaul, ispp_code]
      icollect <- cut(x = icatch, 
                      breaks = c(-1, subset(x = pollock_collection,
                                            area == iarea)$max_threshold), 
                      labels = subset(x = pollock_collection,
                                      area == iarea)$collect_n)
      bootstrap_otos[ihaul, ispp_code] <- as.integer(paste(icollect))
      
    } ## Loop over hauls -- end 
    
    ## Record results
    bootstrap_oto_main[iyear, , iter] <- colSums(bootstrap_otos)
    bootstrap_total_work[iyear, , iter] <- 
      quantile(x = rowSums(bootstrap_otos),
               probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    
    if(iter%%10 == 0) print(paste0("Finished with ", iter, " of ", n_boots, 
                                   ", Year ", iyear))
  } ## Loop over iterations -- end
  
} ## Loop over years -- end

###############################################################################
####   Save Results
###############################################################################  
save(list = c("bootstrap_oto_main", "bootstrap_total_work", "total_ages"), 
     file = "results/bootstrap_results.RData")
