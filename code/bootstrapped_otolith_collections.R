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
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyr)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Retrive previous RACE data from the sumfish package
##      Only do once and save locally.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# library(sumfish)
# sumfish::setUser()
# data <- sumfish::getRacebase(year = c(2019, 2021), survey = "GOA")
# saveRDS(object = data, file = "data/GOA/2023/GOA_RACEBASE.RDS")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set constants about the data inputs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
region <- c("AI", "GOA")[2]
year <- 2023
database_filename <- "GOA_RACEBASE.RDS"
allocation_filename <- "GOA2023_Station_allocation_520_EW.csv"
collection_requests_filename <- "GOA_core_collections.csv"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import data
##   data: historical specimen, catch, and haul data from past 2 survey years
##   allocation: allocations of stations across stratum
##   collection: collection rules and requested targets for each species
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- 
  readRDS(file = paste0("data/", region, "/", year, "/", database_filename))
specimen <- data$specimen
catch <- data$catch
lengths <- data$length

allocation_table <- 
  read.csv(file = paste0("data/", region, "/", year, "/", allocation_filename))
allocation <- table(allocation_table$stratum)

collection <- read.csv(file = paste0("data/", region, "/", year, 
                                     "/", collection_requests_filename))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Optimizer functions:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("code/functions/core_otolith_fns.R")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    Calculate histocial total otolith collection across species and year  
##    Note: Age data is located in the specimen slot in AI
##      with SPECIMEN_SAMPLE_TYPE == 1. These objects aren't really used, 
##      just for reference about historical collections. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
specimen <- subset(x = data$specimen,
                   subset = SPECIMEN_SAMPLE_TYPE == 1)

total_ages <- table(specimen$SPECIES_CODE, specimen$CRUISE)

spp_age_names <- data$species$COMMON_NAME[match(rownames(total_ages),
                                                data$species$SPECIES_CODE)]
rownames(total_ages) <- spp_age_names
total_ages

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set constants for analysis
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
spp_codes <- paste(collection$species_code)
n_spp <- length(x = spp_codes)

years <- colnames(x = total_ages)
n_years <- length(x = years)

n_boots <- 100
n_haul <- 520

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create result object
##   bootstrap_opt_collect will hold the optimal collection rule for each
##      species and bootstrap iteration. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
bootstrap_opt_collect <- matrix(nrow = n_boots, 
                                ncol = n_spp,
                                dimnames = list(NULL, spp_codes))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create a wide dataset that holds catch data at the haul level
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Subset the catches for iyear for just the interested species
sub_catch <- subset(x = data$catch, 
                    subset = SPECIES_CODE %in% spp_codes & CRUISE %in% years,
                    select = c("CRUISE", "HAULJOIN", 
                               "SPECIES_CODE", "NUMBER_FISH"))

## Widen sub_catch to fill in zeros
catch_wide <- tidyr::spread(data = sub_catch, 
                            value = NUMBER_FISH, 
                            key = SPECIES_CODE, 
                            fill = 0)

## Get stratum and depth information for each haul.
haul_locs <- subset(x = data$haul, 
                    subset = HAULJOIN %in% catch_wide$HAULJOIN,
                    select = c(HAULJOIN, STRATUM, BOTTOM_DEPTH))

catch_wide[, c("STRATUM", "DEPTH")] <- 
  haul_locs[match(catch_wide$HAULJOIN, haul_locs$HAULJOIN), 
            c("STRATUM", "BOTTOM_DEPTH")]

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
  for (ispp in spp_codes) {
    
    ## Set up temporary control variables
    sample_type <- collection$sample_type[collection$species_code == ispp]
    init_collect <- collection$init_collect[collection$species_code == ispp]
    min_val <- collection$min_val[collection$species_code == ispp]
    target <- collection$target[collection$species_code == ispp]
    max_depth <- collection$max_depth[collection$species_code == ispp]
    
    if (sample_type == "random") {
      
      ## Optimize collection quantity
      optim_collect <- optim(par = init_collect,
                             fn = otolith_opt,
                             bootstrap_hauls_ = bootstrap_hauls,
                             species_code = ispp,
                             crit_val = min_val,
                             target_val = target,
                             max_depth = max_depth,
                             return_what = c("obj",
                                             "total_by_cruise")[1],
                             method = "Brent", lower = 0, upper = 100)
      
      ## Record optimal colletion quantity
      bootstrap_opt_collect[iter, ispp] <- optim_collect$par
    }
    
    if (sample_type == "stratified") {
      ## Record "optimal" colletion as just the requested quantity
      bootstrap_opt_collect[iter, ispp] <- init_collect
    }
    
  } ## Loop over species -- end
  
  if(iter%%10 == 0) print(paste0("Finished with ", iter, " of ", n_boots))
} ## Loop over iterations -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot the distribution of optimal collection rules for each species
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mar = c(1, 4, 2, 1), mfrow = c(4, 4))
for (ispp in paste(spp_codes)[collection$sample_type == "random"]) {
  boxplot(bootstrap_opt_collect[, ispp],
          las = 1, 
          main = collection$species[collection$species_code == ispp])
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Recalculate bootstrap species totals using the median collection rule, 
##   round up to the next integer (or the collection rule that you end up 
##   negotiate with the stock assessor)
##
##   bootstrap_oto_main will hold the otoliths collected for a given
##      haul, year, species, and bootstrap iteration.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
bootstrap_oto_main <- 
  array(dim = c(n_haul, n_years, n_spp, n_boots),
        dimnames = list(NULL, years, spp_codes, NULL))
bootstrap_haul_totals <- array(dim = c(n_haul, n_years, n_boots),
                               dimnames = list(NULL, years, NULL))

collection$final_collect <- 
  c(10, # Atka mackerel, Up to 10 random/haul
    5,  # walleye pollock, Up to 5 random/haul 
    4,  # Pacific cod, Up to 4 random/haul
    4,  # northern rock sole, Up to 4 random/haul
    3,  # southern rock sole, If >= 5 then collect 3 random/haul, 0 otherwise
    3,  # rex sole, if >= 10 then collect 2 random/haul, 0 otherwise
    3,  # Dover sole, If >= 3 then collect 3 random/haul, 0 otherwise
    2,  # flathead sole, If >= 5 then collect 2 random/haul, 0 otherwise
    2,  # arrowtooth flounder, If >= 10 then collect 2 random/haul, 0 otherwise
    4,  # rougheye rockfish, Up to 5 random/haul
    10, # blackspotted rockfish, Up to 10 random/haul
    5,  # Pacific ocean perch, If >= 5 then collect 5 random/haul, 0 otherwise
    15, # northern rockfish, Up to 20 random/haul
    20, # shortaker rockfish, Up to 20 random/haul
    10, # dusky rockfish, Up to 10 random/haul
    1,  # silvergray rockfish, 1 per cm/sex/haul
    2,  # harlequin rockfish, 2 per cm/sex/haul
    2)  # sablefish, Up to 2 random/haul for hauls shallower than 200 m

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
  
  for (ispp in spp_codes) {
    
    ## Set up temporary control variables
    sample_type <- collection$sample_type[collection$species_code == ispp]
    final_collect <- collection$final_collect[collection$species_code == ispp]
    min_val <- collection$min_val[collection$species_code == ispp]
    target <- collection$target[collection$species_code == ispp]
    max_depth <- collection$max_depth[collection$species_code == ispp]
    
    if (sample_type == "random") {
      bootstrap_oto_main[, , ispp, iter] <-
        otolith_opt(collect = final_collect,
                    species_code = ispp,
                    bootstrap_hauls_ = bootstrap_hauls,
                    crit_val = min_val,
                    target_val = target,
                    max_depth = max_depth,
                    return_what = "total_by_haul" )
    }
    
    if (sample_type == "stratified") {
      bootstrap_oto_main[, , ispp, iter] <-
        stratified_samp_1cm(collect = final_collect,
                            species_code = ispp,
                            bootstrap_hauls_ = bootstrap_hauls,
                            lengths_ = lengths)$total_by_haul
    }
    
  } 
  
  bootstrap_haul_totals[, , iter] <- 
    apply(X = bootstrap_oto_main[, , , iter],
          MARGIN = c(1, 2),
          FUN = function(x) sum(x))
  
  if(iter%%10 == 0) print(paste0("Finished with ", iter, " of ", n_boots))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate the frequency/proportion of hauls where > 30 otoliths are 
##        collected. We want to see how often we hit this upper soft limit. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
over_30_TF <- apply(X = bootstrap_oto_main,
                    MARGIN = c(1, 2, 4),
                    FUN = function(x) sum(x) > 30 )
over_30 <- apply(X = over_30_TF,
                 MARGIN = 3:2,
                 FUN = function(x) round(x = sum(x) / n_haul, digits = 2))

par(mfrow = c(1, 1))
boxplot(over_30)

summary(apply(bootstrap_oto_main,
              MARGIN = c(1, 2, 4),
              FUN = sum))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Tabulate the total otoliths collected by species
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bootstrap_oto_cruise <- apply(bootstrap_oto_main,
                              MARGIN = c(2, 3, 4),
                              FUN = sum)

par(mar = c(3, 3, 2, 1), mfrow = c(4, 5))
for (ispp in paste(spp_codes)) {
  boxplot( t(bootstrap_oto_cruise[ , ispp, ]), 
           main = collection$species[collection$species_code == ispp], 
           cex.main = 0.75, 
           cex.axis = 0.75,
           names = c(2019, 2021),
           las = 1)
  abline(h = collection$target[collection$species_code == ispp]) 
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save objects
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c("bootstrap_opt_collect", "bootstrap_oto_cruise", 
  "bootstrap_oto_main", "over_30")
