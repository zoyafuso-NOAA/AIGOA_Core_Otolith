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

AI <- readRDS(paste0("C:/Users/zack.oyafuso/Desktop/",
                     "AIGOA_Core_Otolith/AI_RACEBASE.RDS"))
specimen <- AI$specimen
catch <- AI$catch

##################################################
####  Species of Interest (w/codes) based on 2020 list
##################################################
spp_list <- sort(c("Atheresthes stomias", "Atheresthes evermanni",
                   "Reinhardtius hippoglossoides",
                   "Lepidopsetta polyxystra", "Lepidopsetta bilineata",
                   "Gadus macrocephalus", "Gadus chalcogrammus",
                   "Pleurogrammus monopterygius", 
                   "Sebastes aleutianus", "Sebastes melanostictus", 
                   "Sebastes alutus", "Sebastes polyspinis"))
spp_list_code <- 
  AI$species$SPECIES_CODE[match(spp_list, AI$species$SPECIES_NAME)]

##################################################
####  Age data is located in the specimen slot in AI
####      with SPECIMEN_SAMPLE_TYPE == 1
##################################################
specimen <- subset(x = AI$specimen,
                   SPECIMEN_SAMPLE_TYPE == 1 & SPECIES_CODE %in% spp_list_code)

total_ages <- table(specimen$SPECIES_CODE, specimen$CRUISE)     

spp_age_names <- AI$species$SPECIES_NAME[match(rownames(total_ages), 
                                               AI$species$SPECIES_CODE)]

# final_age_numbers <- total_ages[spp_age, ]
rownames(total_ages) <- spp_age_names

age_haul <- aggregate(SPECIMEN_SAMPLE_TYPE ~ CRUISE + HAULJOIN, data = specimen, FUN = function(x) sum(x == 1) )

aggregate(SPECIMEN_SAMPLE_TYPE ~ CRUISE, data = age_haul, FUN = quantile)

#########################################################

spp_list_code <- AI$species$SPECIES_CODE[match(x = spp_list, AI$species$SPECIES_NAME)]
sub_catch <- subset(AI$catch, subset = SPECIES_CODE %in% spp_list_code)

spp_list[]

nums <- subset(sub_catch, SPECIES_CODE == spp_list_code[12] & CRUISE == "201801")$NUMBER
nums <- c(nums, rep(0, 420 - length(nums) ))

bootstrap_sum <- c()

for (iter in 1:1000) {
  bootstrap_num <- sample(x = nums, replace = TRUE)
  bootstrap_sum[iter] <- sum(sapply(X = bootstrap_num, 
                                    FUN = function(x) ifelse(x >= 5, 5,
                                                             min(x, 0)) ))
}

quantile(bootstrap_sum, probs = c(0.025, 0.50, 0.975))
