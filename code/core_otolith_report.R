##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Post-cruise otolith collection table
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   After the cruise is completed and the data are finalized
##                this table is created to show how many otoliths were
##                collected vs. requested across species.
##
##                COmpleted for AI 2022 BTS but not generalized yet...
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import libraries
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(RODBC)
library(sumfish)
library(readxl)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import collection rules and requests
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
AI_collection <- readxl::read_xlsx(path = "data/AI_core_collections.xlsx", 
                                   sheet = "2022 collection")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import AI cruise data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
AI <- sumfish::getRacebase(year = 2022, survey = "AI")

specimen <- subset(x = AI$specimen,
                   subset = SPECIMEN_SAMPLE_TYPE == 1)
otolith_table <- table(specimen$SPECIES_CODE)
names(otolith_table) <- AI$species$COMMON_NAME[match(names(otolith_table),
                                                     AI$species$SPECIES_CODE)]

otolith_df <- 
  cbind(AI_collection[, 1:4], 
        total_collected = as.numeric(otolith_table[AI_collection$common_name]))

otolith_df$notes <- 
  ifelse(test = with(otolith_df, 
                     (total_collected - total_requested)/
                       total_requested) > 0.1,
         yes = "over", 
         no = ifelse(test = with(otolith_df, 
                                 ((total_collected - total_requested)/
                                   total_requested) > -0.1),
                     yes = "within 10%", 
                     no = "under"))

otolith_df <- otolith_df[with(otolith_df, 
                              order(total_requested - total_collected)), ]

if(!dir.exists("reports/")) dir.create("reports/")
write.csv(x = otolith_df, 
          file = "reports/202201_AI_otolith.csv", 
          row.names = F)
