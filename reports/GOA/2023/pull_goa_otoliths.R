##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       2023 GOA otolith collections
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Pull 2023 GOA otoliths and save locally 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())


library(gapindex)
sql_channel <- gapindex::get_connected()

observed_collection <- 
  RODBC::sqlQuery(channel = sql_channel, 
                  query = "SELECT SPECIES_CODE, COUNT(SPECIES_CODE) FREQ
                        FROM RACEBASE.SPECIMEN
                        WHERE REGION = 'GOA' AND CRUISE = 202301
                        GROUP BY SPECIES_CODE")

observed_collection <- 
  merge(x = observed_collection,
        y = RODBC::sqlQuery(channel = sql_channel, 
                            query = "SELECT SPECIES_CODE, SPECIES_NAME, COMMON_NAME
                        FROM RACEBASE.SPECIES"),
        by = "SPECIES_CODE")

# Save
names(x = observed_collection) <- tolower(names(x = observed_collection))
saveRDS(object = observed_collection, 
        file = "reports/GOA/2023/observed_otoliths.RDS")
