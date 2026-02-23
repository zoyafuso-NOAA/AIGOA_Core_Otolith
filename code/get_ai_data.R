##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Create input data for shiny app
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Pull specimen data for AI otolith shiny app
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

library(tidyr)
library(gapindex)
library(xlsx)
channel <- gapindex::get_connected(check_access = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Get Upcoming Station Allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# station_allocation <- 
#   RODBC::sqlQuery(channel = channel,
#                   query = "SELECT * FROM AI.STATION_ALLOCATION
#                            WHERE YEAR = 2024
#                            AND STATION_TYPE != 'bonus_stn'")
# station_allocation <- 
#   data.frame(STRATUM = names(x = table(station_allocation$STRATUM)),
#              N_HAUL = as.integer(x = table(station_allocation$STRATUM)))

station_allocation <- 
  xlsx::read.xlsx(file = paste0("Y:/RACE_GF/ALEUTIAN/AI 2026/",
                                "Station Allocation/",
                                "ai_2026_station_allocation_400stn.xlsx"),
                  sheetName = "Station Allocation") |>
  aggregate(VESSEL  ~ STRATUM,
            FUN = length) |>
  transform(N_HAUL = VESSEL) |>
  subset(select = c(STRATUM, N_HAUL))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Get haul-level count data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
spp_codes <- c(10110, 10112, 10115, 10120, 10261, 10262,
               21720, 21740, 21921,
               30020, 30051, 30052, 30060, 30152, 30420, 30535, 30576)

ai_data <- gapindex::get_data(year_set = c(2018, 2022, 2024),
                              survey_set = "AI",
                              spp_codes = spp_codes, 
                              pull_lengths = F, 
                              taxonomic_source = "RACEBASE.SPECIES",
                              channel = channel)
ai_cpue <- gapindex::calc_cpue(gapdata = ai_data)
ai_cpue$REGION <- sapply(X = substr(x = ai_cpue$STRATUM,
                                    start = 1,
                                    stop = 1),
                         FUN = function(x)
                           switch(x,
                                  "2" = "Western",
                                  "3" = "Central", "4" = "Central",
                                  "5" = "Eastern", "6" = "Eastern",
                                  "7" = "SBS"))

data_wide <- 
  tidyr::spread(data = subset(x = ai_cpue, 
                              select = c(YEAR, REGION, SPECIES_CODE, 
                                         STRATUM, HAULJOIN, COUNT)),
                value = COUNT, key = SPECIES_CODE, fill = 0)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(object = ai_data,
        file = "ai_otolith_app/ai_raw_data.RDS")
saveRDS(object = data_wide,
        file = "ai_otolith_app/ai_data.RDS")
saveRDS(object = station_allocation,
        file = "ai_otolith_app/ai_station_allocation.RDS")
