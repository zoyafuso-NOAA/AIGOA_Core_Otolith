##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Create input data for shiny app
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Pull specimen data for AI otolith shiny app
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

library(tidyr)
library(gapindex)
library(StationAllocationAIGOA)
library(terra)

channel <- gapindex::get_connected(check_access = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Get Upcoming Station Allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
station_allocation <- 
  StationAllocationAIGOA::goa_allocate_stations(n = 520)[["ms_allocation"]]
names(x = station_allocation) <- c("STRATUM", "REP_AREA", "N_HAUL")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Get haul-level count data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
spp_codes <- c(10110, 10120, 10130, 10180, 10200, 10261, 10262,
               20510, 21720, 21740, 21921,
               30020, 30051, 30052, 30060, 30100, 
               30152, 30420, 30430, 30475, 30535)

## Pull GOA catch and effort
goa_data <- 
  gapindex::get_data(year_set = c(2019, 2023),
                     survey_set = "GOA",
                     spp_codes = spp_codes, 
                     pull_lengths = F, 
                     channel = channel, 
                     taxonomic_source = "GAP_PRODUCTS.TAXONOMIC_CLASSIFICATION")

## Reclassify stratum of hauls to the 2025 restratification
goa_strata <- terra::vect(x = "data/GOA/shapefiles/goa_strata_2025.gpkg")
goa_pts <- 
  terra::vect(x = goa_data$haul, 
              geom = c("START_LONGITUDE", "START_LATITUDE"), 
              crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
goa_pts <- terra::project(x = goa_pts, "EPSG:3338")

## Intersect the 2023 stations with the 2025 strata
goa_pts_2025_strata <-terra::intersect(x = goa_pts[, c("HAULJOIN")], 
                                       y = goa_strata[, c("STRATUM")])

## Reassign 2023 strata
goa_data$haul <- merge(x = subset(x = goa_data$haul, select = -STRATUM),
                       y = goa_pts_2025_strata, 
                       by = "HAULJOIN")

## Zero fill hauls
goa_cpue <- gapindex::calc_cpue(gapdata = goa_data)

## Spread out data
data_wide <- tidyr::spread(data = subset(x = goa_cpue, 
                                         select = c(YEAR, SPECIES_CODE, 
                                                    STRATUM, HAULJOIN, COUNT)),
                           value = COUNT, key = SPECIES_CODE, fill = 0)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!dir.exists(paths = "goa_otolith_app/")) 
  dir.create(path = "goa_otolith_app/")

saveRDS(object = goa_data,
        file = "goa_otolith_app/goa_raw_data.RDS")
saveRDS(object = data_wide,
        file = "goa_otolith_app/goa_data.RDS")
saveRDS(object = station_allocation,
        file = "goa_otolith_app/goa_station_allocation.RDS")
