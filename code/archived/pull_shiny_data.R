##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Data prep for Shiny App
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Subset relevant historical data for otolith shiny app
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

## Import gapindex package
library(gapindex)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import historical datasets
##   Calculate and zero-fill cpue
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df <- readRDS("C:/Users/zack.oyafuso/Desktop/gapindex_data/GOA_DATA.RDS")
# cpue <- gapindex::calc_cpue(racebase_tables = df)
# saveRDS(object = cpue, file = "data/GOA/2023/GOA_CPUE.RDS")
cpue <- readRDS("data/GOA/2023/GOA_CPUE.RDS")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Extract objects from df
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
specimen <- df$specimen
cruise <- df$cruise
species <- df$species
lengths <- df$size

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Attach year and common names to specimen, cpue, and lengths, 
##   Subset to only years >= minimum year
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
min_year <- 2017

specimen <- merge(x = specimen,
                  y = cruise[, c("CRUISEJOIN", "YEAR")],
                  by = "CRUISEJOIN")
specimen <- merge(x = specimen,
                  y = species[, c("SPECIES_CODE", "COMMON_NAME")],
                  by = "SPECIES_CODE")
specimen <- subset(x = specimen, 
                   subset = YEAR >= min_year)

cpue <- subset(x = cpue,
               subset = SPECIES_CODE %in% unique(specimen$SPECIES_CODE) &
                 YEAR >= min_year, 
               select = c(SURVEY, YEAR, STRATUM, SPECIES_CODE, COUNT))
cpue <- merge(x = cpue,
              y = species[, c("SPECIES_CODE", "COMMON_NAME")],
              by = "SPECIES_CODE")

lengths <- merge(x = lengths,
                 y = cruise[, c("CRUISEJOIN", "YEAR")],
                 by = "CRUISEJOIN")
lengths <- subset(x = lengths,
                  subset = SPECIES_CODE %in% unique(specimen$SPECIES_CODE) &
                    YEAR >= min_year, 
                  select = c(YEAR, SPECIES_CODE, SEX, LENGTH, FREQUENCY))
save(list = c("cpue", "lengths","specimen", "species"), 
     file = "data/GOA/2023/GOA_shiny_data.RData")
