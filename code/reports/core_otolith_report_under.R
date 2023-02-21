##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For species that we did not reach the target requested:
##                was this because we left fish "on the table" but not 
##                collecting them because of our sampling
##
##                COmpleted for AI 2022 BTS but not generalized yet...
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import libraries
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(RODBC)
library(getPass)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Connect to VPN!
##   Connect to oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
assign(x = "channel",
       value = RODBC::odbcConnect("AFSC",
                                  uid = getPass::getPass("uid"),
                                  pwd = getPass::getPass("pwd")),
       envir = .GlobalEnv)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import data
##   AI_collection are the 2021 collection rules
##   otolith_report is the report table of totals (requested vs collected)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
AI_collection <- readxl::read_xlsx(path = "data/AI_core_collections.xlsx", 
                                   sheet = "2022 collection")

otolith_report <- read.csv("reports/202201_AI_otolith.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Query length and specimen data for the species that were under-colected
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
under_spp <- otolith_report$species_code[otolith_report$notes == "under"]
under_spp_txt <- 
  paste0("(", 
         paste0(under_spp[-length(under_spp)], sep = ", ", collapse = ""),
         under_spp[length(under_spp)], ")")

obs_len <- RODBC::sqlQuery(
  channel = channel, 
  query = paste("select * from RACEBASE.LENGTH",
                "where", 
                "REGION = 'AI' and",
                "CRUISE = 202201 and",
                "SPECIES_CODE IN", under_spp_txt,
                collapse = ""))

obs_age <- RODBC::sqlQuery(
  channel = channel, 
  query = paste("select * from RACEBASE.SPECIMEN",
                "where", 
                "REGION = 'AI' and",
                "CRUISE = 202201 and",
                "SPECIES_CODE IN", under_spp_txt,
                collapse = ""))

under_lens <- with(obs_len, table(SPECIES_CODE))
saveRDS(object = under_lens, file = "reports/under_lens.rds")
