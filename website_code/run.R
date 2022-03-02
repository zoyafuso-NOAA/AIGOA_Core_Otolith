library(rmarkdown)
library(readxl)
library(glue)

###############################################################################
####   Load bootstrap results and collection spreadsheet
###############################################################################
load("results/bootstrap_results.RData")
collection <- readxl::read_xlsx("data/AI_core_collections.xlsx", 
                                sheet = "2022 collection")

##################################################
####   Render each species page
##################################################  
for (irow in 1:(nrow(collection) - 1)) {
  species_name <- collection$common_name[irow]
  rmarkdown::render(input = "website_code/core_otolith_reports.Rmd", 
                    output_file = paste0(species_name, ".html"), 
                    output_dir = "docs/", 
                    output_format = "html_document")
}

##################################################
####   Render Main Page
##################################################  
rmarkdown::render(input = "website_code/index.Rmd", 
                  output_file = "index.html", 
                  output_dir = "docs/", 
                  output_format = "html_document")
