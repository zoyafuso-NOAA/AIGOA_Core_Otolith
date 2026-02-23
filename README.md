
# **Aleutian Islands/Gulf of Alaska Otolith Planning Tool**

This repo holds the code that maintains an RShiny application that allows 
users to tune the otolith collection rules (how many to collect, when to 
collect) to meet an acquired target total for a given species. 
Recent historical survey data (e.g., last three survey years) are bootstrapped
and the user-defined otolith collection rules are applied to the observed
numbers of fish caught. The summary output is a boxplot time series showing 
the distribution of total collected otoliths across bootstrapped samples. 

## Pre-survey season maintenance

1) Create the RShiny data: Update the code/get_ai_data.R or code/get_goa_data.R
script to pull from the three most recent years of survey data. Make sure the
resultant data outputs live in the goa_otolith_app/ or ai_otolith_app directory
depending on the regioin. These data objects include haul-level count data for
the set of species and the station allocation across strata. 

2) Update code in the app.R scripts change the file names for the survey
data and station allocation to alter aspects of the shiny app itself.

3) Publish the apps to the shinyapp.io dashboard. We currently don't have an
afsc-gap-products account so currently it is through @zoyafuso-NOAA 's
GitHub account. 

4) Create a GitHub release after otolith negotiations are completed to freeze the
code. 
