
# **Aleutian Islands/Gulf of Alaska Otolith Planning Tool**

## **For Users**
This repo holds the code that maintains an RShiny application that allows 
users to tune the otolith collection rules (how many to collect, when to 
collect) to meet an acquired target total for a given species. 
Recent historical survey data (e.g., last three survey years) are bootstrapped
and the user-defined otolith collection rules are applied to the observed
numbers of fish caught. The summary output is a boxplot time series showing 
the distribution of total collected otoliths across bootstrapped samples. 

### [Otolith Planning RShiny Tool](https://zoyafuso-noaa.shinyapps.io/aigoa_otolith_planning_tool/)

## **For Repo Maintainer**

### Pre-survey season maintenance:

These are the steps to maintain the tool every year:

1) Create the RShiny data: Run the code/get_ai_data.R script for a given 
region for the three most recent years of survey data. Make sure the resultant 
data output lives in the ai_otolith_app/ directory. These data objects include
haul-level count data for the set of species and the station allocation across
strata. 

3) In shiny/app.R change the file names for the survey data and station 
allocation to whatever is in the shiny/ directory. 

### Pre-season report

Use the script in code/bootstrapped_otolith_collections.R to gauge three 
aspects of otolith planning.

1)  Given the requested collection rules, how well could the survey
    achieve the requested target across species?

2)  How many otolith pairs are collected during each haul?

3)  Given 30 otolith pairs are a soft maximum haul-total, how often are hauls
    over this threshold?
