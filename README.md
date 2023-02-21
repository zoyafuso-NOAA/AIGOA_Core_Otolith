
# AI/GOA Core otolith request analysis

## [Core Otolith Requests](#top)

The goal is to assess the feasibility of the core otolith request using
haul data from the past five survey years. Requests for each species
include collection rules and a total otolith pair requested. The main
questions that drive the feasibility of the requeted core otolith
requests are:

1)  Given the requested collection rules, how well could the survey
    achieve the requested target across species?
2)  How many otolith pairs are collected during each haul?
3)  Given 30 otolith pairs are a soft maximum haul-total, how often are hauls
    over this threshold?

For the last two survey years, hauls were sampled with replacement
until a bootstrapped sample of 520 hauls were obtained. The collection
rule was applied to each haul and the total number of otoliths collected
was totaled. This was repeated 1000 times to obtain a distribution of
total otoliths collected conditional on the collection rule requested
this year.

## Repo structure

data/ : contains the station allocation by stratum, collection requests, and historical survey data for either AI or GOA.
code/ : contains the scripts used for the analysis
reports/ : contains the RMD scripts used to generate post-cruise reports for the stock assessors. 

## Useful links

link to GOA otolith roundup

link to CORE otolith requests

links to 

