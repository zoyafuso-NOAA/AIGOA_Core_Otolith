##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Otolith Negotiation Planning Shiny Application
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   This shiny application is used as a tool for stock assessment
##                requestors and survey planners to  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

## Import shiny library
library(shiny)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import data. Make sure file names match those in the shiny/ directory
##   Import stratum allocations for a given year, modify for given year/survey
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
region <- c("GOA", "AI")[1]
load("GOA_shiny_data.RData")
allocation_table <- read.csv(file = "GOA2023_Station_allocation_520_EW.csv")
allocation <- table(allocation_table$stratum)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_codes <- sort(x = unique(x = specimen$SPECIES_CODE))
common_name <- 
  sort(x = species$COMMON_NAME[match(x = species_codes,
                                     table = species$SPECIES_CODE)])
years <- sort(x = unique(x = specimen$YEAR))
strata <- names(x = allocation)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   User Interface
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui <- shiny::fluidPage(
  
  ## Main Title
  shiny::titlePanel(title = paste0(region, " Bottom Trawl Survey Otolith ",
                                   "Collection Planner")),
  
  ## All the user inputs will live on the left sidebar of the app
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::selectInput(inputId = "species", 
                         label = "Choose a Species", 
                         choices = sort(x = common_name), 
                         selected = "Pacific ocean perch"),
      shiny::textInput(inputId = "target", 
                       label = "Target Total Number of Otoliths", 
                       value = 1000),
      shiny::textInput(inputId = "collection_number", 
                       value = 5,
                       label = "Collection Number per Haul (Random Sample)"),
      shiny::textInput(inputId = "crit_val", 
                       label = paste0("Threshold of Collection (i.e, only ",
                                      "collect if the catch number is >= ",
                                      "this number)"), 
                       value = 5),
      shiny::textInput(inputId = "n_boot", 
                       value = 10,
                       label = "Number of Bootstrap Replicates"),
      
      ## Submit button
      shiny::actionButton(inputId = "goButton", 
                          label = "Update", 
                          class = "btn-success"),
      width = 6
    ),
    
    ## The resulting bootstrap boxplots will live on the right side of the app
    shiny::mainPanel(
      shiny::plotOutput(outputId = "summary"),
      width = 5
    )
  )
)

server <- function(input, output, session) {
  
  ## Update user input
  input_species <- 
    shiny::eventReactive(eventExpr = input$goButton, 
                         valueExpr = input$species)
  input_collection <- 
    shiny::eventReactive(eventExpr = input$goButton, 
                         valueExpr = as.integer(x = input$collection_number))
  input_crit_val <- 
    shiny::eventReactive(eventExpr = input$goButton, 
                         valueExpr = as.integer(x = input$crit_val))
  input_nboot <- 
    shiny::eventReactive(eventExpr = input$goButton, 
                         valueExpr = as.integer(x = input$n_boot))
  input_target <- 
    shiny::eventReactive(eventExpr = input$goButton, 
                         valueExpr = as.integer(x = input$target))
  
  ## Subset `specimen` and `cpue` dfs according to the input species
  queried_specimen <-
    shiny::reactive(x = subset(x = specimen,
                               subset = COMMON_NAME == input_species(),
                               select = c(YEAR, REGION, SPECIES_CODE, AGE)))
  
  queried_catch <-
    shiny::reactive(x = subset(x = cpue,
                               subset = COMMON_NAME == input_species(),
                               select = c(YEAR, SURVEY, STRATUM, 
                                          SPECIES_CODE, COUNT)))
  
  ## For each iteration, bootstrap hauls according to the strata allocations,
  ## apply the collection rules specified in the user input, and plot the
  ## distribution of total collected otoliths. 
  bootstrap_df <- shiny::reactive({
    
    bootstrap_df <- matrix(nrow = input_nboot(),
                           ncol = length(years))
    
    for (iboot in 1:input_nboot()){ ## Loop over iterations -- start
      
      #Empty vector of indices
      boot_idx <- c()
      
      for (iyear in years) { ## Loop over years -- start
        for (istratum in strata) { ## Loop over strata -- start
          
          boot_idx <- 
            c(boot_idx,
              sample(x = which(queried_catch()[, "YEAR"] == iyear &
                                 queried_catch()[, "STRATUM"] == istratum),
                     size = allocation[istratum],
                     replace = TRUE))
          
        } ## Loop over strata -- end
      } ## Loop over years -- end
      
      ## Optimize collection quantity.
      bootstrap_df[iboot, ] <-
        tapply(X = queried_catch()[boot_idx, "COUNT"],
               INDEX = queried_catch()[boot_idx, "YEAR"],
               FUN = function(x) 
                 sum(sapply(X = x, 
                            FUN = function(y) 
                              ifelse(test = y < input_crit_val(),
                                     yes = 0,
                                     no = min(y, input_collection()))),
                     na.rm = TRUE)
        )
    } ## Loop over iterations -- end
    bootstrap_df
  })
  
  output$summary <- shiny::renderPlot({
    par(mar = c(6, 6, 3, 1))
    boxplot(bootstrap_df(), names = years, col = "white",
            las = 1, xlab = "Year", ylab = "Total Collected Otoliths")
    abline(h = input_target(), lty = "dotted", col = "darkgrey", lwd = 2)
  }, 
  width = 500)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Launch application
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shiny::shinyApp(ui, server)
