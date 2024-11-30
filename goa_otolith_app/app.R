##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Gulf of Alaska Otolith Negotiation Planning Shiny App
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   This shiny application is used as a tool for stock assessment
##                requestors and survey planners to set collection rules
##                for otolith specimens. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

## Import shiny library
library(shiny)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import data. Make sure file names match those in the shiny/ directory
##   Import stratum allocations for a given year, modify for given year/survey
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_raw_data <- readRDS(file = "goa_raw_data.RDS")
goa_data <- readRDS(file = "goa_data.RDS")
allocation <- readRDS(file = "goa_station_allocation.RDS")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_codes <- sapply(X = as.character(x = goa_raw_data$species$SPECIES_CODE),
                        FUN = list)
common_names <- goa_raw_data$species$COMMON_NAME
names(x = species_codes) <- paste0(species_codes, ": ", common_names)
years <- goa_raw_data$survey$YEAR
strata <- sort(x = unique(x = goa_raw_data$haul$STRATUM))
strata <- strata[strata < 500]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   User Interface
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui <- shiny::fluidPage(
  
  ## Main Title
  shiny::titlePanel(title = paste0("GOA Bottom Trawl Survey Otolith ",
                                   "Collection Planner")),
  
  ## All the user inputs will live on the left sidebar of the app
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::selectInput(inputId = "species", 
                         label = "Choose a Species", 
                         choices = species_codes, 
                         selected = "30060"),
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
                       value = 0),
      shiny::textInput(inputId = "n_boot", 
                       value = 100,
                       label = "Number of Bootstrap Replicates"),
      shiny::textInput(inputId = "seed_number", 
                       value = 2025,
                       label = "Set a random seed"),
      
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
  input_seed <- 
    shiny::eventReactive(eventExpr = input$goButton, 
                         valueExpr = input$seed_number)
  
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
  queried_catch <-
    shiny::reactive(x = subset(x = goa_data,
                               select = c("YEAR", "STRATUM",
                                          input_species() )))
  
  ## For each iteration, bootstrap hauls according to the strata allocations,
  ## apply the collection rules specified in the user input, and plot the
  ## distribution of total collected otoliths. 
  
  bootstrap_df <- shiny::reactive({
    set.seed(seed = input_seed() )
    bootstrap_df <- matrix(nrow = input_nboot(),
                           ncol = length(years))
    
    for (iboot in 1:input_nboot()){ ## Loop over iterations -- start
      
      #Empty vector of indices
      boot_idx <- c()
      
      for (iyear in years) { ## Loop over years -- start
        for (istratum in strata) { ## Loop over strata -- start
          
          boot_idx <- 
            c(boot_idx,
              sample(x = which(
                queried_catch()[, "YEAR"] == iyear &
                  queried_catch()[, "STRATUM"] == istratum),
                size = allocation$N_HAUL[allocation$STRATUM == istratum],
                replace = TRUE))
          
        } ## Loop over strata -- end
      } ## Loop over years -- end
      
      total_collection <- 
        do.call(
          what = rbind,
          args = lapply(
            X = split(x = queried_catch()[boot_idx, ],
                      f = list(queried_catch()[boot_idx, "YEAR"])),
            FUN = function(df) {
              
              data.frame(
                YEAR = unique(x = df$YEAR),
                TOTAL = sum(sapply(X = df[, input_species()],
                                   FUN = function(y)
                                     ifelse(test = y < input_crit_val(),
                                            yes = 0,
                                            no = min(y, input_collection()))), 
                            na.rm = TRUE))
              
              
            }))
      
      bootstrap_df[iboot, ] <- total_collection$TOTAL
      
    } ## Loop over iterations -- end
    bootstrap_df
  })
  
  observeEvent(input$goButton, {
    output$summary <- shiny::renderPlot({
      par(mar = c(6, 6, 3, 1))
      boxplot(bootstrap_df(), names = years, col = "white", pch = 16,
              las = 1, xlab = "Year", ylab = "Total Collected Otoliths",
              main = goa_raw_data$species$COMMON_NAME[
                goa_raw_data$species$SPECIES_CODE == input_species()
              ] )
      abline(h = input_target(), lty = "dotted", col = "darkgrey", lwd = 2)
    }, 
    width = 500)
  })
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Launch application
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shiny::shinyApp(ui, server)
