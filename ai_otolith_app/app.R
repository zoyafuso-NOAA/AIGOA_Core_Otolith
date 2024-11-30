##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Aleutian Islands Otolith Negotiation Planning Shiny App
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
# ai_raw_data <- readRDS(file = "ai_otolith_app/ai_raw_data.RDS")
# ai_data <- readRDS(file = "ai_otolith_app/ai_data.RDS")
# allocation <- readRDS(file = "ai_otolith_app/ai_station_allocation.RDS")

ai_raw_data <- readRDS(file = "ai_raw_data.RDS")
ai_data <- readRDS(file = "ai_data.RDS")
allocation <- readRDS(file = "ai_station_allocation.RDS")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_codes <- as.character(x = ai_raw_data$species$SPECIES_CODE)
names(x = species_codes) <- ai_raw_data$species$COMMON_NAME
years <- ai_raw_data$survey_design$YEAR
strata <- ai_raw_data$strata$STRATUM

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   User Interface
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui <- shiny::fluidPage(
    
    ## Main Title
    shiny::titlePanel(title = paste0("AI Bottom Trawl Survey Otolith ",
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
            
            shiny::textInput(inputId = "collection_number_wai", 
                             value = 5,
                             label = "Western Aleutian Islands: Collection Number per Haul (Random Sample)"),
            shiny::textInput(inputId = "collection_number_cai", 
                             value = 5,
                             label = "Central Aleutian Islands: Collection Number per Haul (Random Sample)"),
            shiny::textInput(inputId = "collection_number_eai", 
                             value = 5,
                             label = "Eastern Aleutian Islands: Collection Number per Haul (Random Sample)"),
            shiny::textInput(inputId = "collection_number_sbs", 
                             value = 5,
                             label = "Southern Bering Sea: Collection Number per Haul (Random Sample)"),
            
            shiny::textInput(inputId = "crit_val", 
                             label = paste0("Threshold of Collection (i.e, only ",
                                            "collect if the catch number is >= ",
                                            "this number)"), 
                             value = 0),
            shiny::textInput(inputId = "n_boot", 
                             value = 10,
                             label = "Number of Bootstrap Replicates"),
            shiny::textInput(inputId = "seed_number", 
                             value = 2024,
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
    
    input_collection_wai <- 
        shiny::eventReactive(eventExpr = input$goButton, 
                             valueExpr = as.integer(x = input$collection_number_wai))
    input_collection_cai <- 
        shiny::eventReactive(eventExpr = input$goButton, 
                             valueExpr = as.integer(x = input$collection_number_cai))
    input_collection_eai <- 
        shiny::eventReactive(eventExpr = input$goButton, 
                             valueExpr = as.integer(x = input$collection_number_eai))
    input_collection_sbs <- 
        shiny::eventReactive(eventExpr = input$goButton, 
                             valueExpr = as.integer(x = input$collection_number_sbs))
    
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
        shiny::reactive(x = subset(x = ai_data,
                                   select = c("YEAR", "REGION", "STRATUM",
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
            
            collection_by_region <- 
                do.call(
                    what = rbind,
                    args = lapply(
                        X = split(x = queried_catch()[boot_idx, ],
                                  f = list(queried_catch()[boot_idx, "YEAR"], 
                                           queried_catch()[boot_idx, "REGION"])),
                        FUN = function(df) {
                            input_collection <- 
                                switch(unique(x = df$REGION),
                                       "SBS" = input_collection_sbs(),
                                       "Eastern" = input_collection_eai(),
                                       "Central" = input_collection_cai(),
                                       "Western" = input_collection_wai())
                            
                            data.frame(
                                REGION = unique(x = df$REGION),
                                YEAR = unique(x = df$YEAR),
                                TOTAL = sum(sapply(X = df[, input_species()],
                                                   FUN = function(y)
                                                       ifelse(test = y < input_crit_val(),
                                                              yes = 0,
                                                              no = min(y, input_collection))), 
                                            na.rm = TRUE))
                            
                            
                        }))
            
            bootstrap_df[iboot, ] <-
                tapply(X = collection_by_region$TOTAL, 
                       INDEX = collection_by_region$YEAR, 
                       FUN = sum)
        } ## Loop over iterations -- end
        bootstrap_df
    })
    
    output$summary <- shiny::renderPlot({
        par(mar = c(6, 6, 3, 1))
        boxplot(bootstrap_df(), names = years, col = "white", pch = 16,
                las = 1, xlab = "Year", ylab = "Total Collected Otoliths",
                main = ai_raw_data$species$COMMON_NAME[
                    ai_raw_data$species$SPECIES_CODE == input_species()
                ] )
        abline(h = input_target(), lty = "dotted", col = "darkgrey", lwd = 2)
    }, 
    width = 700)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Launch application
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shiny::shinyApp(ui, server)
