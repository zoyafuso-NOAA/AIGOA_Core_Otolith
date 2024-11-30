otolith_opt <- function(collect,
                        species_code,
                        bootstrap_hauls_,
                        crit_val,
                        target_val,
                        max_depth,
                        return_what = c("obj",
                                        "total_by_haul",
                                        "total_by_cruise")[1]) {
  
  ## Result output to put ototliths collected
  bootstrap_otos <- vector(length = nrow(bootstrap_hauls_))
  
  for (ihaul in 1:nrow(bootstrap_hauls_)) { ## Loop over hauls -- start
    
    ## How many ispp_code was caught in ihaul?
    icatch <- bootstrap_hauls_[ihaul, species_code]
    
    ## Should an otolith subsample take place?
    # collect_otos <- TRUE
    collect_otos <-
      icatch >= crit_val & bootstrap_hauls_[ihaul, "DEPTH"] < max_depth
    
    ## If collect_otos is TRUE, collect otoliths based on the threshold
    bootstrap_otos[ihaul] <- ifelse(test = collect_otos == TRUE,
                                    no = 0,
                                    yes = min(icatch, collect))
  } ## Loop over hauls -- end

  total_by_haul <- matrix(data = bootstrap_otos,
                          ncol = length(unique(bootstrap_hauls_$CRUISE)),
                          dimnames = list(NULL, unique(bootstrap_hauls_$CRUISE)))
  
  ## Calculate sum of the absolute difference between the total cruise
  ## collection and the target
  total_by_cruise <- tapply(X = bootstrap_otos,
                            INDEX = bootstrap_hauls_$CRUISE,
                            FUN = sum)
  
  obj <- sum(abs(total_by_cruise - target_val)^2)
  
  return(get(return_what))
}

stratified_samp_1cm <- function(collect,
                                bootstrap_hauls_,
                                species_code,
                                lengths_) {
  
  ## Result output to put ototliths collected
  bootstrap_otos <- vector(length = nrow(bootstrap_hauls_))
  
  for (ihaul in 1:nrow(bootstrap_hauls_)) { ## Loop over hauls -- start
    
    ## How many ispp_code was caught in ihaul?
    icatch <- bootstrap_hauls_[ihaul, species_code]
    
    if (icatch == 0) bootstrap_otos[ihaul] <- 0
    if (icatch > 0) {
      table_sex_length <-
        with(subset(x = lengths_,
                    subset = CRUISE == bootstrap_hauls_$CRUISE[ihaul] &
                      HAULJOIN == bootstrap_hauls_$HAULJOIN[ihaul] &
                      SPECIES_CODE == species_code),
             tapply(X = FREQUENCY, 
                    INDEX = list(SEX, LENGTH), 
                    FUN = sum, 
                    na.rm = TRUE ))
      
      collect_otos <- apply(X = table_sex_length,
                            MARGIN = 1:2,
                            FUN = function(x) ifelse(test = x > 0,
                                                     yes = min(collect, x),
                                                     no = 0))
      
      ## If collect_otos is TRUE, collect otoliths based on the threshold
      bootstrap_otos[ihaul] <- sum(collect_otos, na.rm = TRUE)
    }
  } ## Loop over hauls -- end
  
  total_by_cruise <- tapply(X = bootstrap_otos,
                            INDEX = bootstrap_hauls_$CRUISE,
                            FUN = sum)
  
  total_by_haul <- matrix(data = bootstrap_otos,
                          ncol = length(unique(bootstrap_hauls_$CRUISE)),
                          dimnames = list(NULL, 
                                          unique(bootstrap_hauls_$CRUISE)))
  
  return(list(total_by_cruise = total_by_cruise,
              total_by_haul = total_by_haul))
}
