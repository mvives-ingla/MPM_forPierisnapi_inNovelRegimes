# Builds a bootstrapped transition matrix from raw data (also bootstrapping juvenile and adult longevities)
# INPUT:
## data_juv: data used to estimate juvenile survival and eclosion probability
## data_repr: data used to estimate adult fecundity
## data_ad: data_used to estimate adult survival
## brood_period: simulated phenological period of the population (from brood 1 to 2, or onwards)
# OUTPUT: list with
## A: bootstrapped transition matrix
## dims: the dimensions of the transition matrix (m: max juvenile age, a: max adult age)
## fitted_juv: simulated survival and eclosion probabilities of juveniles
## sample_ref: bootstrap initial sample of juvenile cases

simul_matrix <- function(data_juv, data_repr, data_ad, brood_period) {
  # data_juv <- data_exp
  # data_repr <- data_art
  # data_ad <- data_adsurv
  # brood_period <- "g2g3"
  
  
  ## setting the parameters of the simulation depending on the brood period
  if (brood_period == "g1g2") {
    exp_treat <- 20 # thermal treatment of the experiments
    development <- "diap" # developmental path of the adults (diapausing vs directly developed)
    ad_long <- round(runif(1, 20, 30)) # maximum adult longevity (taken from adult fecundity data)
    juv_long <- round(runif(1, 25, 45)) # maximum juvenile longevity (taken from juvenile data)
    f_surv <- n ~ poly(day, 13) # survival function of the juveniles (the same used in raw data)
    f_repr <- daily_eggs ~ poly(day, 3) + 
      (1 | source) +
      (1 | mat_system) # fecundity function (the same used in raw data)
  } else {
    exp_treat <- 25
    development <- "direct"
    ad_long <- round(runif(1, 20, 30))
    juv_long <- round(runif(1, 25, 45))
    f_surv <- n ~ poly(day, 5)
    f_repr <- daily_eggs ~ poly(day, 5) +
      (1 | source) +
      (1 | mat_system)
  }

  ir_check <- "character" # Check parameter of the irreductibility, primitivity & ergocidity of the matrix
  
  while (ir_check == "character") {
    control_surv <- "try-error" # check parameter of the survival model
    control_ecl <- "try-error" # # check parameter of the eclosion model
  
    while(any(control_surv == "try-error" | control_ecl == "try-error")){
      # step 1: sampling from juvenile data
      sample_juv <- data_juv %>% 
        filter(Treatment == exp_treat) %>% 
        slice_sample(prop = 1, replace = T) %>% 
        rownames_to_column(var = "new_id")
      
      # step 2.1: building the GLM of juvenile survival
      count_data <- sample_juv %>% 
        mutate(alive = 1,
               day = 1) %>%
        # I don't distinguish individuals that die from those that arrive to adult stage
        # I assume that pupal_cycle includes the first day when an adults has been found (as larval cycle includes the first day when a pupa has been found)
        split(.$new_id) %>% 
        map(~ add_row(.x,
                      new_id = rep(.$new_id, .$total_cycle-1),
                      unique_id = rep(.$unique_id, .$total_cycle-1),
                      total_cycle = rep(.$total_cycle, .$total_cycle-1),
                      alive = rep(1, .$total_cycle-1),
                      day = 2:.$total_cycle)) %>% # a row per eah day that each individual is alive
        bind_rows() %>% 
        group_by(day) %>% 
        summarise(n = sum(alive)) # Number of  still alife individuals per day
        
      mod_surv <- try(glm(formula = f_surv,
                      family = "poisson",
                      data =  count_data),
                      silent = F)
      
      control_surv <- class(mod_surv)
      
      # step 2.2: building the GLM of juvenile eclosion
      ecl_data <- sample_juv %>% 
        filter(PtoA)
      
      mod_ecl <- try(glm(formula = (total_cycle-1)~1,
                         family = "poisson",
                         data = ecl_data),
                     silent = F)
      
      control_ecl <- class(mod_ecl)
      }
    
    
    # step 3: predicting juvenile survival and eclosion
    fitted_juv <- count_data %>% 
      mutate(fitted = mod_surv$fitted.values,
             fitted_next = lead(fitted),
             surv = fitted_next/fitted,
             lambda = predict(mod_ecl,
                              newdata = slice_head(ecl_data),
                              type = "response"),
             ecl = dpois(day, lambda = lambda)) %>% 
      filter(!is.na(surv))
    
    ## adding more days in the predicted juvenile parameters
    ## in case the bootstrap sample had less days than the real one
    maxjuv <- max(fitted_juv$day)
    if(maxjuv < juv_long) { 
      fitted_juv <- fitted_juv %>% 
        add_row(day = ((maxjuv+1):juv_long),
                surv = rep(fitted_juv$surv[maxjuv],
                           times = length((maxjuv+1):juv_long)),
                lambda = rep(fitted_juv$lambda[maxjuv],
                             times = length((maxjuv+1):juv_long)),
                ecl = dpois(x = (maxjuv+1):juv_long,
                            lambda = fitted_juv$lambda[maxjuv]))
    } else {
      fitted_juv <- fitted_juv %>% 
        filter(day <= juv_long)
    }
    
    
    
    # step 4: sampling from reproductive data
    repr_threshold <- 1
    while (repr_threshold > 0) {
      
      control_repr <- "try-error"
      while(control_repr == "try-error") {
        
        suppressMessages(
          sample_repr <- data_repr %>% 
            filter(dvt == development) %>% 
            distinct(article_fig, dvt) %>% 
            slice_sample(prop = 1, replace = T) %>% # fer sampling de les figures font de dades?
            mutate(source = row_number()) %>% 
            left_join(data_repr, relationship = "many-to-many")
        )
      
        # step 5: building the model
        mod_repr <- try(glmmTMB(formula = f_repr,
                                family = "nbinom2",
                                data = sample_repr),
                        silent = T)
        
        control_repr <- class(mod_repr)
        }
      
      # step 6: predicting reproductive output
      predicted_repr <- predict(mod_repr,
                                newdata = data.frame(day = 1:ad_long,
                                                     mat_system = NA,
                                                     source = NA),
                                type = "response")
      maxad <- max(sample_repr$day)
      if(maxad < ad_long) {
        predicted_repr[maxad:ad_long] <- predicted_repr[maxad]
        }
      
      fitted_repr <- data.frame(day = 1:ad_long,
                                eggs = round(predicted_repr))
      
      ## avoiding fecundities > 100 eggs/day
      repr_threshold <- fitted_repr %>% 
        filter(eggs > 100) %>% 
        nrow()
    }

    
    # step 7: bootstrapping adult survival 
    
    
    ### from values of P. oleracea taken from Fig. S2 of Kerr 2020
    fitted_ad <- data_ad %>%
      filter(day <= ad_long) %>%
      mutate(mean_logit = log(mean/(1-mean)),
             lowint_logit = log((mean-ci)/(1-mean+ci)),
             ci_logit = mean_logit-lowint_logit,
             ad_surv_logit = rnorm(ad_long,
                                   mean = mean_logit,
                                   sd = ci_logit/2),
             ad_surv = 1/(1+exp(-ad_surv_logit)))
  
    # step 8: building the matrix
    A <- build_matrix(juvsurvs = fitted_juv$surv,
                      juvecls = fitted_juv$ecl,
                      adfecs = fitted_repr$eggs,
                      adsurvs = fitted_ad$ad_surv)
    
    
    
    # step 9: checking irreductibility
    ir_check <- class(A)[1]
    }
  
  return(list(A = A, # the simulated transition matrix
              dims = c(m = juv_long, a = ad_long), # the dimensions of the transition matrix (m: max juvenile age, a: max adult age)
              fitted_juv = fitted_juv, # simulated survival and eclosion probabilities
              sample_ref = sample_juv)) # bootstrap sample of juvenile cases
}
