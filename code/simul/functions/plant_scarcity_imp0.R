# Bootstrapping the effects of plant scarcity and adding them on the bootstrapped transition matrix
# Effects on juveniles survival are added as a ratio of survivals between the period with no plant and the reference period
# These ratios are estimated for larvae that were 1 week age when plant scarcity occurred, 2 weeks or 3 weeks
# Effects on juveniles eclosion are added calculating the difference of mean day of larval pupation between the two treatments and substracting it from the reference pupation date
# INPUT
## A: list with the reference transition matrix and the bootstrapped sample
## data_plant: raw dataset with juvenile information during the experiment with plant scarcity
# OUTPUT
## Bootstrapped transition matrix with plant scarcity effects

## Imposing a minimum age for pupal eclosion of 8 days

plant_scarcity <- function(A, data_plant) { 
 
  # data_plant <- data$plsc
  # A <- sims$A[[1]]
  # plant <- sims$sim_df[[1]]$no_more_plant
  
  # Step 1: bootstrap sample and calculation of plant scarcity effects
  sample_plant <- data_plant %>%
    mutate(age_fatal = as.numeric(ymd("20151026")-as.Date(Egg_laying, tz = "UTC")), # what age individuals were when plant scarcity started
           plant_group = case_when(age_fatal < 8 ~ "w1", #it includes individuals that were age 6 and 7
                                   age_fatal < 16 ~ "w2", # ages 12, 13 and 15
                                   age_fatal < 24 ~ "w3", # ages 17, 18 and 19
                                   T ~  NA_character_), # at which week of dvt individuals were when plant scarcity started
           plant_group = factor(plant_group, levels = c("w1", "w2", "w3")),
           period = "plantscarc",
           alive = 1,
           total_cycle = Egg_period + Larval_period) %>% # only considering eggs and larvae, as there's no information on pupation for this experiment
    select(-c(Egg_laying, age_fatal)) %>% 
    filter(!is.na(plant_group), !is.na(total_cycle)) %>% 
    group_by(plant_group) %>% 
    slice_sample(prop = 1, replace = T)
  
  suppressMessages(
    sample <- data.frame(plant_group = c("w1", "w2", "w3")) %>% 
      mutate(ref_df = map(plant_group,
                         ~ mutate(A$sample_ref, # the bootstrapped data for the reference matrix, repeated per each considered larval period
                                  alive = 1,
                                  period = "reference",
                                  total_cycle = Egg_period + Larval_period, # recalculating total cycle only with eggs and larvae
                                  new_id = NULL,
                                  Treatment = NULL,
                                  PtoA = NULL))) %>% 
      unnest(ref_df) %>% 
      bind_rows(sample_plant) %>% 
      select(-c(Egg_period, Larval_period)) %>% 
      rownames_to_column(var = "new_id") %>% 
      nest(df = -(period)) %>% 
      mutate(ecl = map(df, ~ distinct(., new_id, LtoP, total_cycle)), # data for the estimation of larval pupation
             ecl = map(ecl, filter, LtoP),
             ecl = map(ecl, ~ mutate(., total_cycle = total_cycle - 1)),
             ecl = map(ecl, ~ glm(total_cycle ~ 1, family = poisson, data = .)),
             lambda = map(ecl, broom::tidy), # mean day of larval pupation
             lambda = map(lambda, select, estimate),
             lambda = map(lambda, exp), # back-transformation from log distribution
             lambda = map(lambda, as.numeric),
             ecl = NULL) %>% 
      unnest(lambda) %>% 
      unnest(df) %>% 
      mutate(day = map(total_cycle, ~ data.frame(day = 1:.))) %>% # creating a dataset with an entry per individual and day
      unnest(day) %>% 
      mutate(ini = case_when(plant_group == "w1" ~ 6,
                             plant_group == "w2" ~ 12,
                             T ~ 17),
             end = case_when(plant_group == "w1" ~ 14,
                             plant_group == "w2" ~ 21,
                             T ~ 28)) %>% # initial and ending day of the different larval periods in which we will estimate larval survival
                                          # i.e. from the minimum larval age of the group to the end of the preceding week (after this first week we can't guarantee that plant scarcity persisted)
      filter(day >= ini,
             day <= end) %>%
      nest(count_data = -c(period, plant_group, lambda)) %>% 
      mutate(count_data = map(count_data, select, day, alive),
             count_data = map(count_data, group_by, day),
             count_data = map(count_data, summarise, n = sum(alive)), # survival curves (N ~ day)
             model = map(count_data, ~ glm(n ~ day, family = "poisson", data = .)),
             surv = map(model, broom::tidy),
             count_data = NULL,
             model = NULL,
             surv = map(surv, filter, term == "day"), # daily predicted survival
             surv = map(surv, select, estimate),
             surv = map(surv, as.numeric),
             surv = map(surv, exp)) %>% # back-transformation from the log scale
      unnest(surv) %>% 
      pivot_wider(id_cols = plant_group,
                  names_from = period,
                  values_from = c(surv, lambda)) %>%
      mutate(ratio_surv = surv_plantscarc/surv_reference,
             dif_lambda = lambda_reference-lambda_plantscarc)
  )

  # Step 2: Adding bootstrapped plant scarcity effects on bootstrapped transition matrix
  ## maximum juvenile age
  m <- nrow(A$fitted_juv)
  ## maximum adult age
  a <- ncol(A$A) - m
  
  surv_rats <- c(rep(sample$ratio_surv[1], 13),
                 rep(sample$ratio_surv[2], length(14:20)),
                 rep(sample$ratio_surv[3], length(21:m)))
  
  
  fitted_juv_plant <- A$fitted_juv %>%
    mutate(lambda = lambda-sample$dif_lambda[1],
           ecl = dpois(day, lambda = lambda),
           surv = surv*surv_rats,
           surv = if_else(surv > 1, 1, surv)) #when surv > 1, in the matrix I do a correction so that S + E = 1, which is equivalent to surv = 1
                                              # when surv <= 1, it's impossible that S + E > 1
  
  A_plant <- A$A
  for (i in 2:(m+1)) {
    if(i <= m) {
      #surviving and resting as juveniles
      A_plant[i, i-1] <- fitted_juv_plant$surv[i-1]*(1-fitted_juv_plant$ecl[i-1])
    } 
    if(i == (m+1)) {
      #surviving and ecloding
      A_plant[i, 8:(m)] <- fitted_juv_plant$surv[8:m]*fitted_juv_plant$ecl[8:m]
    }
  }
  
  return(A_plant)

}
