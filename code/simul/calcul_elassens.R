#! /usr/bin/env Rscript


# Sensitivity and elasticity analysis of the bootstrapped matrices - calculation

# Packages ----------------------------------------------------------------

library(tidyverse)
library(popdemo)
library(furrr)

source("code/functions/change_juvs_mult.R")

plan(multicore, workers = 40)


# Data --------------------------------------------------------------------

## Simulations
load("data/mpm/sims_output/sims_future_drought_0.2-0.25-0.4-0.6_imp0_2024-12-10.RData")
sims_imp0 <- sims_imp0 %>% 
  filter(!is.na(gr)) %>% 
  select(sim_set, sim, sim_with_plant, sim_df, A)



## Functions ------------------------------------------------------------------------
# Calculate sensitivity and elasticity from the sequence of daily matrices of each simulation
calc_everything <- function(A_trans, # A or A_plant of the simulation
                            mort, # sequence of thermal mortalities
                            rob_test, # survival factor after subtracting extra mortality
                            a, # maximum adult longevity
                            m) { # maximum junveile longevity
  
  # daily sequence of transition matrices
  A_seq <- change_juvs_mult(A = A_trans, mort = mort, larvmort_robtest = rob_test) 
  
  # tibble of elasticity and sensitivity matrices of each daily transition matrix
  sens_elas_df <- tibble(days = 1:70,
                         sens = map(A_seq, sens),
                         elas = map(A_seq, elas))
  
  # removing objects that are not needed anymore
  rm(A_seq)                       
  
  # names of the population stages
  stage_names <- c(paste0("j0", 1:9),
                   paste0("j", 10:m),
                   paste0("a0", 1:9),
                   paste0("a", 10:a))
  
  
  sens_elas_df <- sens_elas_df %>% 
    mutate(sens = map(sens,
                      ~ matrix(.x, # adding row and col names
                               nrow = length(stage_names),
                               ncol = length(stage_names),
                               dimnames = list(
                                 stage_names,
                                 stage_names)))) %>% 
    pivot_longer(cols = sens:elas,
                 names_to = "sens_elas_type",
                 values_to = "sens_elas_value") %>% # long df to manipulate elas and sens at the same time
    mutate(sens_elas_value = map(sens_elas_value, as_tibble), # from matrix to df
           sens_elas_value = map(sens_elas_value, 
                                 ~ mutate(.x,
                                          next_stage = stage_names)), #cols: present stage, rows: next stage
           sens_elas_value = map(sens_elas_value,
                                 ~ pivot_longer(.,
                                                cols = -next_stage,
                                                names_to = "stage",
                                                values_to = "value"))) %>%  # from matrix structure to long df 
    unnest(sens_elas_value) %>% 
    mutate(main_stage = str_sub(stage, start = 1, end = 1),
           main_next_stage = str_sub(next_stage, start = 1, end = 1),
           # classifying transitions by component (juv survival, eclosion, ad survival or fecundity)
           component = case_when(main_stage == "j" & main_next_stage == "j" ~ "juv_surv",
                                 next_stage == "a01" ~ "juv_ecl",
                                 next_stage == "j01" ~ "fec",
                                 T ~ "ad_surv"),
           component = factor(component, levels = c("juv_surv", "juv_ecl", "fec", "ad_surv"))) %>% 
    group_by(days, sens_elas_type, component) %>% 
    summarise(value = sum(value)) %>% 
    ungroup()
  
  return(sens_elas_df)
}




gc() #liberating space


suppressMessages(future_elas_sens <- sims_imp0 %>%  
                   unnest(sim_df) %>% 
                   select(-c(days:tx40)) %>% #deleting columns that are not needed
                   nest(mort = mort) %>%  # a row per simulation, so mortality vectors are nested
                   mutate(mat_type = case_when(str_detect(sim_set, "g1g2") ~ "g1g2",
                                               sim_with_plant == 1 ~ "g2+_plant",
                                               T ~ "g2+"),
                          a = if_else(mat_type == "g1g2", 23, 29),
                          m = if_else(mat_type == "g1g2", 40, 30),
                          mort = map(mort, ~ as.numeric(.$mort))) %>%  #from nested tibble to nested vector
                   #calculating the transition matrix for each day of each simulation, scenario and approach
                   mutate(elas_sens = future_pmap(list(A_trans = A,
                                                       mort = mort,
                                                       rob_test = .95,
                                                       a = a,
                                                       m = m), # arguments of the function calc_everything
                                                  calc_everything,
                                                  .options = furrr_options(chunk_size = 40)),
                          .keep = "unused")) #to remove extra columns that capture a lot of memory (e.g. trans_A)

plan(sequential) # to close the cluster infrastructure

save(future_elas_sens,
     file = paste0("data/mpm/sims_output/elas_sens_future_imp0_", Sys.Date(), ".RData"))
