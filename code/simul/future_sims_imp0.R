#! /usr/bin/env Rscript
#Simulations imposing a minimum age for pupal eclosion

# Packages
library(tidyverse)
library(glmmTMB)
library(popdemo)
library(zoo)
library(furrr)
library(lubridate)
library(broom)


plan(multicore,  workers = 40)


# Functions
## To build transition matrices given survival and eclosion probabilities and fecundities
source("code/simul/functions/build_matrix_imp0.R")
## To simulate a transition matrix by bootstraping from survival, eclosion and fecundity data
source("code/simul/functions/simul_matrix.R")
## To change juvenile mortality and eclosion due to plant scarcity effects
source("code/simul/functions/plant_scarcity_imp0.R")
## To estimate thermal mortality depending on maximum temperature
source("code/simul/functions/thermal_mortality.R")
## To change the parameters of juvenile mortality in the transition matrix given an external mortality pulse
source("code/simul/functions/change_juvs_mult.R")
## To project the dynamics of a population given a transition matrix and scenario
source("code/simul/functions/project_pop_alt.R")
## To perform several simulations of population dynamics given particular conditions
source("code/simul/functions/simul_scenario.R")



# data
## Data of juveniles from the experiment with no plant scarcity
data_exp <- read.csv("data/mpm/growth_chamber_experiments.csv", row.names = 1) %>% 
  filter(!is.na(total_cycle), experiment == "spring") %>%
  dplyr::select(unique_id, total_cycle, Treatment, PtoA, LtoP, Egg_period, Larval_period) #LtoP to Larval_period needed for plant scarcity effects

## Data of juveniles from the experiment with plant scarcity
aut_exp <- read.csv("data/mpm/growth_chamber_experiments.csv", row.names = 1) %>% 
  filter(experiment == "autumn", Treatment == "25") %>% 
  dplyr::select(unique_id, Egg_laying, LtoP, Egg_period, Larval_period)

## Data of adult fecundity (from published papers)
data_art <- read.csv("data/mpm/mean_daily_fec_articles.csv") %>% 
  filter(!is.na(dvt)) %>% 
  mutate(daily_eggs = round(daily_eggs/2)) #just eggs that will become females


## Data of adult survival (P. oleracea, from Fig. S2 of Kerr 2020)
data_adsurv <- read.csv("data/mpm/ad_survival_oleracea.csv", row.names = 1) %>%
  arrange(day)

data_list <- list(juv = data_exp,
                  repr = data_art,
                  ad = data_adsurv,
                  plsc = aut_exp)

## Probability that a thermal extreme event occurs at a MICROCLIMATIC scale
prob_extreme <- read_csv("data/clima/future_extremes_prob.csv") %>%
  filter(GWL != 3,
         Model == "CMIP6",
         !Bias_adjust) %>%
  pivot_wider(id_cols = c(site, period, GWL),
              names_from = TX,
              names_prefix = "prob",
              values_from = prob_occ)

## Drought events probabilities
prob_drought <- data.frame(prob_drought = c(.2, .25, .4, .6),
                           GWL = c(1, 1.5, 2, 4))



sims_imp0 <- expand_grid(site = c("lowl", "mide"),
                           period = c("g1g2", "g2g3", "g3g4"),
                           GWL = c(1, 1.5, 2, 4)) %>%
  mutate(sim_set = paste(site, period, GWL, sep = "-")) %>% 
  filter(!(site == "mide" & period == "g3g4")) %>% 
  left_join(prob_extreme) %>% 
  left_join(prob_drought) %>% 
  mutate(prob_drought = case_when(site == "mide" | period == "g1g2" ~ 0,
                                  period == "g3g4" ~ 1,
                                  T ~ prob_drought)) %>% 
  nest(sim_cond = -sim_set) %>% 
  mutate(simul = map(sim_cond,
                     ~ simul_scenario(nsims = 10000,
                                      ndays = 70,
                                      broods = .$period,
                                      data = data_list,
                                      prob_35 = .$prob35,
                                      prob_40 = .$prob40,
                                      prob_drought = .$prob_drought,
                                      rob_test = .95))) %>% 
  unnest(simul)




save(sims_imp0, file = paste0("data/mpm/sims_output/sims_future_drought_",
                                str_c(prob_drought$prob_drought, collapse = "-"),
                                "_imp0_", Sys.Date(), ".RData"))


