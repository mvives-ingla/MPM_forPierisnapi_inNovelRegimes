#! /usr/bin/env Rscript

# Simulation of the model over true microclimatic series


# packages ----------------------------------------------------------------
library(tidyverse)
library(glmmTMB)
library(popdemo)
library(zoo)
library(lubridate)
library(furrr)

plan(multicore,  workers = 50)


# data preparation --------------------------------------------------------------------

# Microclimatic data until 2020
micro <- read.csv2("data/clima/microclimate_data.csv") %>% 
  mutate(microhabitat = factor(microhabitat, levels = c("C", "OC", "O")),
         sensor = factor(sensor, levels = unique(.$sensor)),
         period = case_when(site == "CJ" & month %in% 5:6 ~ "g1g2",
                            site == "CJ" & month %in% 7:8 ~ "g2g3",
                            site == "AE" & month %in% 4:5 ~ "g1g2",
                            site == "AE" & yday >= 152 & yday < 196 ~ "g2g3",
                            site == "AE" & yday >= 196 & yday < 244 ~ "g3g4",
                            T ~ NA_character_),
         site = if_else(site == "CJ","mide", "lowl")) %>% 
  filter(microhabitat == "OC", !is.na(period), year <= 2018) %>% 
  group_by(sensor, site, year, period, yday) %>% 
  summarise(tx = max(temp, na.rm = T)) %>% 
  group_by(sensor, site, year, period) %>% 
  mutate(min_day = min(yday, na.rm = T),
         max_day = max(yday, na.rm = T),
         ndays = n_distinct(yday),
         needed_days = case_when(period == "g1g2" ~ 61,
                                 site == "mide" ~ 62,
                                 period == "g2g3" ~ 44,
                                 T ~ 48)) %>% 
  filter(needed_days == ndays) %>% #being sure that complete periods are taken
  select(sensor:tx) %>% 
  mutate(tx35 = if_else(tx >= 35, 1, 0),
         tx40 = if_else(tx >= 40, 1, 0))


## Data of juveniles from the experiment with no plant scarcity
data_exp <- read.csv("data/mpm/growth_chamber_experiments.csv") %>% 
  filter(!is.na(total_cycle), experiment == "spring") %>% 
  dplyr::select(unique_id, total_cycle, Treatment, PtoA, LtoP, Egg_period, Larval_period) #LtoP to Larval_period needed for plant scarcity effects

## Data of juveniles from the experiment with plant scarcity
aut_exp <- read.csv("data/mpm/growth_chamber_experiments.csv") %>% 
  filter(experiment == "autumn", Treatment == "25") %>% 
  dplyr::select(unique_id, Egg_laying, LtoP, Egg_period, Larval_period)

## Data of adult fecundity (from published papers)
data_art <- read.csv("data/mpm/mean_daily_fec_articles.csv") %>% 
  filter(!is.na(dvt)) %>% 
  mutate(daily_eggs = round(daily_eggs/2))


## Data of adult survival (P. oleracea, from Fig. S2 of Kerr 2020)
data_adsurv <- read.csv("data/mpm/ad_survival_oleracea.csv", row.names = 1) %>%
  arrange(day)

data_list <- list(juv = data_exp,
                  repr = data_art,
                  ad = data_adsurv,
                  plsc = aut_exp)

treat <- c("g1g2_lowl_F", "g2g3_lowl_F", "g2g3_lowl_T", "g3g4_lowl_T","g1g2_mide_F",
           "g2g3_mide_F") # scenarios to simulate
period <- c("g1g2", "g2g3", "g3g4")
# periods simulated (from brood 1 to brood 2, 2 to 3 & 3 to 4)
site <- c("lowl", "mide") # sites simulated (lowland vs mid-elevation)
plant <- c("T", "F") # whether the scenario simulates plant scarcity (T) or not (F)




# Functions ---------------------------------------------------------------

thermal_mortality <- function(group, tx) {
  
  termor <- as.numeric()
  
  for (i in seq_along(group)) {
    if(group[i] == "00") {
      termor[i] <- exp(-23.1705859770024+0.501511610449519*tx[i])
    } else if (group[i] == "10") {
      termor[i] <- exp(-22.6557794795911+0.484215947292227*tx[i])
    } else if (group[i] == "11") {
      termor[i] <- 1/(1+exp(-(-37.9697676066161+0.873108327565261*tx[i])))
    }
  }
  
  return(termor)
  
}


source("code/simul/functions/build_matrix_imp0.R")
source("code/simul/functions/plant_scarcity_imp0.R")
source("code/simul/functions/project_pop_alt.R")
source("code/simul/functions/simul_matrix.R")



change_juvs_mult <- function(A, mort, larvmort_robtest = 1) {

  # Dimensions of the transition matrix
  ## maximum juvenile age
  m <- length(which(A[1,] == 0))
  ## maximum adult age
  a <- ncol(A) - m
  
  # Sequence of transition matrices
  A_new <- list()
  
  
  for(j in seq_along(mort)) { #for each day of the projection
    
    A_new[[j]] <- A
      
   for (i in 2:(m+a)) { #then, add thermal mortality at juvenile parameters
      if(i <= m) {
        #surviving and resting as juveniles
        A_new[[j]][i, i-1] <- A[i, i-1]*(1-mort[j])*larvmort_robtest
        
      } 
      if(i == (m+1)) {
        A_new[[j]][i, 1:(m)] <- A[i, 1:(m)]*(1-mort[j])*larvmort_robtest
      }
    }
  }
  
  return(A_new)
}



# Calculations ----------------------------------------------------------------

sims_validation_imp0 <- expand_grid(sim = 1:500,
                              site = site,
                              period = period,
                              plant = plant) %>%
  mutate(scenario = paste(period, site, plant, sep = "_"),
         plant = as.logical(plant)) %>% 
  filter(scenario %in% treat) %>% 
  left_join(micro) %>% # adding thermal series per scenario ordered by day
  group_by(scenario, site, period, sensor, sim, year) %>% 
  mutate(days = row_number(),
         mort = thermal_mortality(group = paste0(tx35, tx40), tx = tx)) %>%  #the intensity of the mortality event
  ungroup() %>% 
  nest(sim_df = c(yday, days, tx, tx35, tx40, mort)) %>% 
  # in line with what we observed for R2 in 2017 (plant scarcity from day 156: 5/6)
  mutate(A = future_map(period,
                 ~ simul_matrix(data_juv = data_exp,
                                data_repr = data_art,
                                data_ad = data_adsurv,
                                brood_period = .)),
         # we create a transition matrix per each simulation and return it with the associated bootstrapped samples
         # future_map should be applied to an ungrouped df, as now!
         dims = map(A, pluck, "dims"), #dimensions of the transition matrix (m & a)
         A = future_map_if(.x = A,
                           .p = plant,
                           ~ plant_scarcity(A = .x,
                                        data_plant = aut_exp)), # a bootstrapped transition matrix if plant scarcity occurs
         A = map_if(A, .p = !plant, pluck, "A"),
         A_seq = future_map2(sim_df, A,
                      ~ change_juvs_mult(A = .y,
                                         mort = .x$mort,
                                         larvmort_robtest = .95)),
         # creating a sequence of daily transition matrices depending on thermal mortality and plant scarcity
         growth_rate = future_map2(A_seq, dims,
                            ~ project_pop(
                              A = .x,
                              Aseq = 1:length(.x),
                              vector = c(rep(0, .y["m"]), 1, rep(0, .y["a"] - 1))),
                            .options = furrr_options(chunk_size = 30)),
         # projection of the population with the given sequence of transition matrices
         A_seq = NULL,
         growth_rate = future_map(growth_rate,
                           ~ calculate_growth(vec_transenv = .,
                                              start_name = "a01"))) %>%
  # calculation of the net reproductive rate of the population
  unnest(growth_rate)

save(sims_validation_imp0, file = paste0("data/mpm/sims_output/boot/sims_validation_imp0_",
                                           Sys.Date(),
                                           ".RData"))