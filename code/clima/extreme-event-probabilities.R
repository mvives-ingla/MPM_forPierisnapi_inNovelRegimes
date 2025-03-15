# Estimation of current and future probabilities of heat extremes

# Packages ------
library(tidyverse)
library(lubridate)


# 1. Current probabilities of heat extremes in microclimatic records -----
## Data -------
### Microclimatic data recorded in semi-open microsites where P. napi lays eggs
micro <- read.csv2("data/clima/microclimate_data.csv")

## Estimation ------
probs_micro <- micro %>% 
  mutate(period = case_when(site == "CJ" & month %in% 5:6 ~ "g1g2",
                            site == "CJ" & month %in% 7:8 ~ "g2g3",
                            site == "AE" & month %in% 4:5 ~ "g1g2",
                            site == "AE" & yday >= 152 & yday < 196 ~ "g2g3",
                            site == "AE" & yday >= 196 & yday < 244 ~ "g3g4",
                            T ~ NA_character_)) %>% 
  filter(microhabitat %in% c("OC"),
         !is.na(period)) %>% 
  group_by(sensor, site, microhabitat, date, year, period, yday) %>% 
  summarise(tmax = max(temp, na.rm = T),
            max_35 = tmax >= 35,
            max_40 = tmax >= 40) %>% 
  pivot_longer(cols = max_35:max_40,
               names_to = "TX",
               names_prefix = "max_",
               values_to = "thres_surpassed") %>% 
  nest(df = - TX) %>% 
  mutate(mod = map(df,
                   ~ glm(thres_surpassed ~ site*period,
                         data = .,
                         family = "binomial")),
         prob_occ = map(mod,
                        ~ expand_grid(site = c("CJ", "AE"),
                                      period = c("g1g2", "g2g3", "g3g4"))),
         prob_occ = map2(prob_occ, mod,
                         ~ mutate(.x,
                                  prob_occ = predict(.y,
                                                     newdata = .x,
                                                     type = "response",
                                                     re.form = NA)))) %>% 
  unnest(prob_occ) %>% 
  filter(!(period == "g3g4" & site == "CJ")) %>% 
  select(-(df:mod))


write.csv(probs_micro,
          "data/clima/actual_extreme_prob.csv")


# 2. Future probabilities of extreme events based on IPCC scenarios -----
## More info in: https://interactive-atlas.ipcc.ch/ 

## Macroclimatic probabilities ------
### Baseline scenario -----
historical_ipcc <- paste0("data/clima/interactive_atlas/",
                          paste0(0, 3:9), "/historical/") %>% 
  map(~ list.files(path = ., pattern = "*.csv")) %>% 
  map(~ data.frame(file_name = .)) %>% 
  set_names(., nm = 3:9) %>%
  bind_rows(.id = "Month") %>% 
  mutate(file = paste0("data/clima/interactive_atlas/0",
                       Month,
                       "/historical/",
                       file_name),
         file = map(file, read.csv)) %>% 
  unnest(file) %>% 
  filter(str_detect(Period, "AR6")) %>% 
  separate(col = file_name,
           into = c("Model", "variable", "period", "season"),
           sep = " - ") %>% 
  mutate(TX = str_extract(variable, "\\d+"),
         TX = as.double(TX),
         Bias_adjust = str_detect(variable, "Bias"),
         Period = str_sub(Period, start = -10, end = -2)) %>% 
  select(-c(season, period, Scenario, variable)) %>% 
  rename_with(.fn = ~ paste(., "hist", sep = "_"),
              .cols = Median:P95)

### Future scenario ----
ipcc <- paste0("data/clima/interactive_atlas/0",
               3:9, "/") %>% 
  map(~ list.files(path = ., pattern = "*.csv")) %>% 
  map(~ data.frame(file_name = .)) %>% 
  set_names(., nm = 3:9) %>%
  bind_rows(.id = "Month") %>% 
  mutate(file = paste0("data/clima/interactive_atlas/0",
                       Month,
                       "/",
                       file_name),
         file = map(file, read.csv)) %>% 
  unnest(file) %>% 
  filter(!is.na(Median)) %>% 
  rename(GWL = Period) %>% 
  mutate(file_name = str_sub(file_name, start = 1, end = -34)) %>% 
  separate(col = file_name,
           into = c("Model", "variable", "period", "season"),
           sep = " - ") %>% 
  mutate(TX = str_extract(variable, "\\d+"),
         TX = as.double(TX),
         Bias_adjust = str_detect(variable, "Bias"),
         Period = "1995-2014",
         period = NULL,
         GWL = str_replace(GWL, "Warming ", ""),
         GWL = str_sub(GWL, 1, -3)) %>% 
  select(!(Median:P95), -c(variable, season, Scenario), Median:P95) %>% 
  left_join(historical_ipcc)


## Microclimatic probabilities -----
probs_micro <- read_csv("data/clima/actual_extreme_prob.csv") %>% 
  select(2:5)

probs_extreme_future <- tibble(site = c(rep("CJ", 4), rep("AE", 6)),
                               period = rep(c("g1g2", "g2g3", "g1g2", "g2g3", "g3g4"), each = 2),
                               Month = as.character(c(5:8, 4:7, 7:8))) %>%
  # months included in each generation period
  left_join(ipcc) %>% 
  group_by(site, period, Model, Bias_adjust, GWL, TX) %>% 
  summarise(across(Median:P95_hist, sum)) %>% 
  mutate(prop_change = Median/Median_hist) %>%  # increase factor
  select(site, period, Model, Bias_adjust, GWL, TX,
         Median, prop_change) %>% 
  left_join(probs_micro) %>% 
  rename(actual_prob_occ = prob_occ) %>% 
  mutate(future_prob_occ = actual_prob_occ*prop_change) %>% 
  pivot_longer(cols = actual_prob_occ:future_prob_occ,
               names_to = "trial",
               names_pattern = "(.*)_prob_occ",
               values_to = "prob_occ") %>% 
  # we want to add actual prob occ as a row more with GWL = 1 (or 0.85)
  # but when pivoting, the row of GWL = 1 is repeated for each of the other projected GWL
  mutate(dummy = trial == "actual") %>%
  group_by(site, period, TX, Model, Bias_adjust, trial) %>% 
  mutate(dummy = cumsum(dummy)) %>% 
  filter(dummy <= 1) %>% 
  mutate(GWL = if_else(dummy == 1, 1, as.numeric(GWL)),
         prop_change = if_else(dummy == 1, 1,
                               as.numeric(prop_change)),
         site = if_else(site == "AE", "lowl", "mide")) %>% 
  ungroup() %>% 
  arrange(site, period, Model, Bias_adjust, TX, GWL) %>% 
  select(site, period, Model, Bias_adjust, TX, GWL, prob_occ, prop_change) 


write.csv(probs_extreme_future,
          file = "data/clima/future_extremes_prob.csv")





