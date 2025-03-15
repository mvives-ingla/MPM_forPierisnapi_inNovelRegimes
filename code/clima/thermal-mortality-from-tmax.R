# Code to obtain the model equations of predicted daily thermal mortality
# against the maximum microclimatic temperature recorded during that day(tmax)

# Packages -----
library(tidyverse)


# 1. Predicted thermal mortality with the TDT model -----
## Data -------
### Microclimatic data recorded in semi-open microsites where P. napi lays eggs
micro <- read.csv2("data/clima/microclimate_data.csv")

### Survival time in static assays of heat tolerance with 89 reared larvae of P. napi
#### more info on: https://doi.org/10.1002/ecm.1561
tdt <- read.table("data/clima/tdt_experiments.txt", header = T)

## Functions ------------------
### Functions for calculating the Thermal Death Time curves to thermal survival experiments
### and apply the dynamic model that predicts daily thermal mortality by integrating the effects of all 1-min thermal exposures across the day
### more info on: https://doi.org/10.5061/DRYAD.STQJQ2C1R
### and: https://doi.org/10.5281/ZENODO.7358091
source("code/clima/Thermal_landscape_functions_mod.R")


## 1.1 Build tolerance landscape from static assays of heat tolerance -----
tl <- tolerance.landscape(ta = tdt$SENSOR_mean_temp,
                          time = tdt$aprox_minute_dead)

## 1.2 Predict thermal survival in a sample of 640 days -----
selected_days <- micro %>% 
  filter(month >= 6, month < 9) %>% 
  group_by(sensor, site, microhabitat, date, year, month, yday) %>% 
  summarise(tmax = max(temp, na.rm = T)) %>% 
  mutate(group = case_when(tmax >= 40 ~ "40",
                           tmax >= 35 ~ "35",
                           T ~ "normal")) %>% 
  ungroup() %>% 
  nest(df = -group) %>% 
  mutate(df = map_if(.x = df,
                     .p = group == "normal",
                     .f = ~ slice_sample(.x, n = 240))) %>% 
  unnest(df) %>% 
  select(sensor:date) %>% 
  left_join(micro) %>% 
  group_by(sensor, date, hour) %>% 
  mutate(same_hour = n(),
         dup = row_number(time)) %>% 
  filter(same_hour != 2 |
           (same_hour == 2 & dup == 2)) %>% 
  group_by(sensor, date) %>% 
  mutate(count = n()) %>% 
  filter(count == 24) %>% 
  mutate(sens_date = paste(sensor, date, sep = "_")) %>% 
  split(.$sens_date) %>%
  map(.,
      ~ spline(.$temp, n = (nrow(.)-1)*60)) %>% 
  map(bind_cols) %>%
  map(rename,
      hour = x,
      ta.min = y) %>% 
  map2(names(.),
       ~ mutate(.x,
                sens = str_sub(.y, 1, -11),
                date = str_sub(.y, -10, -1))) %>% 
  map(~ dynamic.landscape.mod(ta = .$ta.min,
                              tolerance.landscape = tl,
                              sens = .$sens,
                              sp = "PN",
                              day = .$date,
                              plot = F))

## 1.3 Prepare a database of predicted thermal mortalities with TDT dynamical model ----
mortalities <- selected_days %>% 
  map(mutate,
      tmax = max(ta)) %>% 
  map(slice_min,
      alive) %>% 
  bind_rows(.id = "sensor_data") %>% 
  mutate(data = str_sub(sensor_data, -10, -1),
         sensor = str_sub(sensor_data, 1, -11),
         sensor_data = NULL,
         mort = 100-alive,
         tmax_group = case_when(tmax >= 40 ~ "40",
                                tmax >= 35 ~ "35",
                                T ~ "non-extreme"))

write.csv(x = mortalities,
          file = ".data/clima/future_mortalities.csv")

# 2. Logistic function to obtain predicted thermal mortality from daily tmax -----
## Data -----
### Obtained thermal mortalities in section 1
mortalities <- read.csv("data/clima/future_mortalities.csv")

## Binomial & logistic fit between predicted thermal mortality and tmax ----
mort_mod <- mortalities %>% 
  mutate(mort = mort/100) %>% 
  nest(df = -tmax_group) %>% 
  mutate(mod_logit = map(df,
                         ~ glm(mort ~ tmax,
                               family = "binomial",
                               data = .)),
         mod_exp = map(df,
                       ~ lm(log(mort) ~ tmax,
                            data = .))) %>% 
  pivot_longer(cols = mod_logit:mod_exp,
               names_to = "model_type",
               values_to = "model") %>% 
  filter((tmax_group == 40 & model_type == "mod_logit") |
           (tmax_group != 40 & model_type != "mod_logit")) %>% 
  mutate(coefs = map(model, broom::tidy)) %>% 
  unnest(coefs) %>% 
  pivot_wider(id_cols = c(tmax_group, model_type, df),
              names_from = term,
              values_from = estimate) %>% 
  rename(int = "(Intercept)") %>% 
  mutate(form = if_else(model_type == "mod_logit",
                        paste0("y ~ 1/(1+exp(-(",
                               int,
                               "+",
                               tmax,
                               "*tmax)))"),
                        paste0("y ~ exp(",
                               int,
                               "+",
                               tmax,
                               "*tmax)")))




mortalities %>% 
  filter(tmax_group != "40") %>% 
  mutate(tmax_group = if_else(tmax_group == "non-extreme", "<35",
                              as.character(tmax_group))) %>% 
  ggplot(aes(x = tmax, y = mort/100)) +
  geom_point(aes(color = as.factor(sensor)), alpha = .5, size = 2) +
  geom_smooth(method = "lm",
              formula = y ~ exp(x)) +
  facet_wrap(vars(tmax_group),
             scales = "free",
             ncol = 1)


mortalities %>% 
  filter(tmax_group == "40") %>% 
  mutate(tmax_group = if_else(tmax_group == "non-extreme", "<35",
                              as.character(tmax_group))) %>% 
  ggplot(aes(x = tmax, y = mort/100)) +
  geom_point(aes(color = as.factor(sensor)), alpha = .5, size = 2) +
  geom_smooth(method = "glm",
              method.args = list(family = binomial),
              se = F) +
  facet_wrap(vars(tmax_group),
             scales = "free",
             ncol = 1)

