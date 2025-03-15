# ----------------------------- #
# -- Supplementary figures --- #
# ---------------------------- #


# Packages -------

library(tidyverse)
library(ggtext)
library(MetBrewer)
library(patchwork)
library(gghalves)
library(ggdist)
library(distributional)
library(zoo)
library(broom)
library(glmmTMB)
library(popbio)
theme_set(theme_classic(base_family = "Arial"))


# FigS1: Mosaic of trends ------

## Data ------
cbms <- readxl::read_excel("data/cbms/annual_abundance.xlsx")
weeks <- readxl::read_excel("data/cbms/NSet_perItinYear.xlsx")

## Plots -----
trends <- cbms %>% 
  left_join(weeks) %>% 
  mutate(std_ab = SUMA_N_INDIV/N_SET) %>% 
  nest(df = -c(ID_ITIN)) %>% 
  mutate(lm = map(df,
                  ~ lm(std_ab ~ YEAR, data = .)),
         summary = map(lm, broom::tidy),
         summary = map(summary,
                       filter,
                       term == "YEAR")) %>% 
  unnest(summary) %>% 
  mutate(trend = case_when(p.value < 0.1 & estimate < 0 ~ "declining",
                           p.value < 0.1 & estimate > 0 ~ "resilient",
                           p.value > 0.1 ~ "non-significant",
                          T ~ NA_character_),
         n_year = map(df,
                      ~ n_distinct(.$YEAR)))


trends %>% 
  filter(!is.na(trend),
         n_year >= 10) %>% 
  group_by(trend) %>% 
  summarize(n_itin = n_distinct(ID_ITIN))


(mosaic <- trends %>% 
  filter(!is.na(trend),
         n_year >= 10) %>% 
  ggplot(aes(x = trend, y = estimate)) +
  geom_boxplot(aes(fill = trend)) +
  labs(y = "slope (weekly counts year<sup>-1</sup>)") +
  geom_hline(aes(yintercept = 0), color = "red") +
    geom_text(data = tibble(trend = c("declining", "non-significant", "resilient"),
                            n = c(4, 16, 4)),
              aes(x = trend,
                  y = 0.4,
                  label = paste("n =", n))) +
  scale_fill_viridis_d(option = "H", begin = .75, end = .25) +
  guides(fill = "none") +
    theme_classic() +
    theme(axis.title.y = element_markdown())) 

ggsave(filename = "figures/paper_figures/supp/figS1.png",
       plot = mosaic,
       device = png,
       height = 10,
       width = 10,
       units = "cm",
       dpi = 600)


# FigS2: Microclimatic conditions during spring and summer broods ------

## Data --------
micro <- read.csv2("data/clima/microclimate_data.csv") %>% 
  mutate(microhabitat = factor(microhabitat, levels = c("C", "OC", "O")),
         sensor = factor(sensor, levels = unique(.$sensor)))

## Plots ------
(micro_cond <- micro %>% 
    filter(hour >= 8,
           hour < 20,
           month > 3,
           month < 9, 
           microhabitat == "OC") %>% 
    mutate(site = factor(site, levels = c("AE", "CJ")),
           gen = case_when(site == "CJ" & month %in% 5:6 ~ "Spring",
                           site == "CJ" & month %in% 7:8 ~ "Summer",
                           site == "AE" & month %in% 4:5 ~ "Spring",
                           site == "AE" & yday >= 152 & yday < 196 ~ "Summer",
                           site == "AE" & yday >= 196 & yday < 244 ~ "Summer",
                           T ~ NA_character_)) %>%
    filter(!is.na(gen)) %>% 
    group_by(site, microhabitat, gen, sensor, year, yday) %>% 
    summarise(mean_diurnal_temp = mean(temp)) %>%  
    ggplot(aes(x = gen, y = mean_diurnal_temp)) +
    geom_boxplot(fill = "grey", outlier.shape = NA) +
    geom_hline(yintercept = c(20, 25), linetype = "dotted") +
    labs(x = "Season",
         y = "Mean diurnal temperature (ºC)") +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 4) +
    theme_classic())

ggsave(filename = "figures/paper_figures/supp/figS2.png",
       device = png,
       dpi = 600,
       height = 8,
       width = 8,
       units = "cm")



# FigS11: Validation -----

## Data -----
abund <- read.csv("data/cbms/cbms_nap_1_9_adultsurv.csv")

gr_real <- abund %>% 
  filter(esp == "piepna",
         anys >= 2014) %>% 
  mutate(period = case_when(mes < 5 ~ "g1g2",
                            mes < 7 ~ "g2g3",
                            (mes == 7 & itin == 1) | itin == 9 ~ "g3g4",
                            T ~ "g4g5")) %>% 
  group_by(itin, anys, period) %>% 
  summarise(abund = sum(ABUND, na.rm = T)) %>% 
  group_by(itin, anys) %>% 
  mutate(gr = lead(abund)/abund,
         site = if_else(itin == 1,
                        "Declining population",
                        "Resilient population")) %>% 
  filter(!is.na(gr)) %>% 
  rename(year = anys)

load("data/mpm/sims_output/sims_validation_imp0_2025-01-21.RData")

## Plots -----

(val_plots <- sims_validation_imp0 %>% 
  filter(year < 2018) %>% 
   mutate(plant = as.factor(plant),
          site = if_else(site == "lowl",
                         "Declining population",
                         "Resilient population")) %>% 
   split(.$site) %>% 
   map2(c(names(.)),
        ~ ggplot(data = .x,
                 aes(x = site, y = log(gr))) +
          geom_boxplot(aes(fill = plant), outlier.shape = NA,
                       position = position_dodge2(preserve = "single"),
                       show.legend = T) +
          geom_hline(aes(yintercept = 0), color = "red") +
          geom_point(data = gr_real %>%
                       filter(year < 2018,
                              site == .y),
                     aes(color = "Observed"),
                     # color = "black",
                     size = 10,
                     shape = "*") +
          facet_grid(. ~ period,
                     labeller = labeller(period = c('g1g2' = 'R<sub>0,1</sub>',
                                                    'g2g3' = 'R<sub>0,2</sub>',
                                                    'g3g4' = 'R<sub>0,3</sub>')),
                     space = "free",
                     scales = "free") +
          labs(title = .y,
               x = "Year",
               y = "ln R<sub>0</sub>") +
          scale_fill_manual(values = met.brewer("Kandinsky",
                                                n = 3)[c(1, 2)],
                            labels = c('TRUE' = "Projected with<br>plant scarcity",
                                       'FALSE' = "Projected with<br>enough plant"),
                            name = NULL,
                            drop = F) +
          scale_color_manual(values = "black",
                             name = NULL,
                             guide = guide_legend(override.aes = list(linetype = 0))) +
          theme_classic() +
          theme(axis.title.y = element_markdown(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank(),
                strip.text = element_markdown(),
                strip.background = element_blank(),
                legend.text = element_markdown())))


(val_together <- (val_plots$`Declining population` + 
    theme(legend.position = "none")) +
    (val_plots$`Resilient population` +
       guide_area() +
       plot_layout(ncol = 2,
                   guides = "collect",
                   width = c(2, 1))) +
    plot_layout(ncol = 1))


ggsave(filename = "figures/paper_figures/supp/FigS11.png",
       device = png,
       plot = val_together,
       width = 15,
       height = 15,
       units = "cm",
       dpi = 600)


# Table S3: Probabilities of extreme event occurrence --------
## Data ------
drought <- expand_grid(site = c("lowl", "mide"),
                       period = c("g1g2", "g2g3", "g3g4"),
                       GWL = c(1, 1.5, 2, 4)) %>% 
  filter(!(site == "mide" & period == "g3g4")) %>% 
  group_by(site, period) %>% 
  mutate(prob_occ = c(.2, .25, .4, .6),
         prob_occ = case_when(site == "mide" |
                                (site == "lowl" & period == "g1g2") ~ 0,
                              site == "lowl" & period == "g3g4" ~ 1,
                              T ~ prob_occ),
         prop_change = c(1, 1.25, 2, 3),
         prop_change = if_else(site == "lowl" & period == "g2g3",
                               prop_change,
                               1)) %>% 
  ungroup()

extreme_prob <- read.csv("data/clima/future_extremes_prob.csv") %>% 
  filter(Model == "CMIP6", !Bias_adjust, GWL != 3) %>% 
  select(TX, site, period, GWL, prob_occ, prop_change) %>% 
  add_row(drought)  %>% 
  mutate(TX = if_else(is.na(TX), "drought", as.character(TX)))

##Table ----
table <- extreme_prob %>% 
  select(-prop_change) %>% 
  pivot_wider(id_cols = c(TX, site, period),
              names_from = GWL,
              values_from = prob_occ) %>% 
  arrange(TX, site, period)

write.csv(table,
          file = "figures/paper_figures/supp/TableS3.txt",
          row.names = F)


# FigS10: IPCC projections --------------------
## data ------
historical_ipcc <- paste0("data/clima/interactive_atlas/monthly/",
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

ipcc <- paste0("data/clima/interactive_atlas/0",
               3:9, "/") %>% 
  map(~ list.files(path = ., pattern = "*.csv")) %>% 
  map(~ data.frame(file_name = .)) %>% 
  set_names(., nm = 3:9) %>%
  bind_rows(.id = "Month") %>% 
  mutate(file = paste0("data/clima/interactive_atlas/monthly/0",
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


## plot
(ipcc_projections <- ipcc %>% 
  filter(Model == "CMIP6",
         !Bias_adjust,
         GWL != 3) %>% 
  select(-starts_with("P")) %>% 
  mutate(prop_change = Median/Median_hist,
         prop_change = if_else(Median_hist == 0, NaN, prop_change)) %>% 
  pivot_longer(cols = starts_with("Med"),
               names_to = "prediction",
               values_to = "abs_month_freq") %>% 
  group_by(Month, Model, TX, Bias_adjust) %>% 
  mutate(dummy = prediction == "Median_hist",
         dummy = cumsum(dummy)) %>% 
  filter(prediction == "Median" | dummy == 1) %>% 
  mutate(GWL = if_else(dummy == 1 & prediction == "Median_hist", .85, as.numeric(GWL)),
         prop_change = if_else(GWL == .85, 1, prop_change)) %>% 
  select(-c(dummy, prediction)) %>% 
  pivot_longer(cols = c(prop_change, abs_month_freq),
               names_to = "var",
               values_to = "val") %>% 
  ggplot(aes(x = Month, y = val, color = as.factor(GWL))) +
  geom_line(aes(group = GWL)) +
  geom_point(size = 2) +
  scale_color_viridis_d(option = "F", begin = .8, end = .2) +
  guides(color = guide_legend(title = "GWL")) +
  scale_y_continuous(breaks = c(0, 5, 10)) +
  geom_text(data = data.frame(Month = 3, y = Inf, lab = c("(a)", "(b)", "(c)", "(d)"),
                              var = rep(c("abs_month_freq", "prop_change"), each = 2),
                              TX = rep(c(35, 40), times = 2)),
            aes(x = Month, y = y, label = lab),
            color = "black",
            hjust = 3, vjust = 2) +
  facet_grid(var ~ TX,
             scales = "free",
             labeller = labeller(var = c(prop_change = "Multiplicative factor<br>(N<sub>x,GWL</sub>/N<sub>x,0.85</sub>)",
                                         abs_month_freq = "Absolute frequency<br>(N<sub>x,GWL</sub>)"),
                                 TX = c('35' = "Extreme hot day<br>(T<sub>max</sub> > 35 °C)",
                                        '40' = "LLHI hot day<br>(T<sub>max</sub> > 40 °C)")),
             switch = "y") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        strip.text.y.left = element_markdown(),
        strip.text.x = element_markdown(),
        strip.background = element_blank(),
        strip.placement = "outside"))

ggsave(filename = "figures/paper_figures/supp/FigS10.png",
       device = png,
       plot = ipcc_projections,
       width = 15,
       height = 12,
       units = "cm",
       dpi = 600)



# FigS13: Macroclimate vs microclimate ------
## Data --------------------------------------------------------------------

## Projected R0's
load("data/mpm/sims_output/sims_future_drought_0.2-0.25-0.4-0.6_MACROCLIM_imp0_2025-01-22.RData")
load("data/mpm/sims_output/sims_future_drought_0.2-0.25-0.4-0.6_imp0_2024-12-10.RData")

## observed growth rates
abund <- read_csv("data/cbms/cbms_nap_1_9_adultsurv.csv")


brood_abund <- abund %>%
  mutate(brood = case_when(mes < 5 ~ 1,
                           mes < 7 ~ 2,
                           mes < 8 & itin == 1 ~ 3,
                           itin == 9 ~ 3,
                           TRUE ~ 4)) %>% 
  group_by(itin, esp, anys, brood) %>% 
  summarise(ABUND = sum(ABUND)) %>% 
  group_by(itin, esp, anys) %>% 
  mutate(next_abund = lead(ABUND),
         net_repr_rate = next_abund/ABUND,
         trend = if_else(net_repr_rate > 1, "inc", "dec"),
         period = case_when(brood == 1 ~ "g1g2",
                            brood == 2 ~ "g2g3",
                            brood == 3 ~ "g3g4",
                            T ~ NA_character_),
         site = if_else(itin == 1, "lowl", "mide")) %>% 
  filter(brood < 4,
         !(brood == 3 & itin == 9)) %>% 
  mutate(GWL = 0) %>% 
  ungroup() %>% 
  select(GWL, period, anys, site, net_repr_rate)

(macro_micro <- sims_imp0_macro %>% 
  select(sim_set, sim_cond, gr) %>% 
  bind_rows(select(sims_imp0, sim_set, sim_cond, gr), .id = "micro") %>% 
  mutate(micro = as.logical(as.numeric(micro)-1)) %>% 
  unnest(sim_cond) %>% 
  group_by(micro, sim_set) %>% 
  mutate(q1 = quantile(log(gr), probs = .25, na.rm = T),
         q3 = quantile(log(gr), probs = .75, na.rm = T),
         iqr = q3-q1,
         max_whisk = q3 + 1.5*iqr,
         min_whisk = q1 - 1.5*iqr,
         gr = if_else(log(gr) > max_whisk | log(gr) < min_whisk, NA_real_, gr)) %>% #filtering outliers 
  filter(!is.na(gr)) %>%
  bind_rows(rename(brood_abund, gr = net_repr_rate) %>% 
              mutate(GWL = 1),
            .id = "origin") %>% 
  mutate(origin = case_when(origin == "1" & micro ~ "Microclimatic<br>projections",
                            !micro ~ "Macroclimatic<br>projections",
                            T ~ "Observed"),
         origin = factor(origin,
                         levels = c("Observed",
                                    "Microclimatic<br>projections",
                                    "Macroclimatic<br>projections"))) %>% 
  filter(origin != "Observed") %>% 
  ggplot(aes(x = as.factor(GWL), y = log(gr), fill = origin)) +
  geom_boxplot(outlier.shape = NA,
               position = position_dodge(preserve = "single")) +
  facet_grid(site ~ period, 
             labeller = labeller(site = c("lowl" = "Declining population",
                                          "mide" = "Non-declining population"),
                                 period = c("g1g2" = "R<sub>0,1</sub>",
                                            "g2g3" = "R<sub>0,2</sub>",
                                            "g3g4" = "R<sub>0,3</sub>"))) +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(values = c(#met.brewer("Kandinsky", n = 3)[3],
                               met.brewer("Renoir")[c(5, 9)])) +
  labs(x = "Global warming level (°C)",
       y = "ln R<sub>0</sub>") +
  scale_x_discrete(labels = c("+0.85", "+1.5", "+2", "+4")) +
  theme_classic() +
  theme(legend.text = element_markdown(size = 11),
        legend.title = element_blank(),
        axis.title.y = element_markdown(),
        strip.text = element_markdown(),
        strip.background = element_blank()))
  

ggsave("figures/paper_figures/supp/FigS13.png",
       macro_micro,
       width = 18,
       height = 15,
       units = "cm",
       dpi = 600,
       device = png)


# jx and ex parameters ----
## Data -----------------
### Data of juveniles from the experiment with no plant scarcity
data_exp <- read.csv("data/mpm/growth_chamber_experiments.csv") %>% 
  filter(!is.na(total_cycle), experiment == "spring") %>% 
  dplyr::select(unique_id, total_cycle, Treatment, PtoA, LtoP, Egg_period, Larval_period) #LtoP to Larval_period needed for plant scarcity effects

### Data of juveniles from the experiment with plant scarcity
aut_exp <- read.csv("data/mpm/growth_chamber_experiments.csv") %>% 
  filter(experiment == "autumn", Treatment == "25") %>% 
  dplyr::select(unique_id, Egg_laying, LtoP, Egg_period, Larval_period)

## Parameters feeding Aspring -------
count_data <- data_exp %>% 
  filter(Treatment == 20) %>% 
  mutate(alive = 1,
         day = map(total_cycle, ~ data.frame(day = 1:.))) %>% 
  unnest(day) %>% 
  group_by(day) %>% 
  summarise(n = sum(alive))


mod_surv <- glm(n ~ poly(day, 13), family = "poisson", data = count_data)


ecl_data <- data_exp %>% 
  filter(Treatment == 20, PtoA)

mod_ecl <- glm(formula = (total_cycle-1)~1,
               family = "poisson",
               data = ecl_data)

fitted_juv_g1g2 <- count_data %>% 
  mutate(fitted = mod_surv$fitted.values,
         fitted_next = lead(fitted),
         surv = fitted_next/fitted,
         lambda = predict(mod_ecl,
                          newdata = slice_head(ecl_data),
                          type = "response"),
         ecl = dpois(day, lambda = lambda)) %>% 
  filter(!is.na(surv))


## Parameters feeding Asummer ----
count_data <- data_exp %>% 
  filter(Treatment == 25) %>% 
  mutate(alive = 1,
         day = map(total_cycle, ~ data.frame(day = 1:.))) %>% 
  unnest(day) %>% 
  group_by(day) %>% 
  summarise(n = sum(alive))

mod_surv <- glm(n ~ poly(day, 5), family = "poisson", data = count_data)

ecl_data <- data_exp %>% 
  filter(Treatment == 25, PtoA)

mod_ecl <- glm(formula = (total_cycle-1)~1,
               family = "poisson",
               data = ecl_data)

fitted_juv_g2g3 <- count_data %>% 
  mutate(fitted = mod_surv$fitted.values,
         fitted_next = lead(fitted),
         surv = fitted_next/fitted,
         lambda = predict(mod_ecl,
                          newdata = slice_head(ecl_data),
                          type = "response"),
         ecl = dpois(day, lambda = lambda)) %>% 
  filter(!is.na(surv))

## Parameters feeding Adrought --------
data_plant <- aut_exp %>%
  mutate(age_fatal = as.numeric(ymd("20151026")-as.Date(Egg_laying, tz = "UTC")),
         # what age individuals were when plant scarcity started
         plant_group = case_when(age_fatal < 8 ~ "w1", #it includes individuals that were age 6 and 7
                                 age_fatal < 16 ~ "w2", # ages 12, 13 and 15
                                 age_fatal < 24 ~ "w3", # ages 17, 18 and 19
                                 T ~  NA_character_),
         # at which week of dvt individuals were when plant scarcity started
         plant_group = factor(plant_group, levels = c("w1", "w2", "w3")),
         period = "plantscarc",
         alive = 1,
         total_cycle = Egg_period + Larval_period) %>%
  # only considering eggs and larvae, as there's no information on pupation for this experiment
  select(-c(Egg_laying, age_fatal)) %>% 
  filter(!is.na(plant_group), !is.na(total_cycle))

ref_data <- data_exp %>% 
  filter(Treatment == 25)


compar <- data.frame(plant_group = c("w1", "w2", "w3")) %>% 
  mutate(ref_df = map(plant_group,
                      ~ mutate(ref_data,
                               # the bootstrapped data for the reference matrix, repeated per each considered larval period
                               alive = 1,
                               period = "reference",
                               total_cycle = Egg_period + Larval_period,
                               # recalculating total cycle only with eggs and larvae
                               Treatment = NULL,
                               PtoA = NULL))) %>% 
  unnest(ref_df) %>% 
  bind_rows(data_plant) %>% 
  select(-c(Egg_period, Larval_period)) %>%  
  nest(df = -(period)) %>% 
  mutate(ecl = map(df, ~ distinct(., unique_id, LtoP, total_cycle)),
         # data for the estimation of larval pupation
         ecl = map(ecl, filter, LtoP),
         ecl = map(ecl, ~ mutate(., total_cycle = total_cycle - 1)),
         ecl = map(ecl, ~ glm(total_cycle ~ 1, family = poisson, data = .)),
         lambda = map(ecl, tidy), # mean day of larval pupation
         lambda = map(lambda, select, estimate),
         lambda = map(lambda, exp), # back-transformation from log distribution
         lambda = map(lambda, as.numeric),
         ecl = NULL) %>% 
  unnest(lambda) %>% 
  unnest(df) %>% 
  mutate(day = map(total_cycle, ~ data.frame(day = 1:.))) %>%
  # creating a dataset with an entry per individual and day
  unnest(day) %>% 
  mutate(ini = case_when(plant_group == "w1" ~ 6,
                         plant_group == "w2" ~ 12,
                         T ~ 17),
         end = case_when(plant_group == "w1" ~ 14,
                         plant_group == "w2" ~ 21,
                         T ~ 28)) %>%
  # initial and ending day of the different larval periods in which we will estimate larval survival
  # i.e. from the minimum larval age of the group to the end of the preceding week (after this first week we can't guarantee that plant scarcity persisted)
  filter(day >= ini,
         day <= end) %>%
  nest(count_data = -c(period, plant_group, lambda)) %>% 
  mutate(count_data = map(count_data, select, day, alive),
         count_data = map(count_data, group_by, day),
         count_data = map(count_data, summarise, n = sum(alive)), # survival curves (N ~ day)
         model = map(count_data, ~ glm(n ~ day, family = "poisson", data = .)),
         surv = map(model, tidy),
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


m <- nrow(fitted_juv_g2g3)
surv_rats <- c(rep(compar$ratio_surv[1], 13),
               rep(compar$ratio_surv[2], length(14:20)),
               rep(compar$ratio_surv[3], length(21:m)))


fitted_juv_plant <- fitted_juv_g2g3 %>%
  mutate(lambda = lambda-compar$dif_lambda[1],
         ecl = dpois(day, lambda = lambda),
         surv = surv*surv_rats,
         surv = if_else(surv > 1, 1, surv))

## Table S1------
jxex_param <- select(fitted_juv_g1g2, day, surv, ecl) %>% 
  left_join(select(fitted_juv_g2g3, day, surv, ecl), by = "day", suffix = c("_Aspring", "_Asummer")) %>% 
  left_join(select(fitted_juv_plant, day, surv, ecl), by = "day")

write_excel_csv2(jxex_param, file = "figures/paper_figures/supp/TableS1.txt", na = "0")

 
## FigS14: extreme events effect on elasticities ----------
## Data ----------------------
## Original matrices (i.e. without bootstrapping)
load("data/mpm/A_g1g2_figS2_imp0.RData")
load("data/mpm/A_g2g3_figS2_imp0.RData")
load("data/mpm/A_plant_figS2_imp0.RData")
source("code/simul/functions/change_juvs_mult.R")
### Original matrices with 5% extra mortality
mort_A <- list(A_g1g2_mort = A_g1g2, A_g2g3_mort = A_g2g3, A_plant_mort = A_plant) %>%
    map(~ change_juvs_mult(A = .,
                         mort = 0,
                         larvmort_robtest = 0.95)) %>%
  map(pluck, 1)

elas_calc <- function(A, m, a) {
  # Main eigen results from popbio package
  eigs <- eigen.analysis(A, zero = T)
  
  # Elasticities by matrix general component
  eigs$elas_comp <- matrix(nrow = 2, ncol = 2)
  eigs$elas_comp[1,1] <- apply(as.matrix(apply(eigs$elasticities[1:m, 1:m], 1, sum)), 2, sum)
  eigs$elas_comp[2,1] <- apply(as.matrix(apply(eigs$elasticities[(m+1):(m+a), 1:m], 1, sum)), 2, sum)
  eigs$elas_comp[1,2] <- apply(as.matrix(apply(eigs$elasticities[1:m, (m+1):(m+a)], 1, sum)), 2, sum)
  eigs$elas_comp[2,2] <- apply(as.matrix(apply(eigs$elasticities[(m+1):(m+a), (m+1):(m+a)], 1, sum)), 2, sum)
  
  return(eigs)
}


## Theoretical mortality in juveniles
all_mat_juvmort <- tibble(mat_type = names(mort_A),
                          mat = mort_A) %>% 
  mutate(mort = map(mat, ~seq(0, .99, by = .05)),
         a = if_else(str_detect(mat_type, "g1g2"), 23, 29),
         m = if_else(str_detect(mat_type, "g1g2"), 40, 30)) %>% 
  unnest(mort) %>% 
  mutate(mat = map2(mat, mort,
                    ~change_juvs_mult(A = .x, mort = .y)[[1]]),
         eigs = pmap(.l = list(A = mat, m = m, a = a),
                     .f = elas_calc)) %>% 
  unnest_wider(eigs)


## Theoretical mortality in adults
change_ad_mort <- function(A, mort) {
  
  # Dimensions of the transition matrix
  ## maximum juvenile age
  m <- length(which(A[1,] == 0))
  ## maximum adult age
  a <- ncol(A) - m
  
  # Sequence of transition matrices
  A_new <- list()
  
  for(j in seq_along(mort)) { #for each day of the projection
    A_new[[j]] <- A
    
    for (i in (m+2):(m+a)) { #then, add thermal mortality at adult parameters
      #surviving and resting as juveniles
      A_new[[j]][i, i-1] <- A[i, i-1]*(1-mort[j])
    }
  }
  
  return(A_new)
}


all_mat_admort <- tibble(mat_type = names(mort_A),
                         mat = mort_A) %>% 
  mutate(mort = map(mat, ~seq(0, .99, by = .05)),
         a = if_else(str_detect(mat_type, "g1g2"), 23, 29),
         m = if_else(str_detect(mat_type, "g1g2"), 40, 30)) %>% 
  unnest(mort) %>% 
  mutate(mat = map2(mat, mort,
                    ~change_ad_mort(A = .x, mort = .y)[[1]]),
         eigs = pmap(.l = list(A = mat, m = m, a = a),
                     .f = elas_calc)) %>% 
  unnest_wider(eigs)

## Theoretical reduction in adult fecundity

change_ad_fec <- function(A, factor) {
  # Dimensions of the transition matrix
  ## maximum juvenile age
  m <- length(which(A[1,] == 0))
  ## maximum adult age
  a <- ncol(A) - m
  
  # Sequence of transition matrices
  A_new <- list()
  
  for(j in seq_along(factor)) { #for each day of the projection
    A_new[[j]] <- A
    A_new[[j]][1, (m+1):(m+a)] <- A[1, (m+1):(m+a)]*(1-factor[j])
  }
  return(A_new)
}


all_mat_adfec <- tibble(mat_type = names(mort_A),
                        mat = mort_A) %>% 
  mutate(mort = map(mat, ~seq(0, .99, by = .05)),
         a = if_else(str_detect(mat_type, "g1g2"), 23, 29),
         m = if_else(str_detect(mat_type, "g1g2"), 40, 30)) %>% 
  unnest(mort) %>% 
  mutate(mat = map2(mat, mort,
                    ~change_ad_fec(A = .x, factor = .y)[[1]]),
         eigs = pmap(.l = list(A = mat, m = m, a = a),
                     .f = elas_calc)) %>% 
  unnest_wider(eigs)

## Plot -----------------
(extremes_on_elasticities <- all_mat_juvmort %>% 
   bind_rows(all_mat_admort, all_mat_adfec, .id = "affected_param") %>% 
   select(mat_type, mort, a, m, elas_comp, affected_param) %>% 
   mutate(elas_comp = map(elas_comp, as.tibble, rownames = "next_stage"),
          affected_param = case_when(affected_param == "1" ~ "juv",
                                     affected_param == "2" ~ "ad_surv",
                                     T ~ "fec"),
          mort = mort*100) %>% 
   unnest(elas_comp) %>% 
   pivot_longer(cols = V1:V2,
                names_to = "stage",
                names_prefix = "V",
                values_to = "elas") %>% 
   mutate(comp = case_when(stage == "1" & next_stage == "1" ~ "juv_surv",
                           stage == "1" ~ "juv_ecl",
                           next_stage == "1" ~ "ad_fec",
                           T ~ "ad_surv")) %>% 
   ggplot(aes(x = mort, y = elas, color = comp)) +
   geom_smooth(se = F, alpha = .1) +
   # geom_point() +
   # geom_line() +
   scale_color_discrete(labels = c("Fecundity", "Adult survival",
                                   "Juvenile transition", "Juvenile stasis"),
                        name = "Matrix<br>component") +
   facet_grid(. ~ affected_param,
              labeller = labeller(
                affected_param = c("ad_surv" = "Adult survival",
                                   "fec" = "Fecundity",
                                   "juv" = "Juvenile survival"))) +
   labs(x = "Extreme event impact<br>(percentual reduction in the transition parameter)",
        y = "Elasticity") +
   ylim(c(0,1)) +
   theme_classic() +
   theme(strip.background = element_blank(),
         legend.title = element_markdown(hjust = .5),
         legend.position = "bottom",
         axis.title.x = element_markdown()))

ggsave("figures/paper_figures/supp/FigS14.png",
       extremes_on_elasticities,
       width = 18,
       height = 12,
       units = "cm",
       dpi = 600,
       device = png)


## FigS3: Plots experiments ----------------------

(count_plot <- data_exp %>% 
  filter(Treatment == 20) %>% 
  mutate(alive = 1,
         day = map(total_cycle, ~ data.frame(day = 1:.))) %>% 
  unnest(day) %>% 
  group_by(day) %>% 
  summarise(n = sum(alive)) %>% 
  ggplot(aes(x = day, y = n)) +
  geom_point(aes(color = '20 ºC'), size = 1.5) +
  geom_smooth(method = "glm", formula = y~poly(x,13), se = T,
              method.args = list(family = "poisson"),
              color = "black", fill = "darkorchid") +
  geom_point(data = data_exp %>% filter(Treatment == 25) %>% 
              mutate(alive = 1,
                     day = map(total_cycle, ~ data.frame(day = 1:.))) %>% 
              unnest(day) %>% 
              group_by(day) %>% 
              summarise(n = sum(alive)),
            aes(color = '25 ºC'),
            size = 1.5) +
  scale_color_manual(values = c("darkorchid", "forestgreen"),
                     name = "Treatment") +
  geom_smooth(data = data_exp %>% filter(Treatment == 25) %>% 
                mutate(alive = 1,
                       day = map(total_cycle, ~ data.frame(day = 1:.))) %>% 
                unnest(day) %>% 
                group_by(day) %>% 
                summarise(n = sum(alive)),
              method = "glm",
              formula = y ~ poly(x, 5),
              method.args = list(family = "poisson"),
              se = T,
              color = "black",
              fill = "forestgreen") +
  labs(x = "Age in days (x)",
       y = "Number of juveniles (N<sub>x</sub>)") +
  theme(axis.title.y = element_markdown()))


lambda_se <- data_exp %>% 
  filter(Treatment != 15, PtoA) %>% 
  nest(df = -Treatment) %>% 
  mutate(mod = map(df, ~glm(total_cycle-1 ~ 1, data = .)),
         lambda = map2(mod, df,
                       ~predict(.x, newdata = slice_head(.y),
                                type = "response", se.fit = T)),
         lambda = map(lambda, flatten)) %>% 
  unnest_wider(lambda) %>%
  select(-c(mod, df)) %>% 
  rename('lambda' = "1") 


(ecl_plot <- data_exp %>% 
  filter(Treatment != 15, PtoA) %>% 
  ggplot(aes(x = as.factor(Treatment), y = total_cycle-1)) +
  geom_dotsinterval(aes(fill = as.factor(Treatment),
                        side = "both"),
                    slab_alpha = .5,
                    slab_linewidth = 0,
                    slab_color = NA) +
  geom_pointrange(data = lambda_se, aes(y = lambda,
                                        ymax = lambda+se.fit,
                                        ymin = lambda-se.fit,
                                        color = "predicted"),
                  size = .5) +
  scale_fill_manual(values = c("darkorchid", "forestgreen"),
                    name = "Treatment") +
  scale_color_manual(values = "black", name = NULL) +
  guides(fill = "none") +
  labs(x = "Treatment (ºC)",
       y = "Age of pupal eclosion (days)"))

(juv_exp_plot <- count_plot + ecl_plot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")"))

ggsave("figures/paper_figures/supp/FigS3.png",
       juv_exp_plot,
       width = 24,
       height = 12,
       units = "cm",
       device = png,
       dpi = 600)


## FigS8: Plant scarcity -------
### data ---------------
data_plant <- aut_exp %>%
  mutate(age_fatal = as.numeric(ymd("20151026")-as.Date(Egg_laying, tz = "UTC")),
         # what age individuals were when plant scarcity started
         plant_group = case_when(age_fatal < 8 ~ "w1", #it includes individuals that were age 6 and 7
                                 age_fatal < 16 ~ "w2", # ages 12, 13 and 15
                                 age_fatal < 24 ~ "w3", # ages 17, 18 and 19
                                 T ~  NA_character_),
         # at which week of dvt individuals were when plant scarcity started
         plant_group = factor(plant_group, levels = c("w1", "w2", "w3")),
         period = "plantscarc",
         alive = 1,
         total_cycle = Egg_period + Larval_period) %>%
  # only considering eggs and larvae, as there's no information on pupation for this experiment
  select(-c(Egg_laying, age_fatal)) %>% 
  filter(!is.na(plant_group), !is.na(total_cycle))

ref_data <- data_exp %>% 
  filter(Treatment == 25)


surv_compar <- data.frame(plant_group = c("w1", "w2", "w3")) %>% 
  mutate(ref_df = map(plant_group,
                      ~ mutate(ref_data,
                               # the bootstrapped data for the reference matrix, repeated per each considered larval period
                               alive = 1,
                               period = "reference",
                               total_cycle = Egg_period + Larval_period,
                               # recalculating total cycle only with eggs and larvae
                               Treatment = NULL,
                               PtoA = NULL))) %>% 
  unnest(ref_df) %>% 
  bind_rows(data_plant) %>% 
  select(-c(Egg_period, Larval_period)) %>% 
  mutate(day = map(total_cycle, ~ data.frame(day = 1:.))) %>%
  # creating a dataset with an entry per individual and day
  unnest(day) %>% 
  mutate(ini = case_when(plant_group == "w1" ~ 6,
                         plant_group == "w2" ~ 12,
                         T ~ 17),
         end = case_when(plant_group == "w1" ~ 14,
                         plant_group == "w2" ~ 21,
                         T ~ 28)) %>%
  # initial and ending day of the different larval periods in which we will estimate larval survival
  # i.e. from the minimum larval age of the group to the end of the preceding week (after this first week we can't guarantee that plant scarcity persisted)
  filter(day >= ini,
         day <= end) %>%
  nest(count_data = -c(period, plant_group)) %>% 
  # nest(count_data = -c(period, plant_group, lambda)) %>% 
  mutate(count_data = map(count_data, select, day, alive),
         count_data = map(count_data, group_by, day),
         count_data = map(count_data, summarise, n = sum(alive))) %>% 
  unnest(count_data)


ecl_compar <- ref_data %>% 
  mutate(total_cycle = Egg_period + Larval_period,
         period = "reference",
         Treatment = NULL,
         PtoA = NULL) %>% 
  bind_rows(data_plant) %>% 
  select(-c(Egg_period, Larval_period)) %>%  
  nest(df = -(period)) %>% 
  mutate(ecl = map(df, ~ distinct(., unique_id, LtoP, total_cycle)),
          # data for the estimation of larval pupation
         ecl = map(ecl, filter, LtoP),
          ecl = map(ecl, ~ mutate(., total_cycle = total_cycle - 1)),
          mod = map(ecl, ~ glm(total_cycle ~ 1, family = poisson, data = .)),
          lambda = map(mod, tidy), # mean day of larval pupation
          lambda = map(lambda, select, estimate, std.error),
           lambda = map(lambda, exp), # back-transformation from log distribution
           lambda = map(lambda, as.numeric)) %>% 
  mutate(Treatment = if_else(period == "reference", "Fresh<br>host plant", "Scarce<br>host plant"))
  

### plot ------
(surv_plot <- surv_compar %>% 
  mutate(period = factor(period,
                         levels = c("reference", "plantscarc"))) %>% 
  ggplot(aes(x = day, y = n, color = period)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "glm",
              method.args = list(family = "poisson"),
              se = T,
              color = "black",
              aes(fill = period)) +
  facet_wrap(. ~ plant_group,
             labeller = labeller(plant_group = c('w1' = "1-week old individuals<br>when host-plant<br>scarcity started",
                                            'w2' = "2-weeks old individuals<br>when host-plant<br>scarcity started",
                                            'w3' = "3-weeks old individuals<br>when host-plant<br>scarcity started"))) +
  scale_color_manual(values = met.brewer("Isfahan1", n = 7)[c(5, 3)],
                     labels = c('plantscarc' = "Scarce<br>host plant",
                                'reference' = "Fresh<br>host plant"),
                     name = "Treatment",
                     aesthetics = c("fill", "color")) +
  labs(x = "Age in days (x)",
       y = "Number of juveniles (N<sub>x</sub>)") +
  theme(strip.text = element_markdown(),
        strip.background = element_blank(),
        axis.title.y = element_markdown(),
        legend.text = element_markdown()))

(ecl_plot <- ecl_compar %>% 
    unnest(ecl) %>% 
    ggplot(aes(x = Treatment, y = total_cycle)) +
    geom_dotsinterval(aes(fill = as.factor(Treatment),
                          side = "both"),
                      slab_alpha = .5,
                      slab_linewidth = 0,
                      slab_color = NA) +
    geom_pointrange(data = ecl_compar %>% 
                      unnest_wider(lambda, names_sep = "_"),
                    aes(y = lambda_1,
                        ymax = lambda_1+lambda_2,
                        ymin = lambda_1-lambda_2,
                        color = "predicted"),
                    size = .5) +
    scale_fill_manual(values = met.brewer("Isfahan1", n = 7)[c(5, 3)]) +
    scale_color_manual(values = "black", name = NULL) +
    guides(fill = "none") +
    labs(x = "Treatment",
         y = "Age of pupation (days)") +
    theme(axis.text.x = element_markdown()))

(plant_exp_plot <- surv_plot + ecl_plot +
    plot_layout(guides = "collect",
                ncol = 1) +
    plot_annotation(tag_levels = "a",
                    tag_prefix = "(",
                    tag_suffix = ")"))

ggsave("figures/paper_figures/supp/figS8.png",
       plant_exp_plot,
       width = 18,
       height = 18,
       units = "cm",
       device = png,
       dpi = 600)

## FigS4: Juv survival vs predation rate -------------
(pred_factor <- list("g1g2" = fitted_juv_g1g2,
                     "g2g3" = fitted_juv_g2g3,
                     "plant" = fitted_juv_plant) %>% 
   bind_rows(.id = "scenario") %>% 
   expand_grid(pred = (0:10)/100) %>% 
   group_by(scenario, pred) %>% 
   mutate(surv = surv,
          surv = surv*(1-pred), # predatory pressure
          ecl = if_else(day < 8, 0, ecl),
          #imposing 0 probability of eclosion for juveniles < 8 days old
          surv_grow = surv*(1-ecl),
          surv_emer = surv * ecl,
          cum_survgrow = cumprod(surv_grow),
          cum_survgrow = lag(cum_survgrow, default = 1),
          prob_ad = surv_emer * cum_survgrow) %>% 
   summarise(prob_ad = sum(prob_ad)))

ad_surv_probs <- tibble(source = c("Yamamoto 1981 - G1", "Yamamoto 1981 - G2",
                                   "Forsberg 1987 - a", "Forsberg 1987 - b"),
                        prob = c(.14, .03, .1, .064),
                        pred = c(.1, 0, .1, .1),
                        hjust = c(1, 0, 1, 1))

(mort_factor_plot <- ggplot(data = pred_factor, aes(x = pred, y = prob_ad, color = scenario)) +
    geom_point() +
    geom_line() +
    geom_hline(data = ad_surv_probs,
               aes(yintercept = prob), linetype = "dotted") +
    geom_text(data = ad_surv_probs,
              aes(x = pred, y = prob + .005, label = source, color = NULL,
                  hjust = hjust), color = "black") +
    scale_color_manual(values = met.brewer("Isfahan1", n = 7)[c(7, 5, 3)],
                       name = "Matrix type",
                       labels = c("A<sub>spring</sub>",
                                  "A<sub>summer</sub>",
                                  "A<sub>drought</sub>")) +
    scale_x_continuous(breaks = seq(.01, .09, by = .02)) +
    labs(x = "Daily predation",
         y = "Total juvenile survival\n(Probability for an egg to become an adult)") +
    theme_bw(base_family = "Arial") +
    theme(legend.text = element_markdown()))

ggsave(filename = "figures/paper_figures/supp/FigS4.png",
       mort_factor_plot,
       height = 15,
       width = 18,
       units = "cm",
       dpi = 600)


# FigS7: Thermal mortality vs TX -------------------------------------------------
### data --------------------------------------------------------------------
mortalities <- read.csv("data/clima/future_mortalities.csv")

### plot --------------------------------------------------------------------
(thermalmort <- mortalities %>% 
   mutate(mort = mort/100,
          tmax_group = factor(tmax_group,
                              levels = c("non-extreme", "35", "40"))) %>% 
   filter(tmax > 25) %>%
   ggplot(aes(x = tmax, y = mort)) +
   geom_point(aes(color = tmax_group)) +
   geom_smooth(method = "glm",
               method.args = list(family = binomial),
               se = F,
               linewidth = .5,
               color = "grey30") +
   scale_color_manual(values = c("black",
                                 met.brewer(name = "Hokusai1", n = 7)[c(3, 2)]),
                      labels = c('35' = "non-LLHI hot day<br>(35 < T<sub>max</sub> < 40 °C)",
                                 '40' = "LLHI hot day<br>(T<sub>max</sub> > 40 °C)",
                                 'non-extreme' = "Mild day<br>(T<sub>max</sub> < 35 °C)"),
                      name = NULL) +
   labs(x = "Daily micro T<sub>max</sub> (ºC)",
        y = "Thermal mortality (<i>k</i>)") +
   theme(axis.title.x = element_markdown(),
         axis.title.y = element_markdown(angle = 90, vjust = .5),
         legend.text = element_markdown()))

ggsave(filename = "figures/paper_figures/supp/FigS7.png",
       plot = thermalmort,
       height = 12,
       width = 18,
       units = "cm")


#fx and ax parameters ----
## Data of adult fecundity (from published papers) ----
data_art <- read.csv("data/mpm/mean_daily_fec_articles.csv") %>% 
  filter(!is.na(dvt)) %>% 
  mutate(daily_eggs = round(daily_eggs/2))

## Data of adult survival (P. oleracea, from Fig. S2 of Kerr 2020)
data_adsurv <- read.csv("data/mpm/ad_survival_oleracea.csv", row.names = 1) %>%
  arrange(day)

## Parameters feeding Aspring -----
mod_repr <- data_art %>% 
  filter(dvt == "diap") %>% 
  glmmTMB(formula = daily_eggs ~ poly(day, 3) + 
            (1 | article_fig) +
            (1 | mat_system),
          family = "nbinom2",
          data = .)

predicted_repr <- predict(mod_repr,
                          newdata = data.frame(day = 1:23,
                                               mat_system = NA,
      # P                                         article_fig = NA),
                          type = "response")

fitted_repr_g1g2 <- data.frame(day = 1:23,
                               eggs = round(predicted_repr))

fitted_ad_g1g2 <- fitted_repr_g1g2 %>% 
  left_join(data_adsurv)

## Parameters feeding Asummer and Adrought -------
mod_repr <- data_art %>% 
  filter(dvt == "direct") %>% 
  glmmTMB(formula = daily_eggs ~ poly(day, 5) + 
            (1 | article_fig) +
            (1 | mat_system),
          family = "nbinom2",
          data = .)

predicted_repr <- predict(mod_repr,
                          newdata = data.frame(day = 1:29,
                                               mat_system = NA,
                                               article_fig = NA),
                          type = "response")

fitted_repr_g2g3 <- data.frame(day = 1:29,
                               eggs = round(predicted_repr))

fitted_ad <- fitted_repr_g2g3 %>% 
  left_join(data_adsurv) %>% 
  left_join(fitted_ad_g1g2, by = "day", suffix = c("Asummer", "Aspring"))

## Table S2 -------------------
write_excel_csv2(fitted_ad, file = "figures/paper_figures/supp/fxax_param.txt", na = "0")


## FigS5: Plots fecundity ------------
(fec_plot <- data_art %>% 
  filter(dvt == "diap") %>% 
  ggplot(aes(x = day, y = daily_eggs)) +
  geom_point(aes(color = "diapausing")) +
  geom_smooth(method = "glm",
              formula = y ~ poly(x, 3),
              method.args = list(family = "poisson"),
              color = "darkorchid",
              fill = "darkorchid") +
  geom_point(data = data_art %>% 
               filter(dvt == "direct"),
             aes(color = "directly<br>developing")) +
  geom_smooth(data = data_art %>% 
                filter(dvt == "direct"),
              method = "glm",
              formula = y ~ poly(x, 5),
              method.args = list(family = "poisson"),
              color = "forestgreen",
              fill = "forestgreen") +
  scale_color_manual(values = c("darkorchid", "forestgreen"),
                     name = "Female<br>group") +
  labs(x = "Days since adult emergence (x)",
       y = "Number of eggs laid") +
  theme(legend.text = element_markdown(),
        legend.title = element_markdown()))

ggsave("figures/paper_figures/supp/FigS5.png",
       fec_plot,
       width = 15,
       height = 12,
       units = "cm",
       device = png,
       dpi = 600)


# FigS6: A matrices parameters ------
## Data ----------------------
load("data/mpm/sims_output/sims_future_drought_0.2-0.25-0.4-0.6_imp0_2024-12-10.RData")
source("code/simul/functions/change_juvs_mult.R")

matrices <- sims_imp0 %>%
  select(sim_set:dims, -sim_df) %>% 
  unnest(sim_cond) %>% 
  filter(GWL == 1) %>% 
  mutate(A = map(A,
                 ~ change_juvs_mult(A = .,
                                    mort = 0,
                                    larvmort_robtest = 0.95)),
         A = map(A, pluck, 1),
         A = map(A, as_tibble),
         A = map(A,
                 ~ mutate(., "next_stage" = colnames(.))),
         A = map(A,
                 ~ pivot_longer(.,
                               cols = -next_stage,
                               names_to = "stage",
                               values_to = "value"))) %>%
  unnest(A) %>% 
  filter(value != 0 | next_stage == "a01") %>%
  mutate(treat = case_when(str_detect(period, "g1g2") ~ "spring",
                           sim_with_plant == "0" ~ "summer",
                           T ~ "dry summer"),
         day = str_sub(stage, start = 2, end = -1L),
         day = factor(day, levels = c(paste0("0",
                                             as.character(1:9)),
                                      as.character(10:50))),
         next_day = str_sub(next_stage, start = 2, end = -1L),
         next_day = as.numeric(next_day),
         parameter = case_when(next_stage == "j01" ~ "repr",
                               next_stage == "a01" ~ "ecl",
                               str_sub(stage,
                                       start = 1,
                                       end = 1) == "j" ~ "juv_growth",
                               T ~ "ad_growth"))
## Plot ----------------
labs_reduced <- map_chr(1:40,
                        ~if_else((.-1)%%3!=0, "", as.character(.)))

(mat_params <- matrices %>% 
  mutate(parameter = factor(parameter, levels = c("juv_growth", "repr",
                                                  "ecl", "ad_growth")),
         treat = factor(treat,
                         levels = c("spring", "summer", "dry summer"))) %>% 
  ggplot(aes(x = day,
                     y = value)) +
  geom_boxplot(outlier.shape = NA,
               aes(color = treat, fill = treat),
               alpha = .6,
               position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = met.brewer("Isfahan1", n = 7)[c(7, 5, 3)],
                    aesthetics = c("fill", "color"),
                    name = "Matrix type",
                    labels = c("Spring", "Summer", "Dry summer")) +
  labs(x = "Age within stage (days)",
       y = "Transition parameter") +
  scale_x_discrete(labels = labs_reduced) +
  facet_wrap(vars(parameter),
             nrow = 2, scales = "free",
             labeller = labeller(parameter = c(juv_growth = "<b>A<sub>JJ</sub></b>",
                                               ecl = "<b>A<sub>JM</sub></b>",
                                               repr = "<b>A<sub>MJ</sub></b>",
                                               ad_growth = "<b>A<sub>MM</sub></b>"))) +
    theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_markdown(size = 11)))

ggsave("figures/paper_figures/supp/FigS6.png",
       mat_params,
       width = 24,
       height = 18,
       units = "cm",
       device = png,
       dpi = 600)

# FigS15: R0 against fecundity -----------------
## Data --------------------------------------------------------------------

# Simulation results
load("data/mpm/sims_output/sims_future_drought_0.2-0.25-0.4-0.6_imp0_2024-12-10.RData")


## Plot --------------------------------------------------------------------
fecs <- sims_imp0 %>% 
  select(sim_set:dims, gr) %>% 
  mutate(fec_df = map2(A, dims,
                       ~ data.frame(fec = .x[1, (.y[1]+1):(.y[1]+.y[2])],
                                    ad_surv = c(1,
                                                .x[(.y[1]+2):
                                                     (.y[1]+.y[2]),
                                                   (.y[1]+1):
                                                     (.y[1]+.y[2])-1]
                                                [which(.x[(.y[1]+2):
                                                            (.y[1]+.y[2]),
                                                          (.y[1]+1):
                                                            (.y[1]+.y[2])-1] != 0)]))),
         A = NULL,
         dims = NULL) %>% 
  unnest(c(fec_df, sim_cond)) %>% 
  group_by(sim_set, site, period, GWL, sim_with_plant, sim) %>% 
  mutate(cum_surv = cumprod(ad_surv),
         mean_fec = fec*cum_surv) %>% 
  summarise(max_fec = sum(fec),
            mean_fec = sum(mean_fec),
            gr = mean(gr))


(fecs_gr <- fecs %>%
    mutate(mat = case_when(period == "g1g2" ~ "Spring",
                           sim_with_plant == 0 ~ "Summer",
                           T ~ "Dry summer"),
           mat = factor(mat, levels = c("Spring", "Summer", "Dry summer"))) %>% 
    ggplot(aes(x = mean_fec, y = gr,
               color = mat,
               fill = mat)) +
    geom_point(alpha = 0.05, size = 0.4, aes(color = mat)) +
    geom_hline(aes(yintercept = 1), linetype = "dashed") +
    geom_smooth(method = "lm") +
    # facet_wrap(vars(GWL), nrow = 2) +
    facet_grid(site ~ period, labeller = labeller(site = c("lowl" = "Declining population",
                                                           "mide" = "Non-declining population"),
                                                  period = c("g1g2" = "R<sub>0,1</sub>",
                                                             "g2g3" = "R<sub>0,2</sub>",
                                                             "g3g4" = "R<sub>0,3</sub>"))) +
    scale_color_manual(values = met.brewer("Isfahan1", n = 7)[c(7, 5, 3)],
                       name = "Matrix type",
                       aesthetics = c("fill", "colour"),
                       guide = guide_legend(override.aes = list(fill = NA))) +
    labs(y = "Net reproductive rate (R<sub>0</sub>)",
         x = "Realised fecundity (total eggs laid by a female)") +
    theme_classic() +
    theme(strip.text = element_markdown(size = 11),
          strip.background = element_blank(),
          axis.title.y = element_markdown()))

ggsave(filename = "figures/paper_figures/fecundity/FigS15.png",
       device = png,
       plot = fecs_gr,
       width = 21,
       height = 15,
       units = "cm",
       dpi = 600)

(fecs_sc <- ggplot(data = fecs, aes(y = mean_fec,
                                    x = scenario)) +
    geom_point(alpha = 0.05, position = "jitter",
               aes(color = scenario)) +
    geom_boxplot(outlier.shape = NA,
                 aes(fill = scenario)))

ggsave(filename = "figures/paper_figures/FigS15.png",
       plot = fecs_sc,
       width = 15,
       height = 15,
       units = "cm",
       dpi = 600)


# FigS16: Predator pressure -------
load("data/mpm/sims_output/boot/sims_future_drought_0.2-0.25-0.4-0.6_pred_0_imp0_2025-01-23.RData")
sims_pred <- sims_future 
  
load("data/mpm/sims_output/boot/sims_future_drought_0.2-0.25-0.4-0.6_pred_4_imp0_2025-01-24.RData")
sims_pred <- sims_pred %>% 
  bind_rows(sims_future)

load("data/mpm/sims_output/boot/sims_future_drought_0.2-0.25-0.4-0.6_pred_8_imp0_2025-01-25.RData")
sims_pred <- sims_pred %>% 
  bind_rows(sims_future)

load("data/mpm/sims_output/boot/sims_future_drought_0.2-0.25-0.4-0.6_pred_7_imp0_2025-01-25.RData")
sims_pred <- sims_pred %>% 
  bind_rows(sims_future)

load("data/mpm/sims_output/boot/sims_future_drought_0.2-0.25-0.4-0.6_imp0_2024-12-10.RData")

sims_pred <- sims_imp0 %>% 
  select(-A) %>% 
  unnest(sim_cond) %>% 
  mutate(pred_factor = .95) %>% 
  bind_rows(select(sims_pred, -A) %>% 
              unnest(sim_cond))


(pred_pressure <- sims_pred %>% 
  group_by(pred_factor, sim_set) %>% 
  mutate(q1 = quantile(log(gr), probs = .25, na.rm = T),
         q3 = quantile(log(gr), probs = .75, na.rm = T),
         iqr = q3-q1,
         max_whisk = q3 + 1.5*iqr,
         min_whisk = q1 - 1.5*iqr,
         gr = if_else(log(gr) > max_whisk | log(gr) < min_whisk, NA_real_, gr)) %>% #filtering outliers 
  filter(!is.na(gr)) %>% 
  ggplot(aes(x = as.factor(GWL), y = log(gr), fill = as.factor(round((1-pred_factor)*100)))) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(site ~ period,
             labeller = labeller(site = c("lowl" = "Declining population",
                                          "mide" = "Non-declining population"),
                                 period = c("g1g2" = "R<sub>0,1</sub>",
                                            "g2g3" = "R<sub>0,2</sub>",
                                            "g3g4" = "R<sub>0,3</sub>"))) +
  scale_fill_manual(values = met.brewer("Renoir", n = 7),
                     name = "Daily<br>predation (%)") +
  labs(x = "Global warming level (°C)",
       y = "ln R<sub>0</sub>") +
  scale_x_discrete(labels = c("+0.85", "+1.5", "+2", "+4")) +
  theme_classic() +
  theme(legend.title = element_markdown(),
        legend.box.just = "center",
        axis.title.y = element_markdown(),
        strip.text = element_markdown(),
        strip.background = element_blank()))


ggsave("figures/paper_figures/supp/FigS16.png",
       plot = pred_pressure,
       device = png,
       height = 18,
       width = 21,
       units = "cm",
       dpi = 600)


# FigS17: Longevity variation ----

load("data/mpm/sims_output/boot/sims_future_drought_diflong_imp0_2025-01-25.RData")

(long_test <- sims_future %>% 
  select(gr, sim_cond, sim, dims) %>% 
  unnest(c(sim_cond, dims)) %>% 
  mutate(long = rep(c("m", "a"), 20000),
         long_num = case_when(period == "g1g2" & long == "m" ~ 40,
                              period == "g1g2" & long == "a" ~ 23,
                              long == "m" ~ 30,
                              T ~ 29),
         mean_ecl = case_when(period == "g1g2" & long == "m" ~ 29,
                              period == "g2g3" & long == "m" ~ 22,
                              long == "m" ~ 20,
                              T ~ NA_integer_)) %>% 
  ggplot(aes(x = dims, y = log(gr))) +
  geom_point(aes(color = as.factor(GWL)), alpha = .01) +
  geom_smooth(aes(color = as.factor(GWL)), alpha = .5, se = F) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = long_num, linetype = "longevity<br>used"),) +
  scale_linetype_manual(values = c("dotted", "dotted"),
                        name = NULL) +
  scale_color_viridis_d(option = "F", begin = .8, end = .2,
                        name = "GWL (°C)") +
  facet_grid(long ~ site + period,
             scales = "free",
             labeller = labeller(long = c("m" = "Juvenile longevitiy",
                                          "a" = "Adult longevity"),
                                 site = c("lowl" = "Declining<br>population",
                                          "mide" = "Non-declining<br>population"),
                                 period = c("g1g2" = "R<sub>0,1</sub>",
                                            "g2g3" = "R<sub>0,2</sub>",
                                            "g3g4" = "R<sub>0,3</sub>"))) +
  # facet_wrap(~long + site + period, scales = "free", ncol = 5) +
  ylim(-6,NA) +
  scale_x_continuous(breaks = seq(from = 25, to = 45, by = 10)) +
  labs(x = "Stage maximum longevity (days)",
       y = "ln R<sub>0</sub>") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text= element_markdown(),
        axis.title.y = element_markdown(),
        legend.text = element_markdown()))


ggsave("figures/paper_figures/supp/FigS17.png",
       plot = long_test,
       device = png,
       height = 15,
       width = 21,
       units = "cm",
       dpi = 600)



# FigS18: 0-eclosion test ------------
load("data/mpm/sims_output/boot/sims_future_drought_0.2-0.25-0.4-0.6_imp0_2024-12-10.RData")
load("data/mpm/sims_output/boot/sims_future_drought_0.2-0.25-0.4-0.6_2023-07-08.RData")

## observed growth rates
abund <- read_csv("data/cbms/cbms_naprap_1_9_adultsurv.csv")


brood_abund <- abund %>%
  mutate(brood = case_when(mes < 5 ~ 1,
                           mes < 7 ~ 2,
                           mes < 8 & itin == 1 ~ 3,
                           itin == 9 ~ 3,
                           TRUE ~ 4)) %>% 
  group_by(itin, esp, anys, brood) %>% 
  summarise(ABUND = sum(ABUND)) %>% 
  group_by(itin, esp, anys) %>% 
  mutate(next_abund = lead(ABUND),
         net_repr_rate = next_abund/ABUND,
         trend = if_else(net_repr_rate > 1, "inc", "dec"),
         period = case_when(brood == 1 ~ "g1g2",
                            brood == 2 ~ "g2g3",
                            brood == 3 ~ "g3g4",
                            T ~ NA_character_),
         site = if_else(itin == 1, "lowl", "mide")) %>% 
  filter(esp == "piepna",
         brood < 4,
         !(brood == 3 & itin == 9),
         anys > 1988) %>% 
  mutate(GWL = 0) %>% 
  ungroup() %>% 
  select(GWL, period, anys, site, net_repr_rate)

(compar_plot <- sims_future %>% 
   select(sim_set, sim_cond, gr) %>% 
   bind_rows(select(sims_imp0, sim_set, sim_cond, gr), .id = "imp0") %>% 
   mutate(imp0 = as.logical(as.numeric(imp0)-1)) %>% 
   unnest(sim_cond) %>% 
   group_by(imp0, sim_set) %>%
    mutate(q1 = quantile(log(gr), probs = .25, na.rm = T),
           q3 = quantile(log(gr), probs = .75, na.rm = T),
           iqr = q3-q1,
           max_whisk = q3 + 1.5*iqr,
           min_whisk = q1 - 1.5*iqr,
           gr = if_else(log(gr) > max_whisk | log(gr) < min_whisk, NA_real_, gr)) %>% #filtering outliers 
   filter(!is.na(gr)) %>%
   bind_rows(rename(brood_abund, gr = net_repr_rate) %>% 
               mutate(GWL = 1),
             .id = "origin") %>% 
   mutate(origin = case_when(origin == "1" & imp0 ~ "Minimum<br>eclosion at 8<br>(followed approach)",
                             !imp0 ~ "Without<br>minimum<br>eclosion age",
                             T ~ "Observed"),
          origin = factor(origin,
                          levels = c("Observed",
                                     "Without<br>minimum<br>eclosion age",
                                     "Minimum<br>eclosion at 8<br>(followed approach)"))) %>% 
   filter(origin != "Observed") %>% 
   ggplot(aes(x = as.factor(GWL), y = log(gr), fill = origin)) +
   geom_boxplot(outlier.shape = NA,
                position = position_dodge(preserve = "single")) +
   facet_grid(site ~ period, 
              labeller = labeller(site = c("lowl" = "Declining population",
                                           "mide" = "Non-declining population"),
                                  period = c("g1g2" = "R<sub>0,1</sub>",
                                             "g2g3" = "R<sub>0,2</sub>",
                                             "g3g4" = "R<sub>0,3</sub>"))) +
   geom_hline(aes(yintercept = 0)) +
   scale_fill_manual(values = c(
     met.brewer("Renoir")[c(9, 5)])) +
   labs(x = "Global warming level (°C)",
        y = "ln R<sub>0</sub>") +
   scale_x_discrete(labels = c("+0.85", "+1.5", "+2", "+4")) +
   theme_classic() +
   theme(legend.text = element_markdown(size = 11),
         legend.title = element_blank(),
         axis.title.y = element_markdown(),
         strip.text = element_markdown(),
         strip.background = element_blank()))

ggsave(filename = "figures/paper_figures/supp/FigS18.png",
       device = png,
       plot = compar_plot,
       width = 21,
       height = 18,
       dpi = 600,
       units = "cm")
