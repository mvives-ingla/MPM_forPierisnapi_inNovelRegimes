# Sensitivity and elasticity analysis of the bootstrapped matrices

# Packages ----------------------------------------------------------------

library(tidyverse)
library(popdemo)
library(MetBrewer)
theme_set(theme_classic(base_family = "Arial"))


# Data --------------------------------------------------------------------

## Elasticity and sensitivity calculated from bootstrapped matrices
load("data/mpm/sims_output/elas_sens_future_imp0_2025-01-20.RData")

## Original matrices (i.e. without bootstrapping)
load("data/mpm/A_g1g2_figS2_imp0.RData")
load("data/mpm/A_g2g3_figS2_imp0.RData")
load("data/mpm/A_plant_figS2_imp0.RData")

### Original matrices with 5% extra mortality (used predation rate)
mort_A <- list(A_g1g2_mort = A_g1g2, A_g2g3_mort = A_g2g3, A_plant_mort = A_plant) %>%
  map(~ change_juvs_mult(A = .,
                         mort = 0,
                         larvmort_robtest = 0.95)) %>%
  map(pluck, 1)




# Plot --------------------------------------------------------------------
## Calculation of the elasticities and sensitivities of true matrices ----
true_elas_sens <- tibble(mat_type = c("g1g2", "g2+", "g2+_plant",
                                      "g1g2", "g2+", "g2+_plant"),
                         mort_5 = c(F, F, F, T, T, T),
                         trans_A = list(A_g1g2, A_g2g3, A_plant) %>%
                           append(mort_A)) %>%
  mutate(Sensitivity = map(trans_A, sens),
         Elasticity = map(trans_A, elas),
         a = if_else(mat_type == "g1g2", 23, 29),
         m = if_else(mat_type == "g1g2", 40, 30),
         dims = map2(m, a, ~ list(m = .x, a = .y)),
         trans_A = NULL,
         Sensitivity = map2(Sensitivity, dims,
                            ~ matrix(.x,
                                     nrow = .y$m + .y$a,
                                     ncol = .y$m + .y$a,
                                     dimnames = list(c(paste0("j0", 1:9),
                                                       paste0("j", 10:.y$m),
                                                       paste0("a0", 1:9),
                                                       paste0("a", 10:.y$a)),
                                                     c(paste0("j0", 1:9),
                                                       paste0("j", 10:.y$m),
                                                       paste0("a0", 1:9),
                                                       paste0("a", 10:.y$a)))))) %>%
  pivot_longer(cols = Sensitivity:Elasticity,
               names_to = "sens_elas_type",
               values_to = "sens_elas_value") %>%
  mutate(sens_elas_value = map(sens_elas_value, as_tibble),
         sens_elas_value = map2(sens_elas_value, dims,
                                ~ mutate(.x,
                                         next_stage = factor(c(paste0("j0", 1:9),
                                                               paste0("j", 10:.y$m),
                                                               paste0("a0", 1:9),
                                                               paste0("a", 10:.y$a)),
                                                             levels = rev(c(paste0("j0", 1:9),
                                                                            paste0("j", 10:.y$m),
                                                                            paste0("a0", 1:9),
                                                                            paste0("a", 10:.y$a)))))),
         sens_elas_value = map(sens_elas_value,
                               ~ pivot_longer(.,
                                              cols = -next_stage,
                                              names_to = "stage",
                                              values_to = "value")),
         sens_elas_value = map2(sens_elas_value, dims,
                                ~ mutate(.x,
                                         stage = factor(stage,
                                                        levels = c(paste0("j0", 1:9),
                                                                   paste0("j", 10:.y$m),
                                                                   paste0("a0", 1:9),
                                                                   paste0("a", 10:.y$a)))))) %>%
  unnest(sens_elas_value) %>%
  mutate(main_stage = str_sub(stage, start = 1, end = 1),
         main_next_stage = str_sub(next_stage, start = 1, end = 1),
         component = case_when(main_stage == "j" & main_next_stage == "j" ~ "juv_surv",
                               next_stage == "a01" ~ "juv_ecl",
                               next_stage == "j01" ~ "fec",
                               T ~ "ad_surv"),
         component = factor(component, levels = c("juv_surv", "juv_ecl", "fec", "ad_surv"))) %>%
  group_by(mat_type, mort_5, sens_elas_type, component) %>%
  summarise(value = sum(value)) %>%
  ungroup()



## Sensitivity and elasticity by matrix component

(elasens_bootreal_rob <- future_elas_sens %>%
    unnest(elas_sens) %>% 
    mutate(sens_elas_type = if_else(sens_elas_type == "sens", "Sensitivity", "Elasticity")) %>%
    filter(sens_elas_type == "Elasticity") %>% 
    group_by(mat_type, sim_set, sim, sens_elas_type, component) %>% 
    summarise(value = median(value, na.rm = T)) %>% 
    ggplot(aes(x = component, y = value,
               fill = mat_type,
               color = mat_type)) +
    geom_boxplot(outlier.shape = NA,
                 alpha = .6) +
    ylim(0, 1) +
   scale_fill_manual(values = met.brewer("Isfahan1", n = 7)[c(7, 5, 3)],
                     name = "Matrix type",
                     labels = c("Spring",
                                "Summer",
                                "Dry summer"),
                     aesthetics = c("color", "fill")) +
    geom_point(data = true_elas_sens %>%
                 filter(mort_5, sens_elas_type == "Elasticity"),
               aes(shape = "Real value"),
               position = position_dodge2(width = .75),
               size = 1,
               alpha = 1,
               color = "black") +
    scale_shape_manual(name = NULL,
                       values = c("Real value" = 8)) +
    guides(fill = guide_legend(order = 1,
                               override.aes = list(shape = NA)),
           color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    scale_x_discrete(labels = c("Juvenile\nstasis",
                                "Juvenile\ntransition",
                                "Adult\nfecundity",
                                "Adult\nsurvival")) +
    labs(x = "Component",
         y = "Elasticity") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     vjust = 1,
                                     size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 10),
          legend.text = element_text(size  = 10),
          # legend.position = "top",
          legend.title = element_text(size = 10)))

ggsave(filename = "figures/paper_figures/Figure4.png",
       plot = elasens_bootreal_rob,
       width = 14,
       height = 12,
       units = "cm",
       device = png,
       dpi = 600)






