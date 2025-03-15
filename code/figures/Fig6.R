# Separate and compound effects of the stressors


# Packages ----------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(ggtext)
library(zoo)
library(MetBrewer)
library(furrr)
library(popdemo)
library(ggnewscale)
library(grid)
theme_set(theme_classic(base_family = "Arial", base_size = 9))


# Data --------------------------------------------------------------------

load("data/mpm/sims_output/sims_future_drought_0.2-0.25-0.4-0.6_imp0_2024-12-10.RData")


## calculating number of extremes per simulation
r0_data <- sims_imp0 %>% 
  unnest(sim_cond) %>% 
  mutate(n_tx35 = map(sim_df,
                      ~ sum(.$tx35)),
         n_tx40 = map(sim_df,
                      ~ sum(.$tx40)),
         n_tx35_30days = map(sim_df,
                             ~ sum(.$tx35[1:30])),
         n_tx40_30days = map(sim_df,
                             ~ sum(.$tx40[1:30])),
         season = if_else(period == "g1g2", "spring", "summer")) %>% 
  select(-A) %>% 
  unnest(starts_with("n_"))

rm(sims_imp0)
gc()

# median number of extremes per simulation set (sitexperiodxGWL)
med_tx40 <- r0_data %>% 
  group_by(sim_set) %>% 
  summarise(med_tx40_30days = median(n_tx40_30days),
            mean_tx40_30days = mean(n_tx40_30days))



# R0 ~ LLHI events --------------------
### linear models ------
lm_no_plant <- r0_data %>% 
  filter(site == "lowl", period == "g2g3", sim_with_plant == 0) %>% 
  lm(log(gr) ~ n_tx40_30days, data = ., na.action = na.omit) %>% 
  coefficients()

loggr_noheat <- r0_data %>% 
  filter(site == "lowl", period == "g2g3",
         n_tx40_30days == 0) %>% 
  group_by(sim_with_plant) %>% 
  summarise(mean_loggr = mean(log(gr), na.rm = T))


### extreme distribution ------
dist_llhi <- r0_data %>% 
  filter(site == "lowl", period == "g2g3") %>%
  mutate(n_tx40_30days = if_else(n_tx40_30days > 3, 3, n_tx40_30days)) %>% 
  group_by(GWL, n_tx40_30days) %>% 
  summarise(n = n()) %>% 
  group_by(GWL) %>% 
  mutate(prop = n/sum(n))


(llhi_bars <-  r0_data %>% 
  filter(site == "lowl", period == "g2g3") %>% 
  group_by(GWL, n_tx40_30days) %>% 
  summarise(n = n()) %>% 
  group_by(GWL) %>% 
  mutate(prop = n/sum(n),
         GWL = if_else(GWL == 1, .85, GWL)) %>% 
  ggplot(aes(x = n_tx40_30days, y = prop, fill = as.factor(GWL))) +
  geom_col(position = "dodge") +
  scale_x_continuous(breaks = 0:5,
                     name = "LLHI hot days month<sup>-1</sup>",
                     expand = expansion(add = c(.1, .1))) +
  scale_y_continuous(name = "Frequency",
                     breaks = c(0, .4, .8)) +
  scale_fill_grey(name = "GWL (°C)",
                  start = .8, end = .2) +
  theme_classic() +
  theme(axis.title.x = element_markdown(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.key.width = unit(5, "pt")))


r0_data %>% 
  filter(site == "lowl", period == "g2g3") %>%
  mutate(n_nonLLHI = n_tx35_30days - n_tx40_30days) %>% 
  group_by(GWL) %>% 
  summarise(across(.cols = starts_with("n_"),
                   .fns = median))

### lm plot -------
arrows <- tibble(n_tx40_30days = 0:3) %>% 
  mutate(pred_loggr = lm_no_plant[1]+lm_no_plant[2]*n_tx40_30days,
         pred_loggr = if_else(n_tx40_30days == 0,
                              loggr_noheat$mean_loggr[2],
                              pred_loggr),
         y = loggr_noheat$mean_loggr[1])

draw_key_vpath <- function(data, params, size) {
  segmentsGrob(0.5, 0.9, 0.5, 0.1,
               gp = gpar(
                 col = alpha(data$colour %||% data$fill %||% "black", data$alpha),
                 lwd = (data$linewidth %||% 0.5) * .pt,
                 lty = data$linetype %||% 1,
                 lineend = params$lineend %||% "butt"
               ),
               arrow = params$arrow
  )
}

(lm_r02 <- r0_data %>%
  filter(site == "lowl", period == "g2g3") %>%
  ggplot(aes(x = n_tx40_30days, y = log(gr))) +
  geom_smooth(aes(color = as.factor(sim_with_plant)),
              method = "lm", fullrange = T) +
  geom_point(data = arrows,
             aes(x = n_tx40_30days, y = y),
             size = 2) +
  geom_point(data = arrows,
             aes(x = n_tx40_30days, y = pred_loggr),
             size = 2) +
  geom_hline(aes(yintercept = 0), color = "red") +
  scale_color_manual(values = met.brewer("Kandinsky",
                                         n = 3)[1:2],
                     labels = c('0' = "Enough plant",
                                '1' = "Drought-induced<br>plant scarcity"),
                     guide = guide_legend(override.aes = list(fill = NA))) +
  new_scale_color() +
  geom_segment(data = arrows,
               aes(x = n_tx40_30days, xend = n_tx40_30days,
                   y = y, yend = pred_loggr,
                   color = "Effect"),
               arrow = arrow(angle = 15,
                             length = unit(7.5, "pt")),
               key_glyph = draw_key_vpath) +
  scale_color_manual(values = "black",
                     guide = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "LLHI hot days month<sup>-1</sup>",
       y = "ln R<sub>0,2</sub>") +
  scale_x_continuous(position = "top",
                     expand = expansion(add = c(.6,.6)),
                     breaks = 0:5) +
  theme_classic() +
  theme(axis.title.y = element_markdown(size = 10),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10),
        axis.title.x.top = element_markdown(size = 10),
        legend.text = element_markdown(size = 10),
        legend.title = element_blank()))
  




# Median effects -------------
## Data preparation ---------
med_r0_data <- r0_data %>% 
  filter(!is.na(gr)) %>% 
  mutate(extreme_group = paste0(sim_with_plant, as.numeric(n_tx35_30days != 0), as.numeric(n_tx40_30days != 0))) %>% 
  group_by(site, period, extreme_group) %>% 
  mutate(med_G000 = median(log(gr))) %>%
  # for the group without events (000) we can use the same median R0 in all GWL
  # this way we avoid sample sizes issues in the simulation sets with low samples with 000 conditions
  group_by(site, period, GWL, extreme_group) %>% 
  summarise(med_r0 = median(log(gr)),
            med_G000 = mean(med_G000)) %>% 
  mutate(med_r0 = if_else(extreme_group == '000', med_G000, med_r0),
         med_r0 = if_else(extreme_group == '100', med_G000, med_r0)) %>% 
  select(-med_G000) %>% 
  pivot_wider(id_cols = c(site, period, GWL),
              names_from = extreme_group,
              names_prefix = "G",
              values_from = med_r0) %>% 
  mutate(plant_eff = G100 - G000,
         tx35_eff = G010 - G000,
         'tx40+tx35_eff' = G011 - G000,
         comp_eff = G111 - G000) %>% 
  pivot_longer(cols = ends_with("eff"),
               names_to = "event",
               names_pattern = "(.*)_eff",
               values_to = "effect") %>% 
  mutate(event = factor(event,
                         levels = c("tx35", "tx40+tx35", "plant","comp")))

### additive effects ------
add_effects <- med_r0_data %>% 
  filter(event %in% c("plant", "tx40+tx35")) %>% 
  group_by(site, period, GWL) %>% 
  summarise(effect = sum(effect)) %>% 
  mutate(event = "comp",
         GWL = if_else(GWL == 1, .85, GWL))

med_r0_data %>% 
  filter(site == "lowl", period == "g2g3", event == "comp") %>% 
  mutate(GWL = if_else(GWL == 1, .85, GWL)) %>% 
  left_join(add_effects,
            by = join_by(site, period, GWL, event),
            suffix = c("int", "add")) %>% 
  mutate(prop = exp(effectadd-effectint),
         effectint_exp = exp(effectint),
         effectadd_exp = exp(effectadd),
         prop_red = effectadd_exp/effectint_exp)

## bar plot -----
(lowlr2_eff <- med_r0_data %>% 
  filter(site == "lowl", period == "g2g3") %>% 
   mutate(GWL = if_else(GWL == 1, .85, GWL)) %>% 
  ggplot(aes(x = event, y = effect)) +
  geom_col(aes(fill = event), position = position_dodge(),
           color = "black") +
  scale_fill_manual(values = met.brewer(name = "Hokusai1",
                                        n = 7)[c(c(3, 2, 5, 1))],
                    labels = c("tx35" = "non-LLHI heat event<br>(35 < T<sub>max</sub> < 40 °C)",
                               "tx40+tx35" = "Extreme heat event<br>(T<sub>max</sub> > 35 °C)",
                               "plant" = "Drought event",
                               "comp" = "Hot-dry compound<br>event")) +
  geom_col(data = add_effects,
           fill = NA, color = "black", linetype = "dashed") +
  facet_grid(. ~ GWL) +
  scale_x_discrete(name = "GWL (°C)", position = "top") +
   scale_y_continuous(breaks = seq(-1, 0, by = .5)) +
  labs(y = "Effect<br>(extreme ln R<sub>0,2</sub> - non-extreme ln R<sub>0,2</sub>)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10),
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.placement = "inside",
        legend.title = element_blank(),
        legend.text = element_markdown(size = 10),
        legend.key.height = unit(25, "pt"),
        legend.key.width = unit(10, "pt"),
        axis.title.y = element_markdown(size = 10),
        axis.text.y = element_text(size = 10),
        axis.line.x = element_blank()))


(effect_plot <- lowlr2_eff +
  (lm_r02 + llhi_bars +
  plot_layout(ncol = 1,
              heights = c(.35, .15))) +
    plot_layout(ncol = 1,
                heights = c(.2, .3)) +
    plot_annotation(tag_levels = "a",
                    tag_prefix = "(",
                    tag_suffix = ")") &
    theme(plot.tag = element_text(size = 10)))

ggsave("figures/paper_figures/Figure6.png",
       device = png,
       plot = effect_plot,
       width = 18,
       height = 21,
       units = "cm",
       dpi = 600)



# R0 ~ All types of ECEs (Fig. S12)--------

labs_reduced <- map_chr(0:15,
                        ~if_else(.%%3!=0, "", as.character(.)))

(all_ete_boxplot <- r0_data %>%
    mutate(n_tx_nonLLHI = n_tx35_30days-n_tx40_30days) %>%
    pivot_longer(cols = c(n_tx35_30days, n_tx40_30days, n_tx_nonLLHI),
                 names_to = "ETE_type",
                 values_to = "n_events") %>%
    filter(n_events <= 15) %>% 
    mutate(ETE_type = case_when(ETE_type == "n_tx35_30days" ~ "Extreme heat events<br>(T<sub>max</sub> > 35 °C)",
                                ETE_type == "n_tx40_30days" ~ "LLHI heat events<br>(T<sub>max</sub> > 40 °C)",
                                T ~ "non-LLHI heat events<br>(35 < T<sub>max</sub> < 40 °C)"),
           site_period = paste(site, period)) %>%
    ggplot(aes(x = factor(n_events,
                          levels = as.character(0:15)),
               y = log(gr))) +
    geom_hline(aes(yintercept = 0), color = "red") +
    geom_boxplot(aes(fill = as.factor(sim_with_plant)),
                 outlier.shape = NA,
                 linewidth = .25) +
    scale_fill_manual(values = met.brewer("Kandinsky", n = 3)[1:2],
                      labels = c('0' = "Enough plant",
                                 '1' = "Drought-induced plant scarcity")) +
    facet_grid(site_period ~ ETE_type,
               labeller = labeller(site_period = c("lowl g1g2" = "Declining<br>population<br>R<sub>0,1</sub>",
                                                   "lowl g2g3" = "Declining<br>population<br>R<sub>0,2</sub>",
                                                   "lowl g3g4" = "Declining<br>population<br>R<sub>0,3</sub>",
                                                   "mide g1g2" = "Non-declining<br>population<br>R<sub>0,1</sub>",
                                                   "mide g2g3" = "Non-declining<br>population<br>R<sub>0,2</sub>"))) +
    labs(x = "Events month<sup>-1</sup>",
         y = "ln R<sub>0</sub>") +
    scale_x_discrete(labels = labs_reduced) +
    scale_y_continuous(breaks = c(-4, -2, 0, 2)) +
    theme_classic() +
    theme(axis.title.y = element_markdown(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_markdown(size = 11),
          axis.text.x = element_text(size = 11),
          strip.background = element_blank(),
          strip.text = element_markdown(size = 11),
          strip.text.y = element_markdown(size = 11),
          legend.title = element_blank(),
          legend.text = element_markdown(size = 11),
          legend.position = "top") +
    coord_cartesian(ylim = c(-5, 2.5)))


ggsave(filename = "figures/paper_figures/FigureS12.png",
       plot = all_ete_boxplot,
       device = png,
       height = 21,
       width = 18,
       units = "cm",
       dpi = 600)