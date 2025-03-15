# Figure of projected R0


# Packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(ggtext)
library(gghalves)
library(ggdist)
library(distributional)
library(zoo)
library(MetBrewer)
theme_set(theme_classic(base_family = "Arial"))



# Data --------------------------------------------------------------------

## Projected R0's
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


# Plot --------------------------------------------------------------------

(lowl <- sims_imp0 %>% 
  group_by(sim_set) %>% 
  mutate(q1 = quantile(log(gr), probs = .25, na.rm = T),
         q3 = quantile(log(gr), probs = .75, na.rm = T),
         iqr = q3-q1,
         max_whisk = q3 + 1.5*iqr,
         min_whisk = q1 - 1.5*iqr,
         gr = if_else(log(gr) > max_whisk | log(gr) < min_whisk, NA_real_, gr)) %>% #filtering outliers
  filter(!is.na(gr)) %>%
  select(-A) %>% 
  unnest(sim_cond) %>% 
  rename(net_repr_rate = gr) %>% 
  bind_rows(brood_abund, .id = "origin") %>%
  mutate(origin = case_when(origin == 1 & sim_with_plant == 0 ~ "Simulated with\nenough plant",
                            origin == 1 ~ "Simulated with\nplant scarcity",
                            T ~ "Observed")) %>% 
  filter(site == "lowl") %>% 
  ggplot(aes(x = as.factor(GWL), y = log(net_repr_rate))) +
  stat_slab(aes(fill = origin),
            slab_color = "black",
            slab_size = .5,
            orientation = "vertical",
            normalize = "groups",
            slab_alpha = .75,
            scale = .5) +
  scale_fill_manual(values = met.brewer("Kandinsky", n = 3)[c(3, 1, 2)],
                    name = NULL,
                    guide = guide_legend(
                      order = 1,
                      label.hjust = .5,
                      label.vjust = .5)) +
  geom_half_boxplot(center = T, outlier.shape = NA, errorbar.draw = F, nudge = .05, width = .5) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_vline(aes(xintercept = 2.6), linetype = "dashed") +
  geom_text(data = data.frame(x = c(1.5, 3.5, 5),
                              y = c(2, 2, 2),
                              txt = c("Current\nscenario", "Paris\nAgreement", "High-\nemissions"),
                              hjust = c(.5, .5, .25)),
            aes(x = x, y = y, label = txt, hjust = hjust),
            size = 3) +
   geom_text(data = data.frame(x = c(2.45, 2.75),
                               y = c(-4.5, -4.5),
                               txt = c("Validation", "Forecast"),
                               hjust = c(1, 0)),
             aes(x = x, y = y, label = txt, hjust = hjust),
             size = 3) +
  facet_grid(. ~ period, labeller = labeller(period = c(g1g2 = "R<sub>0,1</sub>",
                                                        g2g3 = "R<sub>0,2</sub>",
                                                        g3g4 = "R<sub>0,3</sub>"))) +
  labs(x = "Global warming level (°C)",
       y = "ln R<sub>0</sub>",
       title = "Declining population") +
  scale_x_discrete(labels = c("+0.85", "+0.85", "+1.5", "+2", "+4")) +
  scale_y_continuous(breaks = c(-4, -2, 0, 2)) +
  theme_classic() +
  theme(text = element_text(size = 10),
        axis.title.x = element_markdown(size = 10),
        axis.title.y = element_markdown(size = 10),
        strip.text = element_markdown(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10),
        strip.background = element_blank()))




(mide <- sims_imp0 %>%
    group_by(sim_set) %>% 
    mutate(q1 = quantile(log(gr), probs = .25, na.rm = T),
           q3 = quantile(log(gr), probs = .75, na.rm = T),
           iqr = q3-q1,
           max_whisk = q3 + 1.5*iqr,
           min_whisk = q1 - 1.5*iqr,
           gr = if_else(log(gr) > max_whisk | log(gr) < min_whisk, NA_real_, gr)) %>% #filtering outliers
    filter(!is.na(gr)) %>%
    select(-A) %>% 
    unnest(sim_cond) %>% 
    rename(net_repr_rate = gr) %>% 
    bind_rows(brood_abund, .id = "origin") %>%
    mutate(origin = case_when(origin == 1 & sim_with_plant == 0 ~ "Projected with\nenough plant",
                              origin == 1 ~ "Projected with\nplant scarcity",
                              T ~ "Observed"),
           origin = factor(origin,
                           levels = c("Observed",
                                      "Projected with\nenough plant",
                                      "Projected with\nplant scarcity"))) %>% 
    filter(site == "mide") %>% 
    ggplot(aes(x = as.factor(GWL), y = log(net_repr_rate))) +
    stat_slab(aes(fill = origin),
              slab_color = "black",
              slab_size = .5,
              orientation = "vertical",
              normalize = "groups",
              slab_alpha = .75,
              scale = .5,
              show.legend = T) +
    scale_fill_manual(
      values = met.brewer("Kandinsky", n = 3)[c(3, 1, 2)],
      name = NULL,
      drop = F,
      guide = guide_legend(
        order = 1,
        label.hjust = .5,
        label.vjust = .5)) +
    geom_half_boxplot(center = T, outlier.shape = NA, errorbar.draw = F, nudge = .05, width = .5) +
    geom_hline(aes(yintercept = 0), color = "red") +
    geom_vline(aes(xintercept = 2.6), linetype = "dashed") +
    geom_text(data = data.frame(x = c(1.5, 3.5, 5),
                                y = c(4.5, 4.5, 4.5),
                                txt = c("Current\nscenario", "Paris\nAgreement", "High-\nemissions"),
                                hjust = c(.5, .5, .25)),
              aes(x = x, y = y, label = txt, hjust = hjust),
              size = 3) +
    geom_text(data = data.frame(x = c(2.45, 2.75),
                                y = c(-4.5, -4.5),
                                txt = c("Validation", "Forecast"),
                                hjust = c(1, 0)),
              aes(x = x, y = y, label = txt, hjust = hjust),
              size = 3) +
    facet_grid(. ~ period, labeller = labeller(period = c(g1g2 = "R<sub>0,1</sub>",
                                                          g2g3 = "R<sub>0,2</sub>"))) +
    labs(x = "Global warming level (°C)",
         y = "ln R<sub>0</sub>",
         title = "Non-declining population") +
    scale_x_discrete(labels = c("+0.85", "+0.85", "+1.5", "+2", "+4")) +
    scale_y_continuous(breaks = c(-4, -2, 0, 2, 4)) +
    theme_classic() +
    theme(text = element_text(size = 10),
          axis.title.x = element_markdown(size = 10),
          plot.title = element_text(size = 10),
          axis.title.y = element_markdown(size = 10),
          strip.text = element_markdown(size = 10),
          legend.text = element_text(size = 10, margin = margin(t = 4, b = 4)),
          legend.key.height = unit(.5, units = "cm"),
          strip.background = element_blank()))



(both <- (lowl +
            labs(x = NULL) +
            theme(legend.position = "none")) /
    (mide +
       guide_area() +
       plot_layout(nrow = 1,
                   guides = "collect",
                   widths = c(0.7, 0.33))) &
    theme(legend.box = "vertical",
          legend.box.just = "center",
          legend.margin = margin(b = .5, t = .5, unit = "cm")))


ggsave("figures/paper_figures/Figure5.png",
       device = png,
       plot = both,
       width = 21,
       height = 21,
       units = "cm",
       dpi = 600)



# Calculation of the reductions in R0 -----

reductions <- sims_imp0 %>% 
  group_by(sim_set) %>% 
  mutate(q1 = quantile(log(gr), probs = .25, na.rm = T),
         q3 = quantile(log(gr), probs = .75, na.rm = T),
         iqr = q3-q1,
         max_whisk = q3 + 1.5*iqr,
         min_whisk = q1 - 1.5*iqr,
         gr = if_else(log(gr) > max_whisk | log(gr) < min_whisk, NA_real_, gr)) %>% #filtering outliers
  filter(!is.na(gr)) %>%
  select(-A) %>% 
  unnest(sim_cond) %>% 
  group_by(sim_set, site, period, GWL) %>% 
  summarise(med_gr = median(gr))

reductions %>% 
  pivot_wider(id_cols = c(site, period),
              names_from = GWL,
              values_from = med_gr,
              names_prefix = "gr") %>% 
  rename(gr15 = gr1.5) %>% 
  mutate(across(.cols = gr15:gr4,
                .fns = ~ (1-./gr1)*100,
                .names = "division_{.col}"))

1-0.659

