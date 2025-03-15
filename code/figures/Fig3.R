# Trends in observed growth rates


# Packages ----------------------------------------------------------------

library(tidyverse)
library(broom)
library(MetBrewer)
library(ggtext)
library(patchwork)
theme_set(theme_classic(base_family = "Arial"))



# Data --------------------------------------------------------------------

abund <- read_csv("data/cbms/cbms_nap_1_9_adultsurv.csv")


# Estimating growth rates between consecutive generations
## Taking abundance as sum of counts during the whole generation
brood_abund <- abund %>%
  mutate(brood = case_when(mes < 5 ~ 1,
                           mes < 7 ~ 2,
                           mes < 8 & itin == 1 ~ 3,
                           itin == 9 ~ 3,
                           TRUE ~ 4)) %>% 
  group_by(itin, esp, anys, brood) %>% 
  summarise(brood_abund = sum(ABUND)) %>% 
  group_by(itin, esp, anys) %>% 
  mutate(next_abund = lead(brood_abund),
         growth_rate = next_abund/brood_abund,
         rate = if_else(growth_rate > 1, "Increasing", "Decreasing"),
         broods = case_when(brood == 1 ~ "R<sub>0,1</sub>",
                            brood == 2 ~ "R<sub>0,2</sub>",
                            brood == 3 ~ "R<sub>0,3</sub>",
                            T ~ NA_character_),
         site = if_else(itin == 1, "Declining population", "Non-declining population")) %>% 
  filter(brood < 4,
         !(brood == 3 & itin == 9))


# Estimating growth rates between consecutive years
## Taking abundance as sum of counts during the whole year
ann_abund <- abund %>%
  group_by(itin, esp, anys) %>% 
  summarise(ann_abund = sum(ABUND)) %>% 
  mutate(next_abund = lead(ann_abund),
         growth_rate = next_abund/ann_abund,
         rate = if_else(growth_rate > 1, "Increasing", "Decreasing"),
         site = if_else(itin == 1, "Declining population", "Non-declining population")) %>% 
  filter(esp == "piepna",
         anys > 1988,
         !is.na(growth_rate))


# Plot --------------------------------------------------------------------
## Interbrood net reproductive rates ---------------------------------------
(brood_plots_log <- brood_abund %>% 
    left_join(brood_ab_mod) %>% 
    split(.$site) %>%
    map2(mean_brood_gr,
         ~ ggplot(data = .x,
                  aes(x = anys, y = log(growth_rate))) +
           geom_segment(aes(yend = 0,
                            xend = anys,
                            alpha = rate),
                        color = met.brewer("Kandinsky", n = 3)[3],
                        size = 1.5) +
           geom_hline(data = .y,
                      aes(yintercept = log(geommean_gr),
                          # alpha = rate,
                          color = "Geometric mean"),
                      linetype = "dotted",
                      size = .5) +
           scale_alpha_manual(values = c(1, .4),
                              name = "Rate",
                              guide = guide_legend(title.hjust = .5,
                                                   keywidth = .2,
                                                   keyheight = 1,
                                                   reverse = T)) +
           geom_hline(aes(yintercept = 0),
                      color = "black",
                      size = .25) +
           geom_smooth(aes(#linetype = sign,
             color = "Multiannual trend"),
             linetype = "dashed",
             method = "lm",
             size = .5,
             se = F) +
           scale_color_manual(values = met.brewer("Kandinsky", n = 3)[c(3, 3)],
                              name = "",
                              guide = guide_legend(override.aes = list(linetype = c("dotted", "dashed")))) +
           facet_grid(cols = vars(broods),
                      scales = "free") +
           scale_x_continuous(breaks = seq(from = 1995, to = 2015, by = 10)) +
           ylim(-3.8, 3.5)) %>% 
    map2(names(.),
         ~ . +
           labs(y = "ln R<sub>0</sub>",
                x = "Years",
                subtitle = "") +
           theme_classic() +
           theme(strip.text = element_markdown(size = 10),
                 strip.background = element_blank(),
                 legend.title = element_text(size = 10),
                 legend.text = element_text(size = 10),
                 legend.box.just = "center",
                 axis.text = element_text(size = 10),
                 plot.subtitle = element_text(size = 10),
                 axis.title.y = element_markdown(size = 10),
                 axis.title.x = element_text(size = 10))))

## Annual growth rate ------------------------------------------------------
(ann_plots_log <- ann_abund %>% 
    left_join(ann_ab_mod) %>% 
    split(.$site) %>%
    map2(mean_gr,
         ~ ggplot(data = .x,
                  aes(x = anys, y = log(growth_rate))) +
           facet_wrap(facets = vars(empty_facet)) +
           geom_segment(aes(yend = 0,
                            xend = anys,
                            alpha = rate),
                        color = met.brewer("Kandinsky", n = 3)[3],
                        size = 1.5) +
           geom_hline(data = .y,
                      aes(yintercept = log(geommean_gr),
                          color = "Geometric mean"),
                      linetype = "dotted",
                      size = .5) +
           scale_alpha_manual(values = c(1, .4),
                              name = "Rate",
                              guide = guide_legend(title.hjust = .5,
                                                   keywidth = .2,
                                                   keyheight = 1,
                                                   reverse = T)) +
           geom_hline(aes(yintercept = 0),
                      color = "black",
                      size = .25) +
           geom_smooth(aes(color = "Multiannual trend"),
                       method = "lm",
                       linetype = "dashed",
                       size = .5,
                       se = F) +
           scale_linetype_manual(values = c(0, 1),
                                 guide = "none") +
           scale_color_manual(values = met.brewer("Kandinsky", n = 3)[c(3, 3)],
                              name = "",
                              guide = guide_legend(override.aes = list(linetype = c("dotted", "dashed")))) +
           scale_x_continuous(breaks = seq(from = 1995, to = 2015, by = 10)) +
           ylim(-3.8, 3.5) +
           guides()) %>% 
    map2(c("Declining population", "Non-declining population"),
         ~ .x  +
           labs(y = "ln &lambda;<sub>annual</sub>",
                x = "Years",
                subtitle = .y) +
           theme_classic() +
           theme(strip.background = element_blank(),
                 strip.text = element_markdown(size = 10),
                 legend.title = element_text(size = 10),
                 legend.text = element_markdown(size = 10),
                 axis.text = element_text(size = 10),
                 plot.subtitle = element_text(size = 10),
                 axis.title.x = element_text(size = 10),
                 axis.title.y = element_markdown(size = 10))))


## Together ----------------------------------------------------------------
(ann_log <- (ann_plots_log$`Declining population`+
           theme(legend.position = "none")) /
    (ann_plots_log$`Non-declining population` +
       theme(legend.position = "none")) +
    plot_layout(ncol = 1))


(broods_log <- (brood_plots_log$`Declining population` +
              theme(legend.position = "none")) +
    (brood_plots_log$`Non-declining population` + guide_area() +
       plot_layout(ncol = 2,
                   guides = "collect",
                   widths = c(2, 1))) +
    plot_layout(ncol = 1))


(both_log <- (ann_log | broods_log) +
    plot_layout(widths = c(1, 3)) &
    theme(legend.box.just = "center"))

ggsave(filename = "figures/paper_figures/Figure3.png",
       device = png,
       plot = both_log,
       width = 21,
       height = 15,
       units = "cm",
       dpi = 600)
