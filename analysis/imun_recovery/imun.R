##### libraries ----------------------------------------------------------------

library(tidyverse)
library(wesanderson)
library(gamm4)

set.seed(123)

source("utils/plot.R")

#### data preparation ----------------------------------------------------------

cd4_ch <- readRDS("data_clean/ch/cd4_ch.rds")%>% 
  mutate(id = as.factor(id),
         cohort = as.factor("CH")) %>% 
  select(-region)

cd4_rsa <- readRDS("data_clean/rsa/cd4_rsa.rds") %>% 
  mutate(id = as.factor(id),
         cohort = as.factor("RSA"))

cd4 <- rbind(cd4_ch %>% dplyr::select(-timepoint),
             cd4_rsa %>% dplyr::select(-timepoint)) %>% 
  filter(time_diff >= -60 & time_diff < 360) %>% 
  mutate(
    cd4_trans = sqrt(cd4))

cd4_npres <- cd4 %>% 
  filter(presenting_tb == 0)

cd4_pres <- cd4 %>% 
  filter(presenting_tb == 1)

#### complete record analysis [interaction] ---------------------------------------------------

cd4$interaction_term <- interaction(cd4$cohort, cd4$presenting_tb)

model_cd4_count <- gamm4(cd4_trans ~ cohort + presenting_tb  +
                           s(time_diff, by = interaction_term, bs = "cr", k = 4),
                                            random = ~ (1 | id),
                                            data = cd4)
summary(model_cd4_count$gam)

pred <- pred_trend(model_cd4_count, data = cd4)

trend_cd4 <- ggplot(pred, aes(x = time_diff)) +
  facet_wrap(~presenting_tb) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = cohort), alpha = 0.2) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 
  
trend_cd4

ggsave(trend_cd4, filename = "results/cd4_predEff.png",
       width = 20, height = 12.7, unit = "cm")

#### complete record analysis [separate models] -------------------------------------

### non-presenting ###

model_cd4_count_npres <- gamm4(cd4_trans ~ cohort +
                                s(time_diff, by = cohort, bs = "cr", k = 4),
                              random = ~ (1 | id),
                              data = cd4_npres)

pred_npres <- pred_trend_uni(model_cd4_count_npres, 
                   data = cd4_npres)

ggplot(pred_npres, aes(x = time_diff)) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 600)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP, fill = cohort), alpha = 0.2) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 

### presenting ###

model_cd4_count_pres <- gamm4(cd4_trans ~ cohort +
                           s(time_diff, by = cohort, bs = "cr", k = 4),
                         random = ~ (1 | id),
                         data = cd4_pres)

pred_pres <- pred_trend_uni(model_cd4_count_pres, 
                         data = cd4_pres)

ggplot(pred_pres, aes(x = time_diff)) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = cohort), alpha = 0.2) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 

