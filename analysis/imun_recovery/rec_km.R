#### libraries -----------------------------------------------------------------

library(tidyverse)
library(ggsurvfit)
library(tidycmprsk)
library(cmprsk)
library(ggpubr)

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  ggplot2, # for plotting
  survival,
  survminer,
  gridExtra,
  Epi
)

##### data preparation ---------------------------------------------------------

cd4_ch <- readRDS("data_clean/ch/cd4_ch.rds") %>% select(-region)
cd4_rsa <- readRDS("data_clean/rsa/cd4_rsa.rds")
cd4 <- rbind(cd4_ch, cd4_rsa)
df <- readRDS("data_clean/art.rds") %>% 
  select(id, presenting_tb, cohort, last_persontime, exitdate)

valid_patients <- cd4 %>% # patients without any labdata shouldnt be treated as "followed"
  filter(!is.na(cd4)) %>%
  distinct(id)

above_350 <- cd4 %>% 
  dplyr::select(id, cd4, time_diff) %>% 
  filter(time_diff >= 0 & cd4 >= 350) %>%
  arrange(id, time_diff) %>% 
  group_by(id) %>%
  filter(time_diff == min(time_diff)) %>%  # Select the rows with the smallest time_diff
  ungroup()

above_500 <- cd4 %>% 
  dplyr::select(id, cd4, time_diff) %>% 
  filter(time_diff >= 0 & cd4 >= 500) %>%
  arrange(id, time_diff) %>% 
  group_by(id) %>%
  filter(time_diff == min(time_diff)) %>%  # Select the rows with the smallest time_diff
  ungroup()

aj_350 <- valid_patients %>% 
  left_join(above_350, by = "id") %>% 
  left_join(df, by = "id") %>% 
  mutate(time = case_when(!is.na(time_diff) ~ time_diff,
                          TRUE ~ last_persontime),
         event = as.factor(case_when(!is.na(time_diff) ~ 1,
                                     !is.na(exitdate) ~ 2,
                                     TRUE ~ 0))) %>% 
  filter(time >= 0) %>% 
  mutate(stratum = as.factor(paste0(cohort, "_", presenting_tb))) %>% 
  dplyr::select(id, presenting_tb, time, cd4, event, cohort, stratum)

aj_500 <- valid_patients %>% 
  left_join(above_500, by = "id") %>% 
  left_join(df, by = "id") %>% 
  mutate(time = case_when(!is.na(time_diff) ~ time_diff,
                          TRUE ~ last_persontime),
         event = as.factor(case_when(!is.na(time_diff) ~ 1,
                                     !is.na(exitdate) ~ 2,
                                     TRUE ~ 0))) %>% 
  filter(time >= 0) %>% 
  mutate(stratum = as.factor(paste0(cohort, "_", presenting_tb))) %>% 
  dplyr::select(id, presenting_tb, time, cd4, event, cohort, stratum)

##### aalen-johansen model -----------------------------------------------------

options("ggsurvfit.switch-color-linetype" = TRUE)

### 350 ###

aj_rsa350 <- aj_350 %>% filter(cohort == "RSA")
aj_ch350 <- aj_350 %>% filter(cohort == "CH")

plot_aj_rsa350 <- survfit2(Surv(time, event, type = "mstate") ~ presenting_tb, data = aj_rsa350) |>
  ggcuminc(outcome = c("1"),
           aes(linetype = presenting_tb),
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[1]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[1]) +
  theme_classic() +
  coord_cartesian(xlim = c(0,2000), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0,0), 
                     breaks = function(limits) pretty(limits, n = 5, integer = TRUE)) +  
  scale_y_continuous(expand = c(0,0), 
                     labels = function(y) paste0(y * 100)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

plot_aj_rsa350

plot_aj_ch350 <- survfit2(Surv(time, event, type = "mstate") ~ presenting_tb, data = aj_ch350) |>
  ggcuminc(outcome = c("1"),
           aes(linetype = presenting_tb),
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[2]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[2]) +
  theme_classic() +
  coord_cartesian(xlim = c(0,2000), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0,0), 
                     breaks = function(limits) pretty(limits, n = 5, integer = TRUE)) +  
  scale_y_continuous(expand = c(0,0), 
                     labels = function(y) paste0(y * 100)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

plot_aj_ch350

aj_rec_both350 <- ggarrange(plot_aj_rsa350, plot_aj_ch350, ncol = 2) %>% 
  annotate_figure(bottom = "Days after ART start",
                  left = "Cumulative occurence of cd4 recovery (%)")
aj_rec_both350

ggsave(plot = aj_rec_both350, filename = "results/recovery_350.png", 
       width = 16, height = 11, units = "cm")

### 500 ###

aj_rsa500 <- aj_500 %>% filter(cohort == "RSA")
aj_ch500 <- aj_500 %>% filter(cohort == "CH")

plot_aj_rsa500 <- survfit2(Surv(time, event, type = "mstate") ~ presenting_tb, data = aj_rsa500) |>
  ggcuminc(outcome = c("1"),
           aes(linetype = presenting_tb),
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[1]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[1]) +
  theme_classic() +
  coord_cartesian(xlim = c(0,2000), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0,0), 
                     breaks = function(limits) pretty(limits, n = 5, integer = TRUE)) +  
  scale_y_continuous(expand = c(0,0), 
                     labels = function(y) paste0(y * 100)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

plot_aj_rsa500

plot_aj_ch500 <- survfit2(Surv(time, event, type = "mstate") ~ presenting_tb, data = aj_ch500) |>
  ggcuminc(outcome = c("1"),
           aes(linetype = presenting_tb),
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[2]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[2]) +
  theme_classic() +
  coord_cartesian(xlim = c(0,2000), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0,0), 
                     breaks = function(limits) pretty(limits, n = 5, integer = TRUE)) +  
  scale_y_continuous(expand = c(0,0), 
                     labels = function(y) paste0(y * 100)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

plot_aj_ch500

aj_rec_both500 <- ggarrange(plot_aj_rsa500, plot_aj_ch500, ncol = 2) %>% 
  annotate_figure(bottom = "Days after ART start",
                  left = "Cumulative occurence of cd4 recovery (%)")
aj_rec_both500

ggsave(plot = aj_rec_both500, filename = "results/recovery_500.png", 
       width = 16, height = 11, units = "cm")
