#### libraries -----------------------------------------------------------------

library(tidyverse)
library(ggsurvfit)
library(tidycmprsk)
library(cmprsk)
library(ggpubr)

##### data preparation ---------------------------------------------------------
rna_ch <- readRDS("data_clean/ch/rna_ch.rds") %>% 
  dplyr::select(id, rna, presenting_tb, time_diff)
rna_rsa <- readRDS("data_clean/rsa/rna_rsa.rds") %>% mutate(region = NA) %>% 
  dplyr::select(id, rna, presenting_tb, time_diff)
rna <- rbind(rna_ch, rna_rsa)

df <- readRDS("data_clean/art.rds") %>% 
  dplyr::select(id, presenting_tb, cohort, last_persontime, exitdate, gender, age_at_art_start, region, who_stage, cd4_group) %>% 
  mutate(region = relevel(region, ref = "Europe/Northern America"),
         cohort = relevel(cohort, ref = "RSA"))

valid_patients <- rna %>% # patients without any labdata shouldnt be treated as "followed"
  filter(!is.na(rna)) %>%
  distinct(id)

below_400 <- rna %>% 
  dplyr::select(id, rna, time_diff) %>% 
  filter(time_diff >= 0 & rna <= 400) %>%
  arrange(id, time_diff) %>% 
  group_by(id) %>%
  filter(time_diff == min(time_diff)) %>%  # Select the rows with the smallest time_diff
  ungroup()

custom_breaks <- c(16, 34, 44, 100)

aj_400 <- valid_patients %>% 
  left_join(below_400, by = "id") %>% 
  left_join(df, by = "id") %>% 
  mutate(time = case_when(!is.na(time_diff) ~ time_diff,
                          TRUE ~ last_persontime),
         event = as.factor(case_when(!is.na(time_diff) ~ 1,
                                     !is.na(exitdate) ~ 2,
                                     TRUE ~ 0)),
         suppression = as.factor(ifelse(!is.na(time_diff),1,0)),
         agegroup = cut(age_at_art_start, breaks = custom_breaks, labels = c("16-34","35-44", "45+"), include.lowest = TRUE)) %>% 
  filter(time >= 0) %>% 
  mutate(stratum = as.factor(paste0(cohort, "_", presenting_tb))) %>% 
  rename(age = age_at_art_start)  %>% 
  mutate(gender = fct_drop(gender),
         cd4_group = fct_drop(cd4_group),
         cd4_group = fct_relevel(cd4_group, "350+")) 

aj_400_ch <- aj_350 %>% filter(cohort == "CH")
aj_400_rsa <- aj_350 %>% filter(cohort == "RSA")

##### aalen-johansen model -----------------------------------------------------

options("ggsurvfit.switch-color-linetype" = TRUE)

aj_rsa400 <- aj_400 %>% filter(cohort == "RSA")
aj_ch400 <- aj_400 %>% filter(cohort == "CH")

plot_aj_rsa400 <- survfit2(Surv(time, event, type = "mstate") ~ presenting_tb, data = aj_rsa400) |>
  ggcuminc(outcome = c("1"),
           aes(linetype = presenting_tb),
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[1]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[1]) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = c(0,362), ylim = c(0,1))+
  scale_x_continuous(expand = c(0,0),  breaks = seq(0, 360, 180)) +
  scale_y_continuous(expand = c(0,0), 
                     labels = function(y) paste0(y * 100)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "none",
        legend.title = element_blank(), 
        plot.title = element_text(size = 20, face = "bold"), # Set title properties here
        plot.title.position = "plot",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) 
  labs(title = expression("1b | Time to viral suppression VL < 400 (copies/ml)"))

plot_aj_rsa400

plot_aj_ch400 <- survfit2(Surv(time, event, type = "mstate") ~ presenting_tb, data = aj_ch400) |>
  ggcuminc(outcome = c("1"),
           aes(linetype = presenting_tb),
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[2]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[2]) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = c(0,362), ylim = c(0,1))+
  scale_x_continuous(expand = c(0,0),  breaks = seq(0, 360, 180)) +
  scale_y_continuous(expand = c(0,0), 
                     labels = function(y) paste0(y * 100)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "none",
        legend.title = element_blank(), 
        plot.title = element_text(size = 20, face = "bold"), # Set title properties here
        plot.title.position = "plot",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) 

plot_aj_ch400

aj_sup_both400 <- ggarrange(plot_aj_rsa400, plot_aj_ch400, ncol = 2) %>% 
  annotate_figure(
                  left = text_grob("Proportion (%)",size=20, rot = 90),
                  top = text_grob("2a | Time to viral suppression VL < 400 (copies/ml)",size=20, hjust = 0.77))
aj_sup_both400

ggsave(plot = aj_sup_both400, filename = "results/suppression/suppression_400.png", units = "cm")
