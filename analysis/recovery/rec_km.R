#### libraries -----------------------------------------------------------------

library(tidyverse)
library(ggsurvfit)
library(tidycmprsk)
library(cmprsk)
library(ggpubr)

##### data preparation ---------------------------------------------------------

cd4_ch <- readRDS("data_clean/ch/cd4_ch.rds") 
cd4_rsa <- readRDS("data_clean/rsa/cd4_rsa.rds") %>% mutate(region = NA)
cd4 <- rbind(cd4_ch, cd4_rsa)
df <- readRDS("data_clean/art.rds") %>% 
  dplyr::select(id, presenting_tb, cohort, last_persontime, exitdate, gender, age_at_art_start, region, who_stage, cd4_group) %>% 
  mutate(region = relevel(region, ref = "Europe/Northern America"))

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

custom_breaks <- c(16, 34, 44, 100)

aj_350 <- valid_patients %>% 
  left_join(above_350, by = "id") %>% 
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
         cd4_group = fct_relevel(cd4_group, "350+")) %>% 
  dplyr::select(id, presenting_tb, time, cd4, event, cohort, stratum, gender, age, suppression, region, who_stage, cd4_group, agegroup)
  
aj_350_ch <- aj_350 %>% filter(cohort == "CH")
aj_350_rsa <- aj_350 %>% filter(cohort == "RSA")

##### aalen-johansen model -----------------------------------------------------

options("ggsurvfit.switch-color-linetype" = TRUE)

### 350 ###

aj_rsa350 <- aj_350 %>% filter(cohort == "RSA")
aj_ch350 <- aj_350 %>% filter(cohort == "CH")

plot_aj_rsa350 <- survfit2(Surv(time, event, type = "mstate") ~ presenting_tb, data = aj_rsa350) |>
  ggcuminc(outcome = c("1"),
           aes(linetype = presenting_tb),
           linewidth = 1.5, 
           color = wes_palette("Moonrise2")[1]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[1]) +
  theme_classic(base_size = 20) +
  coord_cartesian(xlim = c(0,362), ylim = c(0,1))+
  scale_x_continuous(expand = c(0,0),  breaks = seq(0, 360, 60)) + 
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
        plot.background = element_rect(fill='transparent', color=NA),
        plot.margin = margin(t = 5, r = 18, b = 5, l = 5, unit = "pt")) 

plot_aj_rsa350

plot_aj_ch350 <- survfit2(Surv(time, event, type = "mstate") ~ presenting_tb, data = aj_ch350) |>
  ggcuminc(outcome = c("1"),
           aes(linetype = presenting_tb),
           linewidth = 1.5, 
           color = wes_palette("Moonrise2")[2]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[2]) +
  theme_classic(base_size = 18) +
  coord_cartesian(xlim = c(0,362), ylim = c(0,1))+
  scale_x_continuous(expand = c(0,0),  breaks = seq(0, 360, 60)) +  
  scale_y_continuous(expand = c(0,0), 
                     labels = function(y) paste0(y * 100)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "none",
        legend.title = element_blank(), 
        plot.title = element_text(size = 18, face = "bold"), # Set title properties here
        plot.title.position = "plot",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        plot.margin = margin(t = 5, r = 18, b = 5, l = 5, unit = "pt")) 


plot_aj_ch350

aj_rec_both350 <- ggarrange(plot_aj_rsa350,plot_aj_ch350, ncol = 2) %>% 
  annotate_figure(left = text_grob("Proportion (%)", size = 18, rot = 90),
                  top = text_grob(expression(bold("B | ") * "Time to immun recovery: CD4 cell count > 350 (cells/µl)"), size=18, hjust = 0.73, vjust = .1))

aj_rec_both350

ggsave(plot = aj_rec_both350, filename = "results/suppression/recovery_350.png", 
       width = 16, height = 11, units = "cm")

#### logistic model -----------------------------------------------------------

### CH ###

m1_ch <- glm(suppression ~ gender + agegroup + who_stage + region, 
             family = "binomial", 
             data = aj_350_ch)

summ(m1_ch, exp = T)

### RSA ###

m1_rsa <- glm(suppression ~ gender + agegroup , 
              family = "binomial", 
              data = aj_350_ch)

summ(m1_rsa, exp = T)

#### fine gray model ----

library(tidyverse)
library(cmprsk)
library(survival)
library(finalfit)

explanatory   <- c("cohort","agegroup", "gender")
dependent_crr <- "Surv(time, event)"

aj_350 <-  aj_350 %>% filter_all(all_vars(!is.na(.)))  

aj_350 %>%
  # Summary table
  summary_factorlist(dependent_crr, explanatory, 
                     column = TRUE, fit_id = TRUE) %>% 
  # Fine and Gray competing risks regression
  ff_merge(
    aj_350 %>%
      crrmulti(dependent_crr, explanatory) %>%
      fit2df(estimate_suffix = " (competing risks multivariable)")) %>% 
  dplyr::select(-fit_id, -index) %>%
  dependent_label(aj_350, "first CD4 count > 350") 
