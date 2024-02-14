##### libraries ----------------------------------------------------------------

library(tidyverse)
library(wesanderson)
library(ggsurvfit)
library(tidycmprsk)
library(cmprsk)

##### data preparation ---------------------------------------------------------

df_surv <- readRDS("data_clean/art_noTB.rds")  %>% 
  mutate(incident_tb = as.numeric(as.character(incident_tb)),
         persontime_days = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")),
           incident_tb == 0 ~ last_persontime),
         persontime_years = persontime_days / 360,
         persontime_death_days = case_when(
           !is.na(exitdate) ~ as.numeric(difftime(exitdate, art_start_date, units = "days")),
           is.na(exitdate) ~ last_persontime),
         cohort = fct_relevel(cohort, "RSA")) %>% 
  dplyr::select(id, cohort, art_start_date, incident_tb, date_tb, last_persontime, persontime_years, persontime_days, persontime_death_days, cd4_group, exitdate, last_fup_date, rna_group) %>% 
  filter(persontime_years > 0) %>% 
  mutate(event_type = as.factor(case_when(
    incident_tb == 1 ~ 1,
    !is.na(exitdate) ~ 2,
    TRUE ~ 0 ))) # Loss to fup is not considered a competing risk

#### aalen-Johansen model [all-cause mortality as competing risk] --------------

#' I use ggcuminc to plot the competing risk survival function, it is mentioned in the documention
#' that both survfit2 (mstate) aswell as cuminc can be used to create the dataframe for plotting.

plot_aj <- survfit2(Surv(persontime_years, 
                         event_type, 
                         type = "mstate") ~ cohort, 
                    data = df_surv) %>% 
  ggcuminc(linewidth = 0.5, outcome = c("1")) +
  add_confidence_interval() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Years after ART start",
       y = "Cumulative incidence of TB (%)") +
  coord_cartesian(xlim = c(0,5), ylim = c(0,0.05))+
  scale_x_continuous(expand = c(0,0), 
                     breaks = function(limits) pretty(limits, n = 5, integer = TRUE)) +  
  scale_y_continuous(expand = c(0,0), 
                     labels = function(y) paste0(y * 100)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2"))

plot_aj

ggsave(plot = plot_aj, filename = "results/incidenceTB/aj.png", 
       width = 16, height = 11, units = "cm")

#### gray's test ---------------------------------------------------------------

cif_by_group <- cmprsk::cuminc(ftime = df_surv$persontime_years, 
                       fstatus = df_surv$event_type, 
                       group = df_surv$cohort, 
                       cencode = 0)

cif_by_group

# plot_aj_log <- survfit2(Surv(persontime_years, event_type, type = "mstate") ~ cohort, data = kaplan) |>
#   ggcuminc(linewidth = 0.5) +
#   add_confidence_interval() +
#   theme_classic() +
#   theme(legend.position = "none") +
#   labs(x = "Years after ART start",
#        y = "Cumulative incidence of TB (%)") +
#   scale_y_log10(expand = c(0,0), labels = function(y) paste0(y * 100))+
#   scale_x_continuous(expand = c(0,0), breaks = function(limits) pretty(limits, n = 10, integer = TRUE)) +  
#   theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
#   scale_color_manual(values = wes_palette("Moonrise2")) +
#   scale_fill_manual(values = wes_palette("Moonrise2")) 
# 
# plot_aj_log
# 
# ggsave(plot = plot_aj_log, filename = "results/incidence/incidence_aj_log.png", 
#        width = 16, height = 11, units = "cm")

## Seperate plots ##

# plot_aj_rsa <- survfit2(Surv(persontime_years, event_tb) ~ 1, data = kaplan_rsa) |>
# ggcuminc(linewidth = 0.5, color = wes_palette("Moonrise2")[1]) +
#   add_confidence_interval(fill = wes_palette("Moonrise2")[1]) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   scale_x_continuous(expand = c(0,0), breaks = function(limits) pretty(limits, n = 10, integer = TRUE)) +  
#   scale_y_continuous(limits = c(0,0.25), expand = c(0,0), labels = function(y) paste0(y * 100)) +  
#   theme(axis.title.y = element_blank(),
#         axis.title.x = element_blank())
# 
# print(plot_aj_rsa)
# 
# plot_aj_ch <- survfit2(Surv(persontime_years, incident_tb, type = "mstate") ~ 1, data = kaplan_ch) |>
#   ggcuminc(linewidth = 0.5, color = wes_palette("Moonrise2")[2]) +
#   add_confidence_interval(fill = wes_palette("Moonrise2")[2]) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   scale_x_continuous(expand = c(0,0), breaks = function(limits) pretty(limits, n = 10, integer = TRUE)) +  
#   scale_y_continuous(limits = c(0,0.008), expand = c(0,0), labels = function(y) paste0(y * 100)) +  
#   theme(axis.title.y = element_blank(),
#         axis.title.x = element_blank())
# 
# plot_aj_ch
# 
# aj_both <- ggarrange(plot_aj_rsa, plot_aj_ch, ncol = 2) %>% 
#   annotate_figure(bottom = "Years after ART start",
#                   left = "Cumulative incidence of TB-Incidence (%)")
#   
# aj_both
# 
# ggsave(plot = aj_both, filename = "results/incidence/incidence_aj_sep.png", 
#        width = 16, height = 11, units = "cm")
