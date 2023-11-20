##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  lubridate, # for date handling
  rlang, # for '!!' inside the functions
  AER, # for dispersiontest
  ggplot2, # for plotting
  haven,
  survival,
  survminer,
  gridExtra,
  Epi,
  cmprsk,
  wesanderson,
  ggsurvfit
)

##### data import/preprocessing ------------------------------------------------

kaplan <- readRDS("data_clean/art_noTB.rds")  %>% 
  mutate(incident_tb = as.numeric(as.character(incident_tb)),
         persontime_years = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")/360),
           incident_tb == 0 ~ last_persontime/360),
         persontime_days = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")),
           incident_tb == 0 ~ last_persontime),
         persontime_death_days = case_when(
           !is.na(exitdate) ~ as.numeric(difftime(exitdate, art_start_date, units = "days")),
           is.na(exitdate) ~ last_persontime),
         cohort = fct_relevel(cohort, "RSA")) %>% 
  dplyr::select(id, cohort, art_start_date, incident_tb, date_tb, last_persontime, persontime_years, persontime_days, persontime_death_days, cd4_group, exitdate, last_fup_date, rna_group) %>% 
  filter(persontime_years > 0) %>% 
  mutate(event_type = as.factor(case_when(
    incident_tb == 1 ~ 1,
    !is.na(exitdate) ~ 2,
    TRUE ~ 0 # Loss to follow-up isnt considered a competing risk
  )))

kaplan_rsa <- kaplan %>% 
  filter(cohort == "RSA")

tabyl(kaplan_rsa$event_type)

kaplan_ch <- kaplan %>% 
  filter(cohort == "CH")

tabyl(kaplan_ch$event_type)

#### Kaplan Meier model (excluding competing risk [death]) ---------------------

p <- survfit2(Surv(persontime_days, incident_tb) ~ cohort, data = kaplan) |>
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_quantile(y_value = 0.6, color = "gray50", linewidth = 0.75) +
  theme_bw()+
  theme(legend.position = "none") +
  labs(x = "Days after ART start",
       y = "TB-free survival") +
  scale_x_continuous(expand = c(0,0))

print(p)

p_death <- survfit2(Surv(persontime_death_days, incident_tb) ~ cohort, data = kaplan) |>
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  theme_bw()+
  theme(legend.position = "none") +
  labs(x = "Days after ART start",
       y = "survival") +
  scale_x_continuous(expand = c(0,0))

print(p_death)

########### Aalen-Johansen model (including competing risk [death]) ------------

## Both in one plot ##

plot_aj <- survfit2(Surv(persontime_years, event_type, type = "mstate") ~ cohort, data = kaplan) |>
  ggcuminc(linewidth = 0.5, outcome = "1") +
  add_confidence_interval() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Years after ART start",
       y = "Cumulative incidence of TB (%)") +
  scale_x_continuous(expand = c(0,0), breaks = function(limits) pretty(limits, n = 10, integer = TRUE)) +  
  scale_y_continuous(expand = c(0,0), labels = function(y) paste0(y * 100)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2"))

plot_aj

ggsave(plot = plot_aj, filename = "results/incidence/incidence_aj.png", 
       width = 16, height = 11, units = "cm")

## log-scaled ## 

plot_aj_log <- survfit2(Surv(persontime_years, event_type, type = "mstate") ~ cohort, data = kaplan) |>
  ggcuminc(linewidth = 0.5) +
  add_confidence_interval() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Years after ART start",
       y = "Log-scaled cumulative incidence of TB-Incidence (%)") +
  scale_y_log10(expand = c(0,0), labels = function(y) paste0(y * 100))+
  scale_x_continuous(expand = c(0,0), breaks = function(limits) pretty(limits, n = 10, integer = TRUE)) +  
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) 

plot_aj_log

ggsave(plot = plot_aj_log, filename = "results/incidence/incidence_aj_log.png", 
       width = 16, height = 11, units = "cm")

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
