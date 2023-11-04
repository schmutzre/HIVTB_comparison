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

##### data import/preprocessing ----

kaplan <- readRDS("data_clean/art.rds") %>% 
  mutate(incident_tb = as.numeric(as.character(incident_tb)),
         persontime_years = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")/360),
           incident_tb == 0 ~ last_persontime/360),
         persontime_days = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")),
           incident_tb == 0 ~ last_persontime),
         persontime_death_days = case_when(
           !is.na(exitdate) ~ as.numeric(difftime(exitdate, art_start_date, units = "days")),
           is.na(exitdate) ~ last_persontime)) %>% 
  dplyr::select(id, cohort, art_start_date, incident_tb, date_tb, last_persontime, persontime_years, persontime_days, persontime_death_days, cd4_group, exitdate, last_fup_date, rna_group) %>% 
  filter(persontime_years > 0) %>% 
  mutate(event_type = as.factor(case_when(
    incident_tb == 1 ~ 1,
    !is.na(exitdate) ~ 2,
    TRUE ~ 0 # Loss to follow-up isnt considered a competing risk
  )))

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

## Method 1 ##

library(cmprsk)

print(cmprisk <- cuminc(ftime = kaplan$persontime_years, 
       fstatus = kaplan$event_type, 
       group = kaplan$cohort,
       cencode=0))

plot(cmprisk)

ggcompetingrisks(
  cmprisk,
  gnames = NULL,
  gsep = " ",
  multiple_panels = TRUE,
  ggtheme = theme_minimal(),
  coef = 1.96,
  conf.int = TRUE,
) +
  labs(title = NULL) 

## Method 2 ## <-- graphics easier to modify here 

plot_aj <- survfit2(Surv(persontime_years, incident_tb, type = "mstate") ~ cohort, data = kaplan) |>
  ggcuminc(linewidth = 0.5) +
  add_confidence_interval() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Years after ART start",
       y = "Cumulative incidence of TB-Incidence (%)") +
  scale_x_continuous(expand = c(0,0), breaks = function(limits) pretty(limits, n = 10, integer = TRUE)) +  # 'pretty' generates nice breaks for integers
  scale_y_continuous(expand = c(0,0), labels = function(y) paste0(y * 100)) +  # multiply by 100 for percentages
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))  # Optional: adjust the margin if needed

plot_aj

ggsave(plot = plot_aj, filename = "results/incidence/incidence_aj.png", 
       width = 16, height = 11, units = "cm")

#'Caption: 

#' The cumulative incidence rates for time to first tuberculosis incidence after ART initiation 
#' were calculated using Aalen-Johansen method, adjusting for all-cause death as a competing risk.
