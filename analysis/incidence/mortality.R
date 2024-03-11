##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  lubridate, # for date handling
  rlang, # for '!!' inside the functions
  AER, # for dispersiontest
  ggplot2, # for plotting
  haven,
  epitools,
  sjPlot,
  gridExtra,
  ggpubr,
  wesanderson,
  forcats
)

#### Data prep -----------------------------------------------------------------

df <- readRDS("data_clean/art.rds") %>% 
  mutate(exit = as.numeric(case_when(!is.na(exitdate) ~ 1,
                          TRUE ~ 0)),
         cohort = fct_relevel(cohort, "RSA"),
         persontime_years = last_persontime/360) %>% 
  dplyr::select(id, cohort, art_start_date, exit, persontime_years) %>% 
  filter(persontime_years > 0)


df_manual <- df %>% 
  group_by(cohort) %>% 
  summarise(sum_exit = sum(exit == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_exit, pt = sum_person_years, conf.level = 0.95))

