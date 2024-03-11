library(tidyverse)
library(epitools)
library(wesanderson)

df <- readRDS("data_clean/art.rds") %>% 
  dplyr::select(id, presenting_tb, cohort, last_persontime, exitdate, gender, age_at_art_start, region, who_stage, cd4_group) %>% 
  mutate(region = relevel(region, ref = "Europe/Northern America"),
         cohort = relevel(cohort, ref = "RSA"))


#### iwhod figure 1 (mortality rate)

df_mortalityRate <- df %>%
  filter(last_persontime > 0) %>% 
  mutate(persontime_years = last_persontime / 365,
         cohort = ifelse(cohort == "RSA", "South Africa", "Swizerland")) %>%
  group_by(cohort, presenting_tb) %>% 
  summarise(sum_exit = sum(!is.na(exitdate)), 
            sum_person_years = sum(persontime_years)/1000)  %>% 
  mutate(pois = pois.exact(x = sum_exit, pt = sum_person_years, conf.level = 0.95)) %>% 
  filter(presenting_tb ==1)

mortality_rate <- df_mortalityRate %>%
  ggplot(aes(y = cohort, x = pois$rate)) +
  geom_point(aes(color = cohort), position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(xmin = pois$lower, xmax = pois$upper, color = cohort), 
                position = position_dodge(width = 0.5), width = 0.2) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20), # Set title properties here
        plot.title.position = "plot",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  labs(x = NULL, y = NULL,
       title = "1b | All-cause Mortality rate for people\n       presenting with TB") + # Added line break here
  scale_y_discrete(labels = NULL) +
  coord_cartesian(xlim = c(0,50))+
  scale_x_continuous(breaks = c(0,10,20,30,40, 50), expand = c(0,0)) +
  guides(color = guide_legend(title = NULL))


mortality_rate


