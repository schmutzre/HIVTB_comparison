library(tidyverse)
library(epitools)
library(wesanderson)

df <- readRDS("data_clean/art.rds") %>% 
  dplyr::select(id, presenting_tb, cohort, last_persontime, exitdate, gender, age_at_art_start, region, who_stage, cd4_group, outcome_tb) %>% 
  mutate(region = relevel(region, ref = "Europe/Northern America"),
         cohort = relevel(cohort, ref = "RSA")) %>% 
  filter(presenting_tb ==1)


#### iwhod figure 1 (mortality rate)

df_mortalityRate <- df %>%
  filter(last_persontime > 0) %>% 
  mutate(persontime_years = last_persontime / 365,
         cohort = as.factor(ifelse(cohort == "RSA", "South Africa", "Switzerland")),
         cohort = fct_relevel(cohort, "Switzerland")) %>%
  group_by(cohort, presenting_tb) %>% 
  summarise(sum_exit = sum(!is.na(exitdate)), 
            sum_person_years = sum(persontime_years)/1000)  %>% 
  mutate(pois = pois.exact(x = sum_exit, pt = sum_person_years, conf.level = 0.95)) %>% 
  filter(presenting_tb ==1)

cohort_colors <- wes_palette("Moonrise2")
names(cohort_colors) <- c("South Africa", "Switzerland")

mortality_rate <- df_mortalityRate %>%
  ggplot(aes(x = cohort, y = pois$rate)) +
  geom_point(aes(color = cohort), position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = pois$lower, ymax = pois$upper, color = cohort), 
                position = position_dodge(width = 0.5), width = 0.2, linewidth = 1.5) +
  scale_color_manual(values = cohort_colors) +
  theme_classic(base_size = 28) +
  theme(legend.position = "none",
        plot.title = element_text(size = 28), # Set title properties here
        plot.title.position = "plot",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  labs(x = NULL, y = NULL,
       title = expression(atop(paste(bold("B |"), " All-cause mortality"), "       (presenting with TB)")))+
  coord_cartesian(ylim = c(0,50))+
  scale_y_continuous(breaks = c(0,10,20,30,40, 50), expand = c(0,0)) +
  guides(color = guide_legend(title = NULL))


mortality_rate


