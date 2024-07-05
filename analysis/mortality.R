library(tidyverse)
library(epitools)
library(wesanderson)

df <- readRDS("data_clean/art.rds") %>% 
  dplyr::select(id, presenting_tb, cohort, fup_time, death, gender, age_at_art_start, region, who_stage, cd4_group, outcome_tb) %>% 
  mutate(region = relevel(region, ref = "Europe/Northern America"),
         cohort = relevel(cohort, ref = "RSA")) %>% 
  filter(presenting_tb ==1)

table(df$cohort)
#### iwhod figure 1 (mortality rate)

df_mortalityRate <- df %>%
  filter(fup_time > 0) %>% 
  mutate(persontime_years = fup_time / 365,
         cohort = as.factor(ifelse(cohort == "RSA", "South Africa", "Switzerland")),
         cohort = fct_relevel(cohort, "Switzerland")) %>%
  group_by(cohort, presenting_tb) %>% 
  summarise(sum_exit = sum(death == 1), 
            sum_person_years = sum(persontime_years)/1000)  %>% 
  mutate(pois = pois.exact(x = sum_exit, pt = sum_person_years, conf.level = 0.95)) %>% 
  filter(presenting_tb ==1)

cohort_colors <- wes_palette("Moonrise2")
names(cohort_colors) <- c("South Africa", "Switzerland")

# mortality_rate <- df_mortalityRate %>%
#   ggplot(aes(x = cohort, y = pois$rate)) +
#   geom_point(aes(color = cohort), position = position_dodge(width = 0.5), size = 3) +
#   geom_errorbar(aes(ymin = pois$lower, ymax = pois$upper, color = cohort), 
#                 position = position_dodge(width = 0.5), width = 0.2, linewidth = 1.5) +
#   scale_color_manual(values = cohort_colors) +
#   theme_classic(base_size = 28) +
#   theme(legend.position = "none",
#         plot.title = element_text(size = 28), # Set title properties here
#         plot.title.position = "plot",
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA),
#         legend.background = element_rect(fill = "transparent", color = NA)) +
#   labs(x = NULL, y = NULL,
#        title = expression(atop(paste(bold("B |"), " All-cause mortality"), "       (presenting with TB)")))+
#   coord_cartesian(ylim = c(0,50))+
#   scale_y_continuous(breaks = c(0,10,20,30,40, 50), expand = c(0,0)) +
#   guides(color = guide_legend(title = NULL))


mortality_rate <- df_mortalityRate %>%
  ggplot(aes(x = cohort, y = pois$rate)) +
  geom_point(aes(color = cohort), position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = pois$lower, ymax = pois$upper, color = cohort),
                position = position_dodge(width = 0.5), width = 0.2, linewidth = 1.5) +
  scale_color_manual(values = cohort_colors) +
theme_bw()+
  labs(x = NULL, y = "Mortality rate (per 1000 person-years)",
       title = "All-cause mortality (presenting with TB)")+
  coord_cartesian(ylim = c(0,50))+
  scale_y_continuous(breaks = c(0,10,20,30,40, 50), expand = c(0,0)) +
  guides(color = guide_legend(title = NULL))


mortality_rate

ggsave("results/mortality_rate.png", mortality_rate, width = 16, height = 12, units = "cm")

#### KM ####

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

# Load and filter data for CH cohort
df_ch <- readRDS("data_clean/art.rds") %>% 
  select(id, presenting_tb, cohort, fup_time, death, gender, age_at_art_start, region, who_stage, cd4_group, outcome_tb) %>% 
  filter(cohort == "CH")

# Fit survival model for CH cohort
fit_ch <- survfit(Surv(time = df_ch$fup_time, event = df_ch$death) ~ presenting_tb, data = df_ch)

# Plot survival curves for CH cohort
p_ch <- ggsurvplot(fit_ch, 
                   data = df_ch, 
                   pval = FALSE, 
                   conf.int = FALSE,
                   risk.table = FALSE,
                   risk.table.col = "strata",
                   xlab = "Time (days)",
                   ylab = "Survival probability",
                   title = "Kaplan-Meier Survival Curve for CH Cohort",
                   legend.title = "Presenting TB",
                   legend.labs = c("No", "Yes"))

p_ch <- p_ch$plot + coord_cartesian(ylim = c(0.8, 1),
                                    xlim = c(0,2500))
# Save the CH cohort plot
ggsave("results/mortality_survival_curve_ch.png", plot = p_ch, width = 16, height = 12, units = "cm")

# Load and filter data for RSA cohort
df_rsa <- readRDS("data_clean/art.rds") %>% 
  select(id, presenting_tb, cohort, fup_time, death, gender, age_at_art_start, region, who_stage, cd4_group, outcome_tb) %>% 
  filter(cohort == "RSA")

# Fit survival model for RSA cohort
fit_rsa <- survfit(Surv(time = df_rsa$fup_time, event = df_rsa$death) ~ presenting_tb, data = df_rsa)

# Plot survival curves for RSA cohort
p_rsa <- ggsurvplot(fit_rsa, 
                    data = df_rsa, 
                    pval = FALSE, 
                    conf.int = FALSE,
                    risk.table = FALSE,
                    risk.table.col = "strata",
                    xlab = "Time (days)",
                    ylab = "Survival probability",
                    title = "Kaplan-Meier Survival Curve for RSA Cohort",
                    legend.title = "Presenting TB",
                    legend.labs = c("No", "Yes"))
p_rsa <- p_rsa$plot + coord_cartesian(ylim = c(0.8, 1),
                                      xlim = c(0, 2500))


# Save the RSA cohort plot
ggsave("results/mortality_survival_curve_rsa.png", plot = p_rsa, width = 16, height = 12, units = "cm")

# Display the plots
print(p_ch)
print(p_rsa)
