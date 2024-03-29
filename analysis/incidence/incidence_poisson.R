#### libraries -----------------------------------------------------------------

library(tidyverse)
library(epitools)
library(fmsb)
library(wesanderson)
library(ggpubr)
library(forcats)
library(broom)

#### data preparation ----------------------------------------------------------

custom_breaks <- c(16, 34, 44, 100)

df <- readRDS("data_clean/art_noTB.rds") %>% 
  mutate(incident_tb = as.numeric(as.character(incident_tb)),
         persontime_years = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")/360),
           incident_tb == 0 ~ last_persontime/360),
         agegroup = cut(age_at_art_start, breaks = custom_breaks, labels = c("16-34","35-44", "45+"), include.lowest = TRUE),
         cohort = fct_relevel(cohort, "RSA"),
         cd4_group = fct_relevel(cd4_group, "350+"),
         pre_2016 = fct_relevel(pre_2016, "1")) %>% 
  dplyr::select(id, agegroup, cohort, art_start_date, incident_tb, date_tb, last_persontime, persontime_years, cd4_group, rna_group, pre_2016) %>% 
  filter(persontime_years > 0)

levels(df$pre_2016) <- c("pre", "post")
  
df_manual <- df %>% 
  group_by(cohort) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95),
         rna_group = "Overall",
         cd4_group = "Overall")

#### IRR [RSA / CH] ------------------------------------------------------------

a <-  df_manual %>% filter(cohort == "RSA") %>% pull(sum_incident_tb)
b <-  df_manual %>% filter(cohort == "CH") %>% pull(sum_incident_tb)
PT1 <- df_manual %>% filter(cohort == "RSA") %>% pull(sum_person_years)
PT0 <- df_manual %>% filter(cohort == "CH") %>% pull(sum_person_years)
irr <- rateratio(a, b, PT1, PT0, conf.level=0.95)

df_irr <- tibble(
  est = irr[["estimate"]],
  lwr = irr[["conf.int"]][1],
  uppr = irr[["conf.int"]][2],
  cohort = "RSA")

saveRDS(df_irr, "results/incidenceTB/irrCohort.rds")

df_incidenceRate <- df %>% 
  group_by(cohort) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95),
         cohort = ifelse(cohort == "RSA", "South Africa", "Swizerland")) 

incidence_rate <- df_incidenceRate %>% 
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
  labs( x = NULL, y = NULL,
       title = "1a | Incident TB rate") +
  coord_cartesian(xlim = c(0,10))+
  scale_x_continuous(expand = c(0,0), limits = c(0,10), n.breaks = 3) +
  scale_y_discrete(labels = NULL) +
  guides(color = guide_legend(title = NULL))

incidence_rate

ggsave(plot = incidence_rate, filename = "results/incidenceTB/rate.png", bg='transparent', height = 6.5, units = "cm")

#### incidence rate per baseline group -----------------------------------------

### RNA ###

df_rna <- df %>% 
  mutate(rna_group = case_when(is.na(rna_group) | rna_group == "NA" ~ "NA",
                               rna_group %in% c("0-999", "1000-9999") ~ "< 10^3",
                               TRUE ~ ">= 10^3")) %>% 
  group_by(cohort, rna_group) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95)) %>%  
  bind_rows(df_manual)

plot_rna <- df_rna %>% 
  ggplot(aes(x = rna_group, y = pois$rate, group = interaction(cohort, rna_group))) +
  geom_point(aes(color = cohort), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = pois$lower, ymax = pois$upper, color = cohort), 
                position = position_dodge(width = 0.5), width = 0.2) +
  labs(x = "Baseline HIV RNA viral load", y = "Incident TB rate per 1,000 person-years") +
  scale_y_continuous() +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_x_discrete(labels = c("NA" = "NA", "< 10^3" = expression("< 10"^3), ">= 10^3" = expression("\u2265 10"^3))) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~cohort, scales = "free_y") 

print(plot_rna)

ggsave(plot = plot_rna, filename = "results/incidenceTB/rnaGrouped.png", 
       width = 16, height = 11, units = "cm")

### CD4 ###

df_cd4 <- df %>% 
  mutate(cd4_group = case_when(is.na(cd4_group) | cd4_group == "NA" ~ "NA",
                               TRUE ~ cd4_group)) %>% 
  group_by(cohort, cd4_group) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95)) %>% 
  bind_rows(df_manual)

plot_cd4 <- df_cd4 %>% 
  ggplot(aes(x = cd4_group, y = pois$rate, group = interaction(cohort, cd4_group))) +
  geom_point(aes(color = cohort), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = pois$lower, ymax = pois$upper, color = cohort), 
                position = position_dodge(width = 0.5), width = 0.2) +
  labs(x = "Baseline CD4 count", y = "Incident TB rate per 1,000 person-years") +
  scale_y_continuous() +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~cohort, scales = "free_y")

print(plot_cd4)

ggsave(plot = plot_cd4, filename = "results/incidenceTB/cd4Grouped.png", 
       width = 16, height = 11, units = "cm")

## Combined ##

plot_rna <- plot_rna +
  ylab("")

plot_cd4 <- plot_cd4 +
  ylab("")

# Create a combined plot with a single y-axis label
plot_both <- ggarrange(plot_cd4, plot_rna, # Add labels to identify the plots
                       ncol = 1, nrow = 2, # Arrange the plots in 2 rows
                       common.legend = TRUE, # Use a common legend for both plots
                       legend = "right") # Place the legend at the bottom

# Add a common y-axis label
plot_both <- annotate_figure(plot_both, 
                             left = text_grob("Incident TB rate per 1000 person-years", 
                                              rot = 90, vjust = 1))

plot_both

# Save the combined plot
ggsave(plot = plot_both, file = "results/incidenceTB/bothGrouped.png", 
       width = 16, height = 11, units = "cm") 

#### incidence rate per pre/post 2016 ------------------------------------------

### RNA ###

df_rna2k16 <- df %>% 
  mutate(rna_group = case_when(is.na(rna_group) | rna_group == "NA" ~ "NA",
                               rna_group %in% c("0-999", "1000-9999") ~ "< 10^3",
                               TRUE ~ ">= 10^3")) %>% 
  group_by(cohort, rna_group, pre_2016) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95)) 

plot_rna2k16 <- df_rna2k16 %>% 
  ggplot(aes(x = rna_group, y = pois$rate, shape = as.factor(pre_2016))) +
  geom_point(aes(color = cohort), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = pois$lower, ymax = pois$upper, color = cohort), 
                position = position_dodge(width = 0.5), width = 0.2) +
  labs(x = "Baseline CD4 count", y = "Incident TB rate per 1,000 person-years") +
  scale_y_continuous() +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  theme_classic() +
  scale_x_discrete(labels = c("NA" = "NA", "< 10^3" = expression("< 10"^3), ">= 10^3" = expression("\u2265 10"^3))) +
  facet_wrap(~cohort, scales = "free_y") 

plot_rna2k16

### CD4 ###

df_cd42k16 <- df %>% 
  mutate(cd4_group = case_when(is.na(cd4_group) | cd4_group == "NA" ~ "NA",
                               TRUE ~ cd4_group)) %>% 
  group_by(cohort, cd4_group, pre_2016) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95)) 

plot_cd42k16 <- df_cd42k16 %>% 
  ggplot(aes(x = cd4_group, y = pois$rate, shape = as.factor(pre_2016))) +
  geom_point(aes(color = cohort), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = pois$lower, ymax = pois$upper, color = cohort), 
                position = position_dodge(width = 0.5), width = 0.2) +
  labs(x = "Baseline CD4 count", y = "Incident TB rate per 1,000 person-years") +
  scale_y_continuous() +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  theme_classic() +
  facet_wrap(~cohort, scales = "free_y")

plot_cd42k16

ggsave(plot = plot_cd42k16, 
       file = "results/incidenceTB/2k16.png", 
       width = 16, height = 11, units = "cm") 

#### IRR compared to 350+ group ------------------------------------------------

pois_ch <- glm(incident_tb ~ cd4_group, 
               offset=log(persontime_years), 
               family="poisson", 
               data = df %>% filter(cohort == "CH",
                                  !is.na(cd4_group) & cd4_group != "NA"))

tidy_pois_ch <- tidy(pois_ch, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(cohort = as.factor("CH"))

pois_rsa <- glm(incident_tb ~ cd4_group, offset=log(persontime_years), 
                family="poisson", 
                data = df %>% 
                  filter(cohort == "RSA",
                         !is.na(cd4_group) & cd4_group != "NA"))

tidy_pois_rsa <- tidy(pois_rsa, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(cohort = as.factor("RSA"))

combined_data <- rbind(tidy_pois_rsa, tidy_pois_ch) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(cohort = factor(cohort, levels = c("RSA", "CH")),
         term = case_when(term == "cd4_group0-99" ~ "0-99",
                          TRUE ~ "100-349")) %>% 
  dplyr::select(term, estimate, conf.low, conf.high, cohort)

reference <- tibble(term = rep("350+", 2), 
                    estimate = rep(1,2), 
                    conf.low = NA, 
                    conf.high = NA, 
                    cohort = c("RSA", "CH"))

combined_data <- rbind(combined_data, reference)

combined_plot <- ggplot(combined_data, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = cohort)) +
  geom_point(position = position_dodge(width = -0.5)) +
  geom_errorbar(aes(group = cohort), width = 0.2, position = position_dodge(width = -0.5)) +
  theme_bw() +
  scale_y_log10() +  # Apply logarithmic scale to y-axis
  labs(x = "CD4 count at ART start", y = "Incident Rate Ratios (IRR)") +
  scale_color_manual(values = wes_palette("Moonrise2")[1:2]) +
  scale_x_discrete(labels = c("350+" = "\u2265 350"))+
  coord_flip() +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4)

combined_plot

ggsave(plot = combined_plot, file = "results/incidenceTB/irr.png", 
       width = 16, height = 11, units = "cm") 
