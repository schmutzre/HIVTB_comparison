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
  wesanderson
  )

#### Data prep -----------------------------------------------------------------

custom_breaks <- c(16, 24, 34, 44, 100)

df <- readRDS("data_clean/art_noTB.rds") %>% 
  mutate(incident_tb = as.numeric(as.character(incident_tb)),
    agegroup = cut(age_at_art_start, breaks = custom_breaks, include.lowest = TRUE),
         agegroup = as.factor(agegroup),
         persontime_years = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")/360),
           incident_tb == 0 ~ last_persontime/360),
    cohort = fct_relevel(cohort, "RSA")) %>% 
  dplyr::select(id, agegroup, cohort, art_start_date, incident_tb, date_tb, last_persontime, persontime_years, cd4_group, rna_group) %>% 
  filter(persontime_years > 0)
  
#### Model ---------------------------------------------------------------------

### Poisson ###

## manual calculation ##
df_manual <- df %>% 
  group_by(cohort) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95),
         rna_group = "Overall",
         cd4_group = "Overall")

print(irr <- df_manual$pois$rate[df_manual$cohort == "RSA"] / df_manual$pois$rate[df_manual$cohort == "CH"])

## fitting poisson model ##
m_pois <- glm(incident_tb ~ cohort, offset=log(persontime_years), family="poisson", data=df)

# Calculate exponentiated coefficients (IRR) along with CI
exp(coef(summary(m_pois)))

# CI
exp(confint(m_pois))

#### Assumptions ---------------------------------------------------------------

#' checking for overdispersion - the Poisson model assumes that the mean and variance of the outcome are equal. 
#' If the variance is greater than the mean, the data are overdispersed.

disp_test <- dispersiontest(m_pois, trafo=1)  # 'trafo = 1' for log-transformed data
print(disp_test)

#### Plots ---------------------------------------------------------------------

###### IRR ######

plot_IRR <- plot_model(m_pois,
                       colors = "blue",
                       vline.color = "red",
                       show.values = TRUE, 
                       value.offset = .3) +
  theme_bw() +
  theme(plot.title = element_blank())

plot_IRR

ggsave(filename = "results/incidence/incidence_IRR.png", plot = plot_IRR, width = 10, height = 8, dpi = 300)

###### Incidence per baseline group ######

## RNA ##

df_rna <- df %>% 
  mutate(rna_group = case_when(is.na(rna_group) | rna_group == "NA" ~ "NA",
                               TRUE ~ rna_group)) %>% 
  group_by(cohort, rna_group) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95))  

df_combined <- bind_rows(df_rna, df_manual)

overall_rsa <- df_combined %>% filter(cohort == "RSA",
                                      rna_group == "Overall")
rsa_rate <- overall_rsa$pois$rate
rsa_lrate <- overall_rsa$pois$lower
rsa_urate <- overall_rsa$pois$upper

overall_ch <- df_combined %>% filter(cohort == "CH",
                                      rna_group == "Overall")

ch_rate <- overall_ch$pois$rate
ch_lrate <- overall_ch$pois$lower
ch_urate <- overall_ch$pois$upper

plot_rna <- df_combined %>% 
  ggplot(aes(x = rna_group, y = pois$rate, group = interaction(cohort, rna_group))) +
  geom_point(aes(color = cohort), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = pois$lower, ymax = pois$upper, color = cohort), 
                position = position_dodge(width = 0.5), width = 0.2) +
  labs(x = "Baseline HIV RNA viral load", y = "Incident TB rate per 1,000 person-years") +
  scale_y_continuous() +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~cohort, scales = "free_y")

print(plot_rna)

ggsave(plot = plot_rna, filename = "results/incidence/incidence_rna.png", 
       width = 16, height = 11, units = "cm")

## CD4 ##

df_cd4 <- df %>% 
  mutate(cd4_group = case_when(is.na(cd4_group) | cd4_group == "NA" ~ "NA",
                               TRUE ~ cd4_group)) %>% 
  group_by(cohort, cd4_group) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95))

df_combined <- bind_rows(df_cd4, df_manual)

plot_cd4 <- df_combined %>% 
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

ggsave(plot = plot_cd4, filename = "results/incidence/incidence_cd4.png", 
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
ggsave(plot = plot_both, file = "results/incidence/incidence_both.png", 
       width = 16, height = 11, units = "cm") 

###### IRR compared to 350+ group ######

df <- df %>%
  mutate(cd4_group = relevel(cd4_group, ref = "350+"))

pois_ch <- glm(incident_tb ~ cd4_group, offset=log(persontime_years), family="poisson", data=df %>% filter(cohort == "CH",
                                                                                                           !is.na(cd4_group) & cd4_group != "NA"))
tidy_pois_ch <- broom::tidy(pois_ch, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(cohort = as.factor("CH"))

pois_rsa <- glm(incident_tb ~ cd4_group, offset=log(persontime_years), family="poisson", data=df %>% filter(cohort == "RSA",
                                                                                                           !is.na(cd4_group) & cd4_group != "NA"))
tidy_pois_rsa <- broom::tidy(pois_rsa, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(cohort = as.factor("RSA"))

combined_data <- rbind(tidy_pois_rsa, tidy_pois_ch) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(cohort = factor(cohort, levels = c("RSA", "CH")),
         term = case_when(term == "cd4_group0-99" ~ "0-99",
                          TRUE ~ "100-349")) %>% 
  dplyr::select(term, estimate, conf.low, conf.high, cohort)

reference <- tibble(term = rep("350+",2), estimate = rep(1,2), conf.low = rep(1,2), conf.high = rep(1,2), cohort = c("RSA", "CH"))
combined_data <- rbind(combined_data, reference)


combined_plot <- ggplot(combined_data, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = cohort)) +
  geom_point(position = position_dodge(width = -0.5)) +
  geom_errorbar(aes(group = cohort), width = 0.2, position = position_dodge(width = -0.5)) +
  theme_bw() +
  scale_y_log10() +  # Apply logarithmic scale to y-axis
  labs(x = "CD4 count at ART start", y = "Incident Rate Ratios (IRR)") +
  scale_color_manual(values = wes_palette("Moonrise2")[1:2]) +
  coord_flip() +
  theme(legend.title = element_blank())

combined_plot

ggsave(plot = combined_plot, file = "results/incidence/irr.png", 
       width = 16, height = 11, units = "cm") 

#### Time to TB ----------------------------------------------------------------

# time_to_tb <- ch %>% 
#   filter(incidence == 1) %>% 
#   mutate(time_to_tb = as.numeric(date_tb - (art_start_date))) %>% 
#   summarise(median = median(time_to_tb, na.rm = TRUE),
#             "25%" = quantile(time_to_tb, 0.25, na.rm = TRUE),
#             "75%" = quantile(time_to_tb, 0.75, na.rm = TRUE))
# 
# time_to_tb