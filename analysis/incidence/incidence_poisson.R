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

df <- readRDS("data_clean/art.rds") %>% 
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
pois_rsa <- glm(incident_tb ~ cd4_group, offset=log(persontime_years), family="poisson", data=df %>% filter(cohort == "RSA",
                                                                                                           !is.na(cd4_group) & cd4_group != "NA"))
irr_ch <- plot_model(pois_ch,
                     colors = wes_palette("Moonrise2")[2],
                     vline.color = "black",
                     show.values = TRUE, 
                     value.offset = .3) +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels=c("", "")) 
  
irr_rsa <- plot_model(pois_rsa,
                     colors = wes_palette("Moonrise2")[1],
                     vline.color = "black",
                     show.values = TRUE, 
                     value.offset = .3) +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels=c("0-99", "100-349")) 

irr_both <- ggarrange(irr_rsa, irr_ch, ncol = 2) %>% 
  annotate_figure(bottom = "IRR",
                  left = "CD4 count at ART start")
  
print(irr_both)

ggsave(plot = irr_both, file = "results/incidence/irr_350.png", 
       width = 16, height = 11, units = "cm") 
#### Time to TB ----------------------------------------------------------------

# time_to_tb <- ch %>% 
#   filter(incidence ==1) %>% 
#   mutate(time_to_tb = as.numeric(date_tb - (art_start_date))) %>% 
#   summarise(median = median(time_to_tb, na.rm = TRUE),
#             "25%" = quantile(time_to_tb, 0.25, na.rm = TRUE),
#             "75%" = quantile(time_to_tb, 0.75, na.rm = TRUE))
# 
# time_to_tb