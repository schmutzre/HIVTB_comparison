#### Libraries -----------------------------------------------------------------

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, 
  ggplot2, 
  epitools,
  sjPlot,
  gridExtra,
  ggpubr,
  wesanderson,
  forcats
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
    cohort = fct_relevel(cohort, "RSA"),
    cd4_group = fct_relevel(cd4_group, "350+")) %>% 
  dplyr::select(id, agegroup, cohort, art_start_date, incident_tb, date_tb, last_persontime, persontime_years, cd4_group, rna_group) %>% 
  filter(persontime_years > 0)
  
df_manual <- df %>% # calculating the incident rates per cohort
  group_by(cohort) %>% 
  summarise(sum_incident_tb = sum(incident_tb == 1), 
            sum_person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = sum_incident_tb, pt = sum_person_years, conf.level = 0.95),
         rna_group = "Overall",
         cd4_group = "Overall")

print(irr <- df_manual$pois$rate[df_manual$cohort == "RSA"] / df_manual$pois$rate[df_manual$cohort == "CH"])

#### Incidence per baseline group ----------------------------------------------

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

ggsave(plot = plot_rna, filename = "results/incidence/incidence_rna.png", 
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

#### IRR compared to 350+ group ------------------------------------------------

pois_ch <- glm(incident_tb ~ cd4_group, 
               offset=log(persontime_years), 
               family="poisson", 
               data=df %>% filter(cohort == "CH",
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


reference <- tibble(term = rep("350+",2), estimate = rep(1,2), conf.low = NA, conf.high = NA, cohort = c("RSA", "CH"))

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

ggsave(plot = combined_plot, file = "results/incidence/irr.png", 
       width = 16, height = 11, units = "cm") 
