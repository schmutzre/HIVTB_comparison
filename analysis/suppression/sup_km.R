##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  ggplot2, # for plotting
  survival,
  survminer,
  gridExtra,
  Epi
)

##### data import/preprocessing ------------------------------------------------

rna_ch <- readRDS("data_clean/ch/rna_ch.rds")
rna_rsa <- readRDS("data_clean/rsa/rna_rsa.rds")
rna <- rbind(rna_ch, rna_rsa)

patients_with_valid_rna <- rna %>% # patients without any labdata shouldnt be treated as "followed"
  filter(!is.na(rna)) %>%
  distinct(id)

rna400 <- rna %>% 
  dplyr::select(id, rna,time_diff) %>% 
  filter(time_diff >= 0 & rna < 400) %>%
  arrange(id, time_diff) %>% 
  group_by(id) %>%
  filter(time_diff == min(time_diff)) %>%  # Select the rows with the smallest time_diff
  ungroup()

rna200 <- rna %>% 
  dplyr::select(id, rna,time_diff) %>% 
  filter(time_diff >= 0 & rna < 200) %>%
  arrange(id, time_diff) %>% 
  group_by(id) %>%
  filter(time_diff == min(time_diff)) %>%  # Select the rows with the smallest time_diff
  ungroup()

combined <- readRDS("data_clean/art.rds")

aj_sup400 <- patients_with_valid_rna %>% 
  left_join(rna400, by = "id") %>% 
  left_join(combined %>% dplyr::select(id, presenting_tb, cohort, last_persontime, exitdate), by = "id") %>% 
  mutate(time = case_when(!is.na(time_diff) ~ time_diff,
                          TRUE ~ last_persontime),
         event = as.factor(case_when(!is.na(time_diff) ~ 1,
                           !is.na(exitdate) ~ 2,
                           TRUE ~ 0))) %>% 
  filter(time >= 0) %>% 
  dplyr::select(id, presenting_tb, time, rna, event, cohort)

aj_sup200 <- patients_with_valid_rna %>% 
  left_join(rna200, by = "id") %>% 
  left_join(combined %>% dplyr::select(id, presenting_tb, cohort, last_persontime, exitdate), by = "id") %>% 
  mutate(time = case_when(!is.na(time_diff) ~ time_diff,
                          TRUE ~ last_persontime),
         event = as.factor(case_when(!is.na(time_diff) ~ 1,
                                     !is.na(exitdate) ~ 2,
                                     TRUE ~ 0))) %>% 
  filter(time >= 0) %>% 
  dplyr::select(id, presenting_tb, time, rna, event, cohort)

##### Aalen-Johansen model -----------------------------------------------------

options("ggsurvfit.switch-color-linetype" = TRUE)

# 400 #

aj_rsa400 <- aj_sup400 %>% filter(cohort == "RSA")
aj_ch400 <- aj_sup400 %>% filter(cohort == "CH")

plot_aj_rsa400 <- survfit2(Surv(time, event) ~ presenting_tb, data = aj_rsa400) |>
  ggcuminc(outcome = "1", 
           aes(linetype = presenting_tb), 
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[1]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[1]) +
  theme_classic() +
  coord_cartesian(xlim = c(0,2000)) +
  theme(axis.title.y = element_blank(),
                axis.title.x = element_blank())+
  scale_ggsurvfit()+
  theme(legend.position = "none")

plot_aj_rsa400

plot_aj_ch400 <- survfit2(Surv(time, event) ~ presenting_tb, data = aj_ch400) |>
  ggcuminc(outcome = "1", 
           aes(linetype = presenting_tb), 
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[2]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[2]) +
  theme_classic() +
  coord_cartesian(xlim = c(0,2000)) +
  theme(axis.title.y = element_blank(),
                 axis.title.x = element_blank()) +
  scale_ggsurvfit() +
  theme(legend.position = "none")

plot_aj_ch400

aj_sup_both400 <- ggarrange(plot_aj_rsa400, plot_aj_ch400, ncol = 2) %>% 
     annotate_figure(bottom = "Days after ART start",
                     left = "Cumulative occurence of viral suppression (%)")
aj_sup_both400

ggsave(plot = aj_sup_both400, filename = "results/sup/sup_aj400.png", 
       width = 16, height = 11, units = "cm")

# 200 #

aj_rsa200 <- aj_sup200 %>% filter(cohort == "RSA")
aj_ch200 <- aj_sup200 %>% filter(cohort == "CH")

plot_aj_rsa200 <- survfit2(Surv(time, event) ~ presenting_tb, data = aj_rsa200) |>
  ggcuminc(outcome = "1", 
           aes(linetype = presenting_tb), 
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[1]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[1]) +
  theme_classic() +
  coord_cartesian(xlim = c(0,2000)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  scale_ggsurvfit()+
  theme(legend.position = "none")

plot_aj_rsa200

plot_aj_ch200 <- survfit2(Surv(time, event) ~ presenting_tb, data = aj_ch200) |>
  ggcuminc(outcome = "1", 
           aes(linetype = presenting_tb), 
           linewidth = 0.7, 
           color = wes_palette("Moonrise2")[2]) +
  add_confidence_interval(fill = wes_palette("Moonrise2")[2]) +
  theme_classic() +
  coord_cartesian(xlim = c(0,2000)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_ggsurvfit() +
  theme(legend.position = "none")

plot_aj_ch200

aj_sup_both200 <- ggarrange(plot_aj_rsa200, plot_aj_ch200, ncol = 2) %>% 
  annotate_figure(bottom = "Days after ART start",
                  left = "Cumulative occurence of viral suppression (%)")
aj_sup_both200

ggsave(plot = aj_sup_both200, filename = "results/sup/sup_aj200.png", 
       width = 16, height = 11, units = "cm")
