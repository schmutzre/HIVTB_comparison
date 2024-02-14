#### Data ----------------------------------------------------------------------

both <- readRDS("data_clean/art.rds")
rsa <- readRDS("data_clean/rsa/art_rsa.rds")
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(biostatUZH)

#### Prep ----------------------------------------------------------------------

interval <- both %>% 
  mutate(interval = case_when(last_persontime < 30 ~ "0-3",
                              last_persontime < 180 ~ "3-6",
                              last_persontime < 270 ~ "6-9",
                              last_persontime < 360 ~ "9-12",
                              TRUE ~ "12+"),
         exit = case_when(!is.na(exitdate) ~ 1, 
                          TRUE ~ 0),
         drop = case_when(is.na(exitdate) ~ 1, 
                          TRUE ~ 0)) %>% 
  filter(last_persontime >= 0) %>% 
  dplyr::select(id, cohort, art_start_date, exitdate, presenting_tb, last_persontime, interval, exit, drop)

interval$interval <- factor(interval$interval, levels = c("0-3", "3-6", "6-9", "9-12"))

rsa_interval <- interval %>% filter(cohort == "RSA")

rsa_exit <- rsa_interval %>% 
  filter(!is.na(exitdate)) 

#### Plots ---------------------------------------------------------------------

## Absolute counts ----

exit_count <- rsa_exit %>% 
  ggplot()+
  geom_bar(aes(x = interval, fill = presenting_tb), color = "black") +
  theme_classic() +
  labs(x = "Months after ATR start",
       y = "Number of deaths") +
  scale_y_continuous(expand=c(0,0), limits = c(0,200)) +
  facet_wrap(~ presenting_tb, labeller = labeller(presenting_tb = c("0" = "not presenting with TB", 
                                                                    "1" = "presenting with TB")))+
  scale_fill_manual(values = wes_palette("Rushmore1"))+
  theme(legend.position = "none") 

ggsave(plot = exit_count, filename = "results/exit_count.png", 
         width = 16, height = 11, units = "cm")

## Relative counts / rel. to total number of deaths per group ----

n_pres <- rsa_interval %>% 
  filter(presenting_tb == 1) %>% 
  nrow()

n_no_pres <- rsa_interval %>% 
  filter(presenting_tb == 0) %>% 
  nrow()

exit_prop <- rsa_interval %>%
  filter(!is.na(interval)) %>% 
  group_by(presenting_tb, interval) %>%
  summarize(exit = sum(exit),
            drop = sum(drop)) %>% 
  mutate(out = exit + drop) %>% 
  ungroup() %>% 
  group_by(presenting_tb) %>%
  mutate(total_out = cumsum(out),
         total_in = sum(out) - total_out,
         total_in_before = lag(total_in)) %>% 
  ungroup() %>% 
  mutate(total_in_before = case_when(presenting_tb == 0 & is.na(total_in_before) ~ n_no_pres,
                                     presenting_tb == 1 & is.na(total_in_before) ~ n_pres,
                                     TRUE ~ total_in_before)) %>%
  rowwise() %>%
  mutate(ci = list(wilson(exit, total_in_before))) %>% 
  ungroup() %>%
  dplyr::select(presenting_tb, interval, exit, drop, total_in_before, ci) %>%
  unnest_wider(ci, names_sep = "_")
  
# Adding the plot
exit_prop_plot <- ggplot(exit_prop, aes(x = interval, y = ci_prop * 100, color = factor(presenting_tb), group = factor(presenting_tb))) +
  geom_point() +
  geom_line() +  
  geom_errorbar(aes(ymin = ci_lower * 100, ymax = ci_upper * 100), width = 0.2) +
  theme_classic() +
  scale_color_manual(values = wes_palette("BottleRocket1")) +
  theme(legend.position = "bottom") +
  labs(y = "Mortality rate (%)",
       x = "Months after ART start",
       color = "Presenting TB") +
  scale_y_log10(labels = function(x) paste0(format(round(x, 1))))

exit_prop_plot

ggsave(plot = exit_prop_plot, filename = "results/exit_prop.png", 
       width = 16, height = 11, units = "cm")

## Stacked density plot ----
facet_labels <- c("CH.0" = "CH - Not Presenting TB",
                  "RSA.0" = "RSA - Not Presenting TB",
                  "CH.1" = "CH - Presenting TB",
                  "RSA.1" = "RSA - Presenting TB")

fup_prop <- interval %>%
  group_by(presenting_tb, last_persontime, cohort) %>% 
  summarize(exit = sum(exit),
            drop = sum(drop)) %>% 
  mutate(out = exit + drop) %>% 
  ungroup() %>% 
  group_by(presenting_tb, cohort) %>% 
  mutate(total_exit = cumsum(exit),
         total_drop = cumsum(drop),
         total_out = cumsum(out),
         total_alive = sum(out) - total_out,
         total_alive_before = lag(total_alive)) %>% 
  ungroup() %>% 
  mutate(sum = total_exit + total_drop + total_alive) %>% 
  dplyr::select(presenting_tb, last_persontime, total_exit, total_drop, total_alive, sum, cohort) %>% 
  pivot_longer(cols = c("total_alive", "total_drop", "total_exit"), 
               names_to = "category", values_to = "value") %>% 
  group_by(last_persontime, presenting_tb, cohort) %>%
  mutate(proportion = value / sum(value)) %>% 
  ungroup() %>% 
  mutate(cohort = factor(cohort, levels = c("CH", "RSA")),
         presenting_tb = factor(presenting_tb, levels = c("0", "1")),
         facet_row = factor(cohort, levels = c("CH", "RSA")),
         facet_col = factor(presenting_tb, levels = c("0", "1"), labels = c("Not Presenting TB", "Presenting TB"))) %>% 
  ggplot(aes(x = last_persontime, y = proportion, fill = category)) +
  geom_area(position = 'stack', alpha = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  facet_grid(facet_row ~ facet_col) +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x = "Days after ART start") +
  scale_fill_manual(values = wes_palette("BottleRocket2"),
                    labels = c("Followed", "Lost to follow-up", "Dead")) +
  coord_cartesian(xlim = c(0,800))

ggsave(plot = fup_prop, filename = "results/fup_prop.png", 
       width = 16, height = 11, units = "cm")

fup_prop

## fup plot from first eda: 

# rsa <- rsa %>% 
#   mutate(exit = case_when(!is.na(exitdate) ~ 1,
#                           TRUE ~ 0))
# 
# ggplot(data = rsa, aes(x = last_persontime, fill = factor(exit))) +
#   geom_histogram(binwidth = 7, color = "black", alpha = 0.7) +
#   labs(
#     x = "last_persontime",
#     y = "Frequency",
#     fill = "exit"
#   )+
#   scale_x_continuous(limits = c(0,360), expand = c(0,0)) +
#   scale_y_continuous(limits = c(0,1750), expand = c(0,0)) +
#   theme_bw()
