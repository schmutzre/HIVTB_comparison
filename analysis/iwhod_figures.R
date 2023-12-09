##### set-up -------------------------------------------------------------------

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  ggplot2, # for plotting
  haven,
  lme4,
  lmerTest,
  mgcv,
  gridExtra, 
  wesanderson,
  ggpubr,
  survival,
  survminer,
  scam
)

source("utils/plot.R")

set.seed(123)

##### Figure 1a ----------------------------------------------------------------

## data ##

cd4_ch <- readRDS("data_clean/ch/cd4_ch.rds")%>% 
  mutate(id = as.factor(id),
         cohort = as.factor("CH"))

cd4_rsa <- readRDS("data_clean/rsa/cd4_rsa.rds") %>% 
  mutate(id = as.factor(id),
         cohort = as.factor("RSA"))

cd4 <- rbind(cd4_ch %>% dplyr::select(-timepoint),
             cd4_rsa %>% dplyr::select(-timepoint)) %>% 
  filter(time_diff >= -90 & time_diff <= 360) %>% 
  mutate(cd4_trans = sqrt(cd4))

## model / plot ##

# datasets #

cd4_npres <- cd4 %>% 
  filter(presenting_tb == 0)

cd4_pres <- cd4 %>% 
  filter(presenting_tb == 1)

# model definition #

m.cd4_npres <- scam(cd4 ~ s(time_diff, k = 3, by = cohort, bs = "cr"), 
                    data = cd4_npres)

m.cd4_pres <- gam(cd4 ~ s(time_diff, k = 3, by = cohort, bs = "cr"), 
                  data = cd4_pres)

pred_data <- expand.grid(cohort = c("RSA", "CH"), time_diff = seq(-90, 365, 1)) %>% arrange(cohort)

cd4_pred_npres <- predict(m.cd4_npres, type = "response", se.fit = T, newdata = pred_data) %>%
  data.frame() %>%
  unnest(cols = c(fit, se.fit)) %>%
  mutate(time_diff = pred_data$time_diff,
         cohort = rep(c("South Africa", "Switzerland"), each = nrow(.) / 2),
         upper = fit + qnorm(.975) * se.fit,
         lower = fit - qnorm(.975) * se.fit)

cd4_pred_pres <- predict(m.cd4_pres, type = "response", se.fit = T, newdata = pred_data) %>%
  data.frame() %>%
  unnest(cols = c(fit, se.fit)) %>%
  mutate(time_diff = pred_data$time_diff,
         cohort = rep(c("South Africa", "Switzerland"), each = nrow(.) / 2),
         upper = fit + qnorm(.975) * se.fit,
         lower = fit - qnorm(.975) * se.fit)

cd4_pred <- rbind(cd4_pred_npres %>% mutate(type = "Without prevalent TB"), 
                  cd4_pred_pres %>% mutate(type = "With prevalent TB"))

cd4_pred_pl <- cd4_pred %>%
  ggplot(mapping = aes(x = time_diff)) +
  facet_wrap(vars(type)) +
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, fill = cohort), alpha = 0.2) +
  geom_line(mapping = aes(y = fit, color = cohort)) +
  #geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,600), breaks = seq(0, 600, 100)) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)"),
       title = "1a | CD4 count showing immunological recovery after ART start") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 

cd4_pred_pl

ggsave(cd4_pred_pl, filename = "results/cd4_pred.png",
       width = 12.7, height = 12.7, unit = "cm")

##### Figure 1b ----------------------------------------------------------------

## data ##

km <- readRDS("data_clean/art.rds") %>% 
  filter(last_persontime >= 0) %>% 
  mutate(exit = case_when(is.na(exitdate) ~ 0, TRUE ~ 1),
         cohort = case_when(cohort == "RSA" ~ 0,
                            TRUE ~1))

## model / plot ##

SurvObj <- Surv(time = km$last_persontime, event = km$exit == 1)

fit <- survfit(SurvObj ~ cohort, data = km)

# Plotting
km_plot <- ggsurvplot(fit, data = km, conf.int = TRUE, 
                       fun = "pct", 
                       palette = wes_palette("Moonrise2"),  # Use your predefined colors
                       legend = "right",
                       legend.title = "",censor = F, risk.table = T)

kaplan_plot_tot <- km_plot$plot + 
  coord_cartesian(xlim=c(0,362), ylim = c(97,100.1)) +
  scale_x_continuous(breaks = seq(0,360,30), expand = c(0,0)) +
  scale_y_continuous(limits = c(97,100.2), breaks = seq(97,100,1), expand = c(0,0))+
  xlab("Days after ART start")+
  theme_classic(base_size = 10) +
  scale_color_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) + 
  scale_fill_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 10,
                                                                  face = "bold"), # Set title properties here
        plot.title.position = "plot") +
  labs(title = "1b | Probability of survival after ART start",
       caption = "Mean as lines, 95%-CI as ribbons")

kaplan_plot_tot

ggsave(kaplan_plot_tot, filename = "results/kaplan.png",
       width = 12.7, height = 12.7, unit = "cm")

##### Figure 1b (weighted) -----------------------------------------------------

library(ipw)

km_f <- km %>% 
  filter(!is.na(cd4_baseline))

weights <- ipwpoint(exposure = cohort, family = "binomial", link = "logit", 
                    denominator = ~ age_at_art_start + gender + cd4_baseline, data = km_f)

SurvObj2 <- Surv(time = km_f$last_persontime, event = km_f$exit == 1)

fit2 <- survfit(SurvObj2 ~ cohort, data = km_f, weights = weights$den.mod$weights)

km_plot2 <- ggsurvplot(fit2, data = km_f, conf.int = TRUE, 
                      fun = "pct", 
                      palette = wes_palette("Moonrise2"),
                      legend = "right",
                      legend.title = "",censor = F,risk.table = T)
km_plot2$table

kaplan_plot_tot2 <- km_plot2$plot + 
  coord_cartesian(xlim=c(0,362), ylim = c(90,100.1)) +
  scale_x_continuous(breaks = seq(0,360,30), expand = c(0,0)) +
  scale_y_continuous(limits = c(90,100.2), breaks = seq(90,100,5), expand = c(0,0))+
  xlab("Days after ART start")+
  theme_classic(base_size = 10) +
  scale_color_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) + 
  scale_fill_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 10,
                                                                  face = "bold"), # Set title properties here
        plot.title.position = "plot") +
  labs(title = "1b | Probability of survival after ART start",
       caption = "Mean as lines, 95%-CI as ribbons")

kaplan_plot_tot2

ggsave(kaplan_plot_tot2, filename = "results/kaplan_weighted.png",
       width = 12.7, height = 12.7, unit = "cm")

##### Join together ------------------------------------------------------------

arranged_plots <- ggpubr::ggarrange(cd4_pred_pl, kaplan_plot_tot2, ncol = 1, common.legend = F) 
arranged_plots
ggsave(arranged_plots, filename = "results/conference.png", 
       width = 12.7, height = 12.7, units = "cm")

