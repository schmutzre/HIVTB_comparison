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
  survminer
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
  filter(time_diff >= -60 & time_diff < 365) %>% 
  mutate(cd4_trans = sqrt(cd4))

## model / plot ##

# datasets #

cd4_npres_rsa <- cd4 %>% 
  filter(presenting_tb == 0,
         cohort == "RSA") 

cd4_npres_ch <- cd4 %>% 
  filter(presenting_tb == 0,
         cohort == "CH")

cd4_pres_rsa <- cd4 %>% 
  filter(presenting_tb == 1,
         cohort == "RSA")

cd4_pres_ch <- cd4 %>% 
  filter(presenting_tb == 1,
         cohort == "CH")

# model definition #

m.cd4_rsa_npres <- gam(cd4 ~ s(time_diff, bs = "ad", k = 8),
                       data = cd4_npres_rsa,
                       method = "REML")
gam

m.cd4_ch_npres <- gam(cd4 ~ s(time_diff, k = 6, bs = "cr"), 
                      data = cd4_npres_ch, 
                      method = "REML")

m.cd4_rsa_pres <- gam(cd4 ~ s(time_diff, k = 6, bs = "cr"), 
                      data = cd4_pres_rsa, 
                      method = "REML")

m.cd4_ch_pres <- gam(cd4 ~ s(time_diff, k = 6, bs = "cr"), 
                     data = cd4_pres_ch, 
                     method = "REML")

# model prediction #

cd4_rsa_np <- pred_trend_single(m.cd4_rsa_npres, cd4_npres_rsa) %>% 
  mutate(tb = factor("Not presenting TB"),
         cohort = factor("RSA"))

cd4_ch_np <- pred_trend_single(m.cd4_ch_npres, cd4_npres_ch) %>% 
  mutate(tb = factor("Not presenting TB"),
         cohort = factor("CH"))

cd4_rsa_p <- pred_trend_single(m.cd4_rsa_pres, cd4_pres_rsa) %>% 
  mutate(tb = factor("Presenting TB"),
         cohort = factor("RSA"))  

cd4_ch_p <- pred_trend_single(m.cd4_ch_pres, cd4_pres_ch) %>% 
  mutate(tb = factor("Presenting TB"),
         cohort = factor("CH"))

# model plotting #

cd4_pred_single <- rbind(cd4_rsa_np, cd4_ch_np, cd4_rsa_p, cd4_ch_p)

cd4_pred1 <- ggplot(cd4_pred_single, aes(x = time_diff)) +
  geom_ribbon(aes(ymin = lwrS_adj, ymax = uprS, fill = cohort), alpha = 0.2) +
  geom_line(aes(y = fit, color = cohort)) +
  theme_classic(base_size = 10) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,700,100)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(-60,360,60)) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        plot.title = element_text(size = 10, face = "bold"), # Set title properties here
        plot.title.position = "plot") +
  scale_color_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) +
  scale_fill_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) +
  coord_cartesian(ylim = c(0,700), xlim = c(-60,380))+
  facet_wrap(vars(tb)) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)"),
       title = "a | CD4 count showing immunological recovery after ART start for patients presenting with/without TB")

cd4_pred1

# rohwerte rsa non-presenting #

#checking
rec.cd4 <- cd4_npres_rsa %>% 
  mutate(
    group_30days = as.integer((time_diff %/% 30) + 1)
  )

summary_stats.cd4 <- rec.cd4 %>% 
  group_by(group_30days) %>% 
  summarise(
    median.cd4 = median(cd4, na.rm = TRUE),
    Q1 = quantile(cd4, 0.25, na.rm = TRUE),
    Q3 = quantile(cd4, 0.75, na.rm = TRUE),
    time_diff_midpoint = mean(time_diff, na.rm = TRUE),
    Qlow = quantile(cd4, 0.05, na.rm = TRUE),
    Qhigh = quantile(cd4, 0.95, na.rm = TRUE),
    mean = mean(cd4, na.rm = TRUE),
    se = sd(cd4, na.rm = TRUE) / sqrt(n()),
    sd = sd(cd4, na.rm = TRUE)
  )

IQR.cd4 <- ggplot(summary_stats.cd4, aes(x = time_diff_midpoint, y = median.cd4)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Q1, ymax = Q3),
    width = 10
  ) +
  labs(x = 'Time Difference (days)', y = 'CD4') +
  theme_bw() +
  geom_hline(yintercept = 350)
IQR.cd4

##### Figure 1b ----------------------------------------------------------------

## data ##

km <- readRDS("data_clean/art.rds") %>% 
  filter(last_persontime >= 0) %>% 
  mutate(exit = case_when(is.na(exitdate) ~ 0, TRUE ~ 1),
         cohort = case_when(cohort == "RSA" ~ 0,
                            TRUE ~1))

## model / plot ##

SurvObj <- Surv(time = km$last_persontime, event = km$exit == 1)

fit2 <- survfit(SurvObj ~ cohort, data = km)

# Plotting
km_plot2 <- ggsurvplot(fit2, data = km, conf.int = TRUE, 
                       fun = "pct", 
                       palette = wes_palette("Moonrise2"),  # Use your predefined colors
                       legend = "right",
                       legend.title = "",censor = F)

kaplan_plot_tot <- km_plot2$plot + 
  coord_cartesian(xlim=c(0,362), ylim = c(95,100.1)) +
  scale_x_continuous(breaks = seq(0,360,30), expand = c(0,0)) +
  scale_y_continuous(limits = c(95,100.2), breaks = seq(95,100,1), expand = c(0,0))+
  xlab("Days after ART start")+
  theme_classic(base_size = 10) +
  scale_color_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) + 
  scale_fill_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 10,
                                                                  face = "bold"), # Set title properties here
        plot.title.position = "plot") +
  labs(title = "b | Probability of survival after ART start",
       caption = "Mean as lines, 95%-CI as ribbons")

kaplan_plot_tot

##### Join together ------------------------------------------------------------

arranged_plots <- ggpubr::ggarrange(cd4_pred1, kaplan_plot_tot, ncol = 1, common.legend = F) 
arranged_plots
ggsave(arranged_plots, filename = "results/conference.png", 
       width = 12.7, height = 12.7, units = "cm", )

