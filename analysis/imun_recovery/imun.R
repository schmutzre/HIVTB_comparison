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

# Set seed for reproducibility
set.seed(123)

#### data ---------------------------------------------------------------

rna_ch <- readRDS("data_clean/ch/rna_ch.rds") %>% 
  mutate(id = as.factor(id),
         cohort = as.factor("CH"))
cd4_ch <- readRDS("data_clean/ch/cd4_ch.rds")%>% 
  mutate(id = as.factor(id),
         cohort = as.factor("CH"))
rna_rsa <- readRDS("data_clean/rsa/rna_rsa.rds")%>% 
  mutate(id = as.factor(id),
         cohort = as.factor("RSA"))

cd4_rsa <- readRDS("data_clean/rsa/cd4_rsa.rds") %>% 
  mutate(id = as.factor(id),
         cohort = as.factor("RSA"))

cd4 <- rbind(cd4_ch %>% dplyr::select(-timepoint),
             cd4_rsa %>% dplyr::select(-timepoint)) %>% 
  filter(time_diff >= 0 & time_diff < 365) %>% 
  mutate(
    group_30days = as.integer((time_diff %/% 30) + 1),
    cd4_trans = sqrt(cd4))

cd4_npres <- cd4 %>% 
  filter(presenting_tb == 0)

cd4_pres <- cd4 %>% 
  filter(presenting_tb == 1)

rna <- rbind(rna_ch %>% dplyr::select(-timepoint),
             rna_rsa %>% dplyr::select(-timepoint)) %>% 
  filter(time_diff >= 0 & time_diff < 365) %>% 
  mutate(group_30days = as.integer((time_diff %/% 30) + 1),
         rna_trans = log10(rna + 1))

rna_npres <- rna %>% 
  filter(presenting_tb == 0)

rna_pres <- rna %>% 
  filter(presenting_tb == 1)

#### models --------------------------------------------------------------------

### cd4 ###
#+ s(id, bs = "re")

m.cd4_npres <- gam(cd4 ~ s(time_diff, k = 3, by = cohort, bs = "cr"), 
                   data = cd4_npres, 
                   method = "REML")

m.cd4_pres <- gam(cd4 ~ s(time_diff, k = 3, by = cohort, bs = "cr"), 
                   data = cd4_pres, 
                   method = "REML")

### rna ###

m.rna_npres <- gam(rna_trans ~ s(time_diff, k = 4, by = cohort, bs = "cr"), 
             data = rna_npres, 
             method = "REML")

m.rna_pres <- gam(rna_trans ~ s(time_diff, k = 4, by = cohort, bs = "cr"), 
             data = rna_pres, 
             method = "REML")

#### data exploration ----------------------------------------------------------

# cd4 #

summary_stats_cd4 <- cd4 %>% 
  group_by(group_30days) %>% 
  summarise(
    md_trans_cd4 = median(cd4_trans, na.rm = TRUE),
    Q1 = quantile(cd4_trans, 0.25, na.rm = TRUE),
    Q3 = quantile(cd4_trans, 0.75, na.rm = TRUE),
    time_diff_midpoint = mean(time_diff, na.rm = TRUE),
    Qlow = quantile(cd4_trans, 0.05, na.rm = TRUE),
    Qhigh = quantile(cd4_trans, 0.95, na.rm = TRUE),
    mean = mean(cd4_trans, na.rm = TRUE),
    se = sd(cd4_trans, na.rm = TRUE) / sqrt(n()),
    sd = sd(cd4_trans, na.rm = TRUE)
  )

iqr_cd4 <- ggplot(summary_stats_cd4, aes(x = time_diff_midpoint, y = md_trans_cd4)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Q1, ymax = Q3),
    width = 10
  ) +
  labs(x = 'Time Difference (days)', y = 'Square root of CD4') +
  theme_bw() +
  geom_hline(yintercept = sqrt(350)) 

iqr_cd4

sd_cd4 <- ggplot(summary_stats_cd4, aes(x = time_diff_midpoint, y = mean)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean - sd, 
        ymax = mean + sd),
    width = 10
  ) +
  labs(x = 'Time Difference (days)', y = 'Square root of CD4') +
  theme_bw()+   
  geom_hline(yintercept = sqrt(350)) 
  
sd_cd4

# rna #

summary_stats_rna <- rna %>% 
  group_by(group_30days) %>% 
  summarise(
    md_trans_rna = median(rna_trans, na.rm = TRUE),
    Q1 = quantile(rna_trans, 0.25, na.rm = TRUE),
    Q3 = quantile(rna_trans, 0.75, na.rm = TRUE),
    time_diff_midpoint = mean(time_diff, na.rm = TRUE),
    Qlow = quantile(rna_trans, 0.05, na.rm = TRUE),
    Qhigh = quantile(rna_trans, 0.95, na.rm = TRUE),
    Mean = mean(rna_trans, na.rm = TRUE),
    se = sd(rna_trans, na.rm = TRUE) / sqrt(n()),
    sd = sd(rna_trans, na.rm = TRUE))

iqr_rna <- ggplot(summary_stats_rna, aes(x = time_diff_midpoint, y = md_trans_rna)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Q1, ymax = Q3),
    width = 10) +
  labs(x = 'Time Difference (days)', y = 'Log-transformed RNA') +
  theme_bw() +
  geom_hline(yintercept = log10(400)) 

sd_rna <- ggplot(summary_stats_rna, aes(x = time_diff_midpoint, y = Mean)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Mean - sd, 
        ymax = Mean + sd),
    width = 10) +
  labs(x = 'Time Difference (days)', y = 'Log-transformed RNA') +
  theme_bw() +
  geom_hline(yintercept = log10(400)) 

iqr_rna

plot_se <- grid.arrange(sd_cd4, sd_rna, ncol = 2)

ggsave(plot = plot_se, file = "results/slopes/se.png")

#### model check ---------------------------------------------------------------

# cd4 #

plot.gam(m.cd4)
summary(m.cd4)
residuals.cd4 = resid(m.cd4)
hist(residuals.cd4, main = "Histogram of Residuals", xlab = "Residuals")

ggplot(rec.cd4, aes(x = time_diff, y = residuals.cd4)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw() +
  labs(x = "Time Difference", y = "Residuals")

# rna #

plot.gam(m.rna)
summary(m.rna)
residuals.rna = resid(m.rna)
hist(residuals.rna, main = "Histogram of Residuals", xlab = "Residuals")

ggplot(rec.rna, aes(x = time_diff, y = residuals.rna)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw() +
  labs(x = "Time Difference", y = "Residuals")

#### plots ---------------------------------------------------------------------

#### Confidence intervals + raw  ----

# I calculate the simultaneous CI instead of point-wise https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

### cd4 ####

## not presenting ----

# fit only #

trend_cd4_fit_np <- plot_trend(m.cd4_npres, cd4_npres) +
  labs(title = "not presenting with TB") +
  geom_hline(yintercept = 350, linetype = "dotted") +
  coord_cartesian(ylim = c(0,600),
                  xlim = c(0,360))

trend_cd4_fit_np

# fit + raw #

trend_cd4_npres <- plot_trend.raw.cd4(m.cd4_npres, cd4_npres) +
  facet_wrap(~cohort) +
  ylab("not presenting with TB")+
  coord_cartesian(ylim = c(0,2000))

trend_cd4_npres

## presenting ----

# fit only #

trend_cd4_fit_p <- plot_trend(m.cd4_pres, cd4_pres) +
  ylab("presenting with TB") +
  geom_hline(yintercept = 350, linetype = "dotted") +
  coord_cartesian(ylim = c(0,600),
                  xlim = c(0,360)) +
  theme(axis.title.x = element_blank())

trend_cd4_fit_p

# fit + raw #

trend_cd4_pres <- plot_trend.raw.cd4(m.cd4_pres, cd4_pres) +
  facet_wrap(~cohort) +
  theme(strip.text.x = element_blank(), strip.background = element_blank())+
  ylab("presenting with TB") +
  coord_cartesian(ylim = c(0,2000))

trend_cd4_pres

## combine ----

# fit only #

fit_cd4 <- ggarrange(trend_cd4_fit_p, trend_cd4_fit_np, ncol = 2)

fit_cd4

ggsave(plot = fit_cd4, filename = "results/slopes/cd4_fit.png", 
       width = 16, height = 11, units = "cm")

# fit + raw #

trend_cd4 <- ggarrange(trend_cd4_npres, 
                             trend_cd4_pres, 
                             ncol = 1) %>% 
  annotate_figure(bottom = "Days after ART start",
                  left = "cd4")
trend_cd4

ggsave(plot = trend_cd4, filename = "results/slopes/cd4.png", 
       width = 16, height = 11, units = "cm")

### rna ####

## not presenting ----

# fit only # 

trend_rna_fit_np <- plot_trend(m.rna_npres, rna_npres) +
  ylab(" not presenting with TB") +
  geom_hline(yintercept = log10(400), linetype = "dotted") +
  coord_cartesian(ylim = c(0,8),
                  xlim = c(0,360)) +
  theme(axis.title.x = element_blank())

trend_rna_fit_np

# fit + raw #

trend_rna_npres <- plot_trend.raw.rna(m.rna_npres, rna_npres) +
  facet_wrap(~cohort) +
  ylab("not presenting with TB")

trend_rna_npres

## presenting

# fit only # 

trend_rna_fit_p <- plot_trend(m.rna_pres, rna_pres) +
  ylab("presenting with TB") +
  geom_hline(yintercept = log10(400), linetype = "dotted") +
  coord_cartesian(ylim = c(0,8),
                  xlim = c(0,360)) +
  xlab("Days since ART start")

trend_rna_fit_p

# fit + raw #

trend_rna_pres <- plot_trend.raw.rna(m.rna_pres, rna_pres) +
  facet_wrap(~cohort) +
  theme(strip.text.x = element_blank(), strip.background = element_blank())+
  ylab("presenting with TB")

trend_rna_pres

## combine

# fit only # 

fit_rna <- ggarrange(trend_rna_fit_p, trend_rna_fit_np, ncol = 1)

fit_rna

ggsave(plot = fit_rna, filename = "results/slopes/rna_fit.png", 
       width = 16, height = 11, units = "cm")+
  
fit_rna

# fit + raw #

trend_rna <- ggarrange(trend_rna_npres, 
                       trend_rna_pres, 
                       ncol = 1) %>% 
  annotate_figure(bottom = "Days after ART start",
                  left = "log10(rna)")

trend_rna

ggsave(plot = trend_rna, filename = "results/slopes/rna.png", 
       width = 16, height = 11, units = "cm")

#### prediction interval ----

### cd4

## CH non-presenting ##

cd4_npres_ch <- cd4_npres %>% 
  filter(cohort == "CH")

m.cd4_npres_ch <- gam(cd4_trans ~ s(time_diff, k = 3), 
             data = cd4_npres_ch, 
             method = "REML")

pred_cd4_npres_ch <- plot_pred.cd4(m.cd4_npres_ch, cd4_npres_ch, 2) 
pred_cd4_npres_ch

## CH presenting ##

cd4_pres_ch <- cd4_pres %>% 
  filter(cohort == "CH")

m.cd4_pres_ch <- gam(cd4_trans ~ s(time_diff, k = 3), 
                      data = cd4_pres_ch, 
                      method = "REML")

pred_cd4_pres_ch <- plot_pred.cd4(m.cd4_pres_ch, cd4_pres_ch, 2) 
pred_cd4_pres_ch

## RSA non-presenting ##

cd4_npres_rsa <- cd4_npres %>% 
  filter(cohort == "RSA")

m.cd4_npres_rsa <- gam(cd4_trans ~ s(time_diff, k = 3), 
                      data = cd4_npres_rsa, 
                      method = "REML")

pred_cd4_npres_rsa <- plot_pred.cd4(m.cd4_npres_rsa, cd4_npres_rsa, 2) 
pred_cd4_npres_rsa

## RSA presenting ##

cd4_pres_rsa <- cd4_pres %>% 
  filter(cohort == "RSA")

m.cd4_pres_rsa <- gam(cd4_trans ~ s(time_diff, k = 3), 
                     data = cd4_pres_rsa, 
                     method = "REML")

pred_cd4_pres_rsa <- plot_pred.cd4(m.cd4_pres_rsa, cd4_pres_rsa, 2) 
pred_cd4_pres_rsa

## combine them ##

pred_cd4 <- ggarrange(pred_cd4_npres_rsa, 
                      pred_cd4_npres_ch, 
                      pred_cd4_pres_rsa,
                      pred_cd4_pres_ch, 
                      ncol = 2) %>% 
  annotate_figure(bottom = "Days after ART start",
                  left = expression(sqrt(CD4)))

#### Figures for the abstract --------------------------------------------------

### Immun recovery (Fig 1a) ------

#CD4

cd4_np <- pred_trend(m.cd4_npres, cd4_npres) %>% 
  mutate(tb = factor("not presenting TB"))

cd4_p <- pred_trend(m.cd4_pres, cd4_pres) %>% 
  mutate(tb = factor("presenting TB"))

cd4_pred <- rbind(cd4_np, cd4_p)

cd4_pred <- ggplot(cd4_pred, aes(x = time_diff)) +
  geom_ribbon(aes(ymin = lwrS_adj, ymax = uprS, fill = cohort), alpha = 0.2) +
  geom_line(aes(y = fit, color = cohort)) +
  theme_classic(base_size = 10) +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  coord_cartesian(ylim = c(0,700))+
  facet_wrap(vars(tb)) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)"),
       labs(title = "a. Probability of survival after ART start"))

cd4_pred

#CD4 seperate models ----------------------------------------

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

m.cd4_rsa_npres <- gam(cd4 ~ s(time_diff, k = 3, bs = "cr"), 
                     data = cd4_npres_rsa, 
                     method = "REML")

m.cd4_ch_npres <- gam(cd4 ~ s(time_diff, k = 3, bs = "cr"), 
                     data = cd4_npres_ch, 
                     method = "REML")

m.cd4_rsa_pres <- gam(cd4 ~ s(time_diff, k = 3, bs = "cr"), 
                       data = cd4_pres_rsa, 
                       method = "REML")

m.cd4_ch_pres <- gam(cd4 ~ s(time_diff, k = 3, bs = "cr"), 
                      data = cd4_pres_ch, 
                      method = "REML")

cd4_rsa_np <- pred_trend_single(m.cd4_rsa_npres, subset(cd4_npres, cohort == "RSA")) %>% 
  mutate(tb = factor("Not presenting TB"),
         cohort = factor("RSA"))

cd4_ch_np <- pred_trend_single(m.cd4_ch_npres, subset(cd4_npres, cohort == "CH")) %>% 
  mutate(tb = factor("Not presenting TB"),
         cohort = factor("CH"))

cd4_rsa_p <- pred_trend_single(m.cd4_rsa_pres, subset(cd4_pres, cohort == "RSA")) %>% 
  mutate(tb = factor("Presenting TB"),
         cohort = factor("RSA"))  

cd4_ch_p <- pred_trend_single(m.cd4_ch_pres, subset(cd4_pres, cohort == "CH")) %>% 
  mutate(tb = factor("Presenting TB"),
         cohort = factor("CH"))

cd4_pred_single <- rbind(cd4_rsa_np, cd4_ch_np, cd4_rsa_p, cd4_ch_p)

cd4_pred1 <- ggplot(cd4_pred_single, aes(x = time_diff)) +
  geom_ribbon(aes(ymin = lwrS_adj, ymax = uprS, fill = cohort), alpha = 0.2) +
  geom_line(aes(y = fit, color = cohort)) +
  theme_classic(base_size = 10) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,360, by = 60))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,700,100)) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        plot.title = element_text(size = 10, face = "bold"), # Set title properties here
        plot.title.position = "plot") +
  scale_color_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) +
  scale_fill_manual(values = wes_palette("Moonrise2"), labels = c("South Africa", "Switzerland")) +
  coord_cartesian(ylim = c(0,700), xlim = c(0,380))+
  facet_wrap(vars(tb)) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)"),
       title = "a | CD4 count showing immunological recovery after ART start for patients presenting with/without TB")

cd4_pred1

# RNA ------

rna_np <- pred_trend(m.rna_npres, rna_npres) %>% 
  mutate(tb = factor("Not presenting TB"))

rna_p <- pred_trend(m.rna_pres, rna_pres) %>% 
  mutate(tb = factor("Presenting TB"))

rna_pred <- rbind(rna_np, rna_p)

ggplot(rna_pred, aes(x = time_diff)) +
  geom_ribbon(aes(ymin = lwrS_adj, ymax = uprS, fill = cohort), alpha = 0.2) +
  geom_line(aes(y = fit, color = cohort)) +
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  coord_cartesian(ylim = c(0,7))+
  facet_wrap(vars(tb)) +
  geom_hline(yintercept = log10(400), linetype = "dotted") +
  labs(x = "Days after ART start",
       y = "CD4 count")

###### seperate plots RNa ########

rna_npres_rsa <- rna %>% 
  filter(presenting_tb == 0,
         cohort == "RSA") 

rna_npres_ch <- rna %>% 
  filter(presenting_tb == 0,
         cohort == "CH")

rna_pres_rsa <- rna %>% 
  filter(presenting_tb == 1,
         cohort == "RSA")

rna_pres_ch <- rna %>% 
  filter(presenting_tb == 1,
         cohort == "CH")

m.rna_rsa_npres <- gam(rna_trans ~ s(time_diff, k = 3, bs = "cr"), 
                       data = rna_npres_rsa, 
                       method = "REML")

m.rna_ch_npres <- gam(rna_trans ~ s(time_diff, k = 4, bs = "cr"), 
                      data = rna_npres_ch, 
                      method = "REML")

m.rna_rsa_pres <- gam(rna_trans ~ s(time_diff, k = 4, bs = "cr"), 
                      data = rna_pres_rsa, 
                      method = "REML")

m.rna_ch_pres <- gam(rna_trans ~ s(time_diff, k = 4, bs = "cr"), 
                     data = rna_pres_ch, 
                     method = "REML")

rna_rsa_np <- pred_trend_single(m.rna_rsa_npres, subset(rna_npres, cohort == "RSA")) %>% 
  mutate(tb = factor("Not presenting TB"),
         cohort = factor("RSA"))

rna_ch_np <- pred_trend_single(m.rna_ch_npres, subset(rna_npres, cohort == "CH")) %>% 
  mutate(tb = factor("Not presenting TB"),
         cohort = factor("CH"))

rna_rsa_p <- pred_trend_single(m.rna_rsa_pres, subset(rna_pres, cohort == "RSA")) %>% 
  mutate(tb = factor("Presenting TB"),
         cohort = factor("RSA"))  

rna_ch_p <- pred_trend_single(m.rna_ch_pres, subset(rna_pres, cohort == "CH")) %>% 
  mutate(tb = factor("Presenting TB"),
         cohort = factor("CH"))

rna_pred_single <- rbind(rna_rsa_np, rna_ch_np, rna_rsa_p, rna_ch_p)

ggplot(rna_pred_single, aes(x = time_diff)) +
  geom_ribbon(aes(ymin = lwrS_adj, ymax = uprS, fill = cohort), alpha = 0.2) +
  geom_line(aes(y = fit, color = cohort)) +
  theme_classic(base_size = 8) +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  geom_hline(yintercept = log10(350), linetype = "dotted") +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  facet_wrap(vars(tb)) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)"))


### Mortality RSA (1b)

### both cohorts ---------------------------------------------------------------

col_rsa<- wes_palette("Moonrise2")[1]
col_ch <- wes_palette("Moonrise2")[2]

# Read and mutate the data first

km <- readRDS("data_clean/art.rds") %>% 
  filter(last_persontime >= 0) %>% 
  mutate(exit = case_when(is.na(exitdate) ~ 0, TRUE ~ 1),
         cohort = case_when(cohort == "RSA" ~ 0,
                            TRUE ~1))

SurvObj <- Surv(time = km$last_persontime, event = km$exit == 1)
fit <- survfit(SurvObj ~ cohort + presenting_tb, data = km)

# Plotting
km_plot <- ggsurvplot(fit, data = km, conf.int = TRUE, 
                 fun = "pct", 
                 palette = c(col_rsa,col_rsa, col_ch,col_ch),  # Use your predefined colors
                 legend = "right",
                 legend.title = "")

kaplan_plot <- km_plot$plot + 
  facet_wrap(~presenting_tb, labeller = labeller(presenting_tb = c("0" = "Not Presenting TB", "1" = "Presenting TB")))+
  coord_cartesian(xlim=c(0,380)) +
  scale_x_continuous(breaks = c(0,180,360), expand = c(0,0)) +
  scale_color_manual(values = c(col_rsa,col_rsa, col_ch,col_ch), guide = "none") + 
  scale_fill_manual(values = c(col_rsa,col_rsa, col_ch,col_ch), guide = "none") +
  xlab("Days after ART start")+
  theme_classic(base_size = 8)

#### without presenting TB strata -----

fit2 <- survfit(SurvObj ~ cohort, data = km)

# Plotting
km_plot2 <- ggsurvplot(fit2, data = km, conf.int = TRUE, 
                      fun = "pct", 
                      palette = c(col_rsa,col_ch),  # Use your predefined colors
                      legend = "right",
                      legend.title = "",censor = F)

kaplan_plot_tot <- km_plot2$plot + 
  coord_cartesian(xlim=c(0,362), ylim = c(95,100.1)) +
  scale_x_continuous(breaks = seq(0,360,30), expand = c(0,0)) +
  scale_y_continuous(limits = c(95,100.2), breaks = seq(95,100,1), expand = c(0,0))+
  xlab("Days after ART start")+
  theme_classic(base_size = 10) +
  scale_color_manual(values = c(col_rsa,col_ch), labels = c("South Africa", "Switzerland")) + 
  scale_fill_manual(values = c(col_rsa,col_ch), labels = c("South Africa", "Switzerland")) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 10,
                          face = "bold"), # Set title properties here
        plot.title.position = "plot") +
  labs(title = "b | Probability of survival after ART start")

kaplan_plot_tot


## weighted 

# library(survminer)
# library(dplyr)
# 
# # Filter the data for missing values first
# rsa_km_filtered <- rsa_km %>% 
#   filter(!is.na(presenting_tb), 
#          !is.na(age_at_art_start), 
#          !is.na(gender), 
#          !is.na(cd4_baseline))
# 
# # Fit the glm model on the filtered data
# fit_glm <- glm(presenting_tb ~ age_at_art_start + gender, 
#                family = "binomial", 
#                data = rsa_km_filtered)
# 
# # Calculate weights
# weights <- 1/predict(fit_glm, type = "response")
# 
# # Create the survival object
# SurvObj <- Surv(time = rsa_km_filtered$last_persontime, event = rsa_km_filtered$exit == 1)
# 
# # Fit the Kaplan-Meier survival curve with weights
# fit_km <- survfit(SurvObj ~ presenting_tb, data = rsa_km_filtered, weights = weights)
# 
# # Now fit_km should be correctly weighted and matched to the data
# 
# ggsurvplot(fit_km, data = rsa_km_filtered , conf.int = TRUE, 
#                       fun = "pct", 
#                       palette = c("darkolivegreen3", "darkolivegreen4"),
#                       legend.title = "", 
#                       legend.labs = c("Not presenting TB", "Presenting TB"))

### Join them together

arranged_plots <- ggpubr::ggarrange(cd4_pred1, kaplan_plot_tot, ncol = 1, common.legend = F) 

arranged_plots

library(ggpubr)
library(grid)

annotated_plots <- annotate_figure(arranged_plots,
                                   bottom = text_grob("Figure: (a) CD4 count showing immunological recovery after ART start for patients presenting with/without TB. \n(b) Probability of survival after ART start.", size = 10,
                                                      vjust = 0.2, hjust = 0, x = 0, lineheight = 1))

annotated_plots <- annotate_figure(arranged_plots,
                                   bottom = grobTree(
                                     textGrob(expression(paste(bold("Figure: (a)"), " CD4 count showing immunological recovery after ART start for patients presenting with/without TB.")), 
                                              x = unit(0, "npc"), y = unit(2.5, "npc"), hjust = 0, vjust = 0, 
                                              gp = gpar(fontsize = 10, lineheight = 0.8)),
                                     textGrob(expression(paste(bold("(b)"), " Probability of survival after ART start.")), 
                                              x = unit(0, "npc"), y = unit(1.0, "npc"), hjust = 0, vjust = 0, 
                                              gp = gpar(fontsize = 10, lineheight = 0.8))
                                   ))


annotated_plots

ggsave(arranged_plots, filename = "results/conference.png", 
       width = 12.7, height = 12.7, units = "cm", )

                  
