##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  lubridate, # for date handling
  AER, # for dispersiontest
  ggplot2, # for plotting
  haven,
  survival,
  survminer,
  gridExtra,
  Epi
)

##### data import/preprocessing ----

kaplan_ch <- readRDS("data_clean/supress_df.rds") %>% 
  filter(persontime_years.suppression >= 0)  #when last follow up date was before art start

kaplan_test <- kaplan_ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA")))

kaplan_sa <- #...
  
kaplan <- kaplan_test #... #join them together 

#### Kaplan Meier model (inc. results/plotting) ----

# Create a survival object with your time and event variables
km.surv_obj <- Surv(kaplan$persontime_years.suppression, kaplan$incidence.sup)

# Use the survfit() function to calculate survival probabilities
km.fit <- survfit(km.surv_obj ~ kaplan$cohort)  

# Generate Kaplan-Meier plot
kaplan_plot <- ggsurvplot(
  km.fit, 
  data = kaplan, 
  pval = TRUE, 
  xlab = "Time after starting ART in years",
  ylab = "Survival without viral suppression",
  risk.table = TRUE,
  conf.int = TRUE,
  palette = c("#FFA07A", "#B0C4DE"),
  legend.labs = 
    c("Switzerland", "South Africa")) 

print(kaplan_plot)

# Create a zoomed in plot
zoom_plot <- kaplan_plot$plot + 
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Survival without viral suppression") +
  theme_bw()

print(zoom_plot)

ggsave(filename = "results/sup/sup-km_zoom.png", plot = zoom_plot, width = 10, height = 8, dpi = 300)

#### log-rank test ----

# Run the log-rank test
log_rank_test <- survdiff(km.surv_obj ~ kaplan$cohort)  

# Print the results
print(log_rank_test)

#p<.05 --> significant differents between the two cohorts. 

#### Kaplan Meier (switzerland only) ----

# Create a survival object with time and event variables
km.surv_objCH <- Surv(kaplan_ch$persontime_years.suppression, kaplan_ch$incidence.sup)

# Use the survfit() function to calculate survival probabilities
km.fitCH <- survfit(km.surv_objCH ~ 1)  

# Generate Kaplan-Meier plot
kaplan_plotCH <- ggsurvplot(
  km.fitCH, 
  data = kaplan_ch, 
  pval = TRUE, 
  xlab = "Time after starting ART in years",
  ylab = "Survival without viral suppression",
  risk.table = TRUE,
  conf.int = TRUE,
  palette = c("#FFA07A", "#B0C4DE"),
  legend.labs = 
    c("Switzerland")) 

# Create a zoomed in plot
zoom_plotCH <- kaplan_plotCH$plot + 
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Probability of being viral suppression free") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))

print(kaplan_plotCH)
print(zoom_plotCH)

ggsave(filename = "results/sup/sup-km_CH.png", plot = zoom_plotCH, width = 10, height = 8, dpi = 300)

#### Aalen Johansen Model ----

# Create the multi-state survival object
aj_fit <- Surv(kaplan$persontime_years.suppression, kaplan$event_type, type = "mstate")

# Fit the survival curves with cohort as a grouping variable
fit <- survfit(aj_fit ~ kaplan$cohort)

# Define the colors and line types
cols <- c("red", "blue")   # Assuming you have two cohorts, adjust if more
lty <- c(1, 2)             # Solid line for the first cohort and dashed for the second

# Plot the CIF for event 1
plotCIF(fit, event = 1, 
        ylim = c(1,0), 
        col = cols, 
        lty = lty,
        xlab = "Years",
        ylab = "Probability of being viral suppression free",
        yaxt = "n",
        main = "Aalen Johansen estimator")

# Add a new y-axis with transformed labels
axis(2, at = seq(0, 1, by = 0.1), 
     labels = formatC(1 - seq(0, 1, by = 0.1), digits = 3))

# Add legend
legend("topright", legend = levels(kaplan$cohort), col = cols, lty = lty)

AJ.both = recordPlot()

png(filename = "results/sup/sup-AJ-both.png", width = 800, height = 600)
replayPlot(AJ.both)
dev.off()

#### Aalen Johansen Model (switzerland only) ----

#prepare the KM estimates
times_km <- km.fitCH$time
surv_km <- km.fitCH$surv

#AJ 
fit.ch <- survfit(aj_fit ~ 1)

# Plot the 1 - CIF
plotCIF(fit.ch, 
        lwd = 1,
        event = 1, 
        ylim = c(1, 0), 
        xlim = c(0, 12),
        xlab = "Years",
        ylab = "Probability of being viral suppression free",
        yaxt = "n",
        main = "")

# Add a new y-axis with transformed labels
axis(2, at = seq(0, 1, by = 0.1), 
     labels = formatC(1 - seq(0, 1, by = 0.1), digits = 3))

lines(times_km, 1 - surv_km, col = "blue", lty = 2, lwd = 1) ## add KM line

legend("topright",           # Position of legend. Adjust as needed.
       legend = c("Aalen Johansen", "Kaplan Meier"),   # Labels
       lty = c(1, 2),        # Line types for AJ and KM
       col = c("black", "blue"), # Colors for AJ and KM
       lwd = 2,              # Line width
       cex = 0.8)           # Adjust size as needed

AJ.ch = recordPlot()
png(filename = "results/sup/sup-AJ-ch.png", width = 750, height = 500)
replayPlot(AJ.ch)
dev.off()

