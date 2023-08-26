##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  lubridate, # for date handling
  rlang, # for '!!' inside the functions
  AER, # for dispersiontest
  ggplot2, # for plotting
  haven,
  survival,
  survminer,
  gridExtra,
  Epi,
  cmprsk
)

##### data import/preprocessing ----

kaplan_ch <- readRDS("data_clean/art_ch.rds") %>% 
  mutate(persontime_years = case_when(
           case_incident_2m == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")/360),
           case_incident_2m == 0 ~ last_persontime/360
         )
  ) %>% 
  filter(persontime_years > 0) %>% 
  dplyr::select(id, art_start_date, case_incident_2m, date_tb, cohort, persontime_years, exitdate, exit_why, last_fup_date) %>% 
  mutate(event_type = case_when(
    case_incident_2m == 1 ~ 1,
    !is.na(exitdate) ~ 2,
    TRUE ~0 # Loss to follow-up isnt considered a competing risk
  ))

kaplan_test <- kaplan_ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA")))

kaplan_sa <- #...
  
kaplan <- kaplan_test #... #join them together 

#### Kaplan Meier model (inc. results/plotting) ----

# Create a survival object with your time and event variables
km.surv_obj <- Surv(kaplan$persontime_years, kaplan$case_incident_2m)

# Use the survfit() function to calculate survival probabilities
km.fit <- survfit(km.surv_obj ~ kaplan$cohort)  

# Generate Kaplan-Meier plot
kaplan_plot <- ggsurvplot(
  km.fit, 
  data = kaplan, 
  pval = TRUE, 
  xlab = "Time after starting ART in years",
  ylab = "TB-free survival",
  risk.table = TRUE,
  conf.int = TRUE,
  palette = c("#FFA07A", "#B0C4DE"),
  legend.labs = 
    c("Switzerland", "South Africa")) 

print(kaplan_plot)

# Create a zoomed in plot
zoom_plot <- kaplan_plot$plot + 
  coord_cartesian(ylim = c(0.95, 1)) +
  labs(y = "TB-free survival") +
  theme_bw()

print(zoom_plot)

ggsave(filename = "results/km_zoom.png", plot = zoom_plot, width = 10, height = 8, dpi = 300)

#### log-rank test ----

# Run the log-rank test
log_rank_test <- survdiff(km.surv_obj ~ kaplan$cohort)  

# Print the results
print(log_rank_test)

#p<.05 --> significant differents between the two cohorts. 

#### Kaplan Meier (switzerland only) ----

# Create a survival object with your time and event variables
km.surv_objCH <- Surv(kaplan_ch$persontime_years, kaplan_ch$case_incident_2m)

# Use the survfit() function to calculate survival probabilities
km.fitCH <- survfit(km.surv_objCH ~ 1)  

# Generate Kaplan-Meier plot
kaplan_plotCH <- ggsurvplot(
  km.fitCH, 
  data = kaplan_ch, 
  pval = TRUE, 
  xlab = "Time after starting ART in years",
  ylab = "TB-free survival",
  risk.table = TRUE,
  conf.int = TRUE,
  palette = c("#FFA07A", "#B0C4DE"),
  legend.labs = 
    c("Switzerland")) 

# Create a zoomed in plot
zoom_plotCH <- kaplan_plotCH$plot + 
  coord_cartesian(ylim = c(0.995, 1)) +
  labs(y = "TB-free survival") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))

print(kaplan_plotCH)
print(zoom_plotCH)

ggsave(filename = "results/km_CH.png", plot = zoom_plotCH, width = 10, height = 8, dpi = 300)

#### Aalen Johansen Model ----
library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)


# Create the multi-state survival object
aj_fit <- Surv(kaplan$persontime_years, kaplan$event_type, type = "mstate")

# Fit the survival curves with cohort as a grouping variable
fit <- survfit(aj_fit ~ kaplan$cohort)

# Define the colors and line types
cols <- c("red", "blue")   # Assuming you have two cohorts, adjust if more
lty <- c(1, 2)             # Solid line for the first cohort and dashed for the second

# Plot the CIF for event 1
plotCIF(fit, event = 1, 
        ylim = c(0.01,0), 
        col = cols, 
        lty = lty,
        xlab = "Years",
        ylab = "TB-free survival",
        yaxt = "n",
        main = "Aalen Johansen estimator")

# Add a new y-axis with transformed labels
axis(2, at = seq(0, 0.01, by = 0.005), 
     labels = formatC(1 - seq(0, 0.01, by = 0.005), digits = 3))

# Add legend
legend("topright", legend = levels(kaplan$cohort), col = cols, lty = lty)

AJ.both = recordPlot()

png(filename = "results/incidence/AJ.both.png", width = 800, height = 600)
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
        ylim = c(0.005, 0), 
        xlim = c(0, 12),
        xlab = "Years",
        ylab = "TB-free survival",
        yaxt = "n",
        main = "Aalen Johansen estimator vs. Kaplan Meier estimator")

# Add a new y-axis with transformed labels
axis(2, at = seq(0, 0.005, by = 0.001), 
     labels = formatC(1 - seq(0, 0.005, by = 0.001), digits = 3))

lines(times_km, 1 - surv_km, col = "blue", lty = 2, lwd = 1) ## add KM line

legend("topright",           # Position of legend. Adjust as needed.
       legend = c("Aalen Johansen", "Kaplan Meier"),   # Labels
       lty = c(1, 2),        # Line types for AJ and KM
       col = c("black", "blue"), # Colors for AJ and KM
       lwd = 2,              # Line width
       cex = 0.8)           # Adjust size as needed

AJ.ch = recordPlot()
png(filename = "results/incidence/AJ.ch.png", width = 750, height = 500)
replayPlot(AJ.ch)
dev.off()

## Gray's test

gray.result <- cuminc(ftime = kaplan$persontime_years, 
                      fstatus = kaplan$event_type, 
                      group = kaplan$cohort,
                      cencode = 0)

## Alternative Plot for Aalen Johansen with both states included (could also be combinded with both cohorts)

install.packages("mstate")
library(mstate)

# Compute CIF
cif <- Cuminc(time = kaplan$persontime_years, 
              status = kaplan$event_type) 

theme_set(theme_bw(base_size = 8))

Aalen2states <- plot(cif, use.ggplot = TRUE,
     conf.int = 0.95)



