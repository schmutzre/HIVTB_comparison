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

kaplan_ch <- readRDS("data_clean/rebound.rds") 

#AJ 

# Create the multi-state survival object
aj_fit <- Surv(kaplan_ch$supp_until_rebound, kaplan$indicator_rebound, type = "mstate")

fit.ch <- survfit(aj_fit ~ 1)

# Plot the 1 - CIF
plotCIF(fit.ch, 
        lwd = 1,
        event = 1, 
        ylim = c(1, 0), 
        xlim = c(0, 12),
        xlab = "Years",
        ylab = "Probability of not rebounding after being virally suppressed",
        yaxt = "n",
        main = "")

# Add a new y-axis with transformed labels
axis(2, at = seq(0, 1, by = 0.1), 
     labels = formatC(1 - seq(0, 1, by = 0.1), digits = 3))



AJ.ch = recordPlot()
png(filename = "results/sup/rebound_ch.png", width = 750, height = 500)
replayPlot(AJ.ch)
dev.off()

