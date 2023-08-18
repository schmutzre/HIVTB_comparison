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
  gridExtra,
  cmprsk,
  sjPlot
)

#### Data ----
custom_breaks <- c(16, 24, 34, 44, 100)
       
cox_ch <- readRDS("data_clean/supress_df.rds") %>% 
         filter(persontime_years.suppression >= 0)  #when last follow up date was before art start
       
cox_test <- cox_ch %>% 
         mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA"))) 
       
cox_sa <- #...
         
cox <- cox_test %>% 
         mutate(agegroup = cut(age_at_ART_start, breaks = custom_breaks, include.lowest = TRUE))
       #... #join them together 
       
cox$sex <- droplevels(cox$sex)
       
       
#### Model ----
       
# Create a survival object with time and event variables
cox.surv_obj <- Surv(cox$persontime_years.suppression, cox$incidence.sup)
       
#Model
model.cox <- coxph(cox.surv_obj ~ cohort + agegroup + sex + cd4_group + rna_group, data = cox)
       
#### results /plot ----
       
summary(model.cox)
       
plot_model(model.cox)
hr <- exp(coef(model.cox))
ci <- exp(confint(model.cox))
       
#plot odds ratio
plot <- plot_model(model.cox, 
                          vline.color = "red",
                          show.values = TRUE, 
                          value.offset = .3,
                          group.terms = c(1, 2, 2, 2, 3, 4, 4, 4, 5, 5, 5)) +
         theme_bw()+
         labs(title = "Hazard Ratios for Factors Associated with viral suppression") + 
         theme(plot.title = element_text(hjust = 0.5))
       
plot
       
ggsave(plot = plot, filename = "results/sup/hazard.png")
       
##### assumptions ----
       
       #A significant p-value (less than 0.05) indicates that the proportional hazards assumption has been violated
       #the smooth line in the plot should be more or less a horizontal line y = 0. 
       cox.fit <- cox.zph(model.cox)
       print(cox.fit)
       plot(cox.fit)
       