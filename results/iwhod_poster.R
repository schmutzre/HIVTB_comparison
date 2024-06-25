library(patchwork)
library(tidyverse)

fig1a <- ggarrange(incidence_rate,mortality_rate, ncol = 2) %>% 
  annotate_figure(bottom = text_grob("Rate per 1,000 person years", 
                                   size=28))

fig1a
ggsave(plot = fig1a, filename = "results/mortality/fig1a.png", bg='transparent', width = 25, height = 12, units = "cm")

#### Figure 2 #####

dummy_data <- data.frame(x = 1, y = 1, 
                         linetype = c("Presenting with TB", "Not presenting with TB"),
                         color = c("South Africa", "Switzerland"))
# Create a dummy plot just for the legend
legend_plot <- ggplot(dummy_data, aes(x = x, y = y, linetype = linetype, color = color)) +
  geom_line(show.legend = TRUE, size = 1.1)+
  scale_linetype_manual(values = c("Presenting with TB" = "dashed", "Not presenting with TB" = "solid")) +
  scale_color_manual(values = c("South Africa" = wes_palette("Moonrise2")[1], "Switzerland" = wes_palette("Moonrise2")[2])) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.spacing.x = unit(0.5, "lines"),
        legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = T))

# Extract just the legend
legend_only <- get_legend(legend_plot)
legend_plot
fig2 <- ggarrange(aj_sup_both400, trend_cd42,legend_only, ncol = 1,
                  heights = c(1, 1,  0.2)) 

fig2
ggsave(fig2, filename= "results/fig2.png",
       width = 25, height = 26, units = "cm")

