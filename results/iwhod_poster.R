library(patchwork)

#### Figure 1 #####

fig1a <- ggarrange(incidence_rate,mortality_rate, ncol = 2, common.legend = TRUE, legend = "right") %>% 
  annotate_figure(bottom = text_grob("Rate per 1,000 person years", 
                                   size=20))



fig1a
ggsave(plot = fig1a, filename = "results/mortality/fig1a.png", bg='transparent', width = 25, height = 6.5, units = "cm")

#### Figure 2 #####


fig2 <- ggarrange(aj_sup_both400, aj_rec_both350, trend_cd42, ncol = 1) %>% 
  annotate_figure(bottom = text_grob("Days after ART start",size=20))

fig2
ggsave(fig2, filename= "results/fig2.png",
       width = 24, height = 30, units = "cm")

