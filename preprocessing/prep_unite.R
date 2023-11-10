#### Data ----------------------------------------------------------------------

art_ch <- readRDS("data_clean/ch/art_ch.rds")
art_rsa <- readRDS("data_clean/rsa/art_rsa.rds")

tb_ch <- readRDS("data_clean/ch/tb_ch.rds")
tb_rsa <- readRDS("data_clean/rsa/tb_rsa.rds")

## Re-order columns in RSA data ##

art_rsa <- art_rsa[ , names(art_ch)]
tb_rsa <- tb_rsa[ , names(tb_ch)]

## Compare types ##
#' This will all be solved more elegantly in the future

compare_types <- data.frame(
  Column = names(art_ch),
  Type_in_ch = sapply(art_ch, class),
  Type_in_rsa = sapply(art_rsa, class),
  stringsAsFactors = FALSE
)

# Display types that don't match
compare_types[compare_types$Type_in_ch != compare_types$Type_in_rsa, ]

art_ch$born <- as.Date(paste0(art_ch$born, "-01-01"), format = "%Y-%m-%d")
art_rsa$who_stage <- as.factor(art_rsa$who_stage)

tb_ch$born <- as.Date(paste0(tb_ch$born, "-01-01"), format = "%Y-%m-%d")
tb_rsa$who_stage <- as.factor(tb_rsa$who_stage)

#### Combine data & save -------------------------------------------------------

art <- rbind(art_ch, art_rsa)
saveRDS(art, "data_clean/art.rds")

tb <- rbind(tb_ch, tb_rsa)
saveRDS(tb, "data_clean/tb.rds")

#### For Lukas to check

library(janitor)
library(flextable)

flextable(
  tabyl(art_ch$treatment, ))

freq_table <- tabyl(art_ch$treatment, show_na = FALSE)
freq_table2 <- tabyl(art_rsa$treatment, show_na = FALSE)
# Sort by proportion and get top 5
top_ch <- freq_table[order(-freq_table$percent), ][1:5, ]
top_rsa <- freq_table2[order(-freq_table2$percent), ][1:5, ]


table(art_ch$regimen_tb)
table(art_ch$regimen_tb_group)

table(art_rsa$regimen_tb)
