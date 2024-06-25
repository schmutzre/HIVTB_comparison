#### Data ----------------------------------------------------------------------

art_ch <- readRDS("data_clean/ch/art_ch.rds")
art_rsa <- readRDS("data_clean/rsa/art_rsa.rds")

#tb_ch <- readRDS("data_clean/ch/tb_ch.rds")
#tb_rsa <- readRDS("data_clean/rsa/tb_rsa.rds")

art_noTB_ch <- readRDS("data_clean/ch/art_noTB_ch.rds")
art_noTB_rsa <- readRDS("data_clean/rsa/art_noTB_rsa.rds")

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
tb_ch$born <- as.Date(paste0(tb_ch$born, "-01-01"), format = "%Y-%m-%d")
art_noTB_ch$born <- as.Date(paste0(art_noTB_ch$born, "-01-01"), format = "%Y-%m-%d")

#### Combine data & save -------------------------------------------------------

art <- rbind(art_ch, art_rsa)
saveRDS(art, "data_clean/art.rds")

tb <- rbind(tb_ch, tb_rsa)
saveRDS(tb, "data_clean/tb.rds")

art_noTB <- rbind(art_noTB_ch, art_noTB_rsa)
saveRDS(art_noTB, "data_clean/art_noTB.rds")
