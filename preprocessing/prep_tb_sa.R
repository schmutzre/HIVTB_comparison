TB_sa <- read_csv("data_raw/RSA/tblTB.csv", show_col_types = FALSE) %>% 
  filter(between(reg_dmy, as.Date("2010-01-01"), as.Date("2022-12-31")))


sum(!is.na(TB_sa$tb_start_dmy))
