library(dplyr)

###age
age_quantiles <- function(data) {
  data %>% 
    summarize("25quantil" = quantile(age_at_ART_start, 0.25),
              "50quantil" = quantile(age_at_ART_start, 0.5),
              "75quantil" = quantile(age_at_ART_start, 0.75))
}

##baseline cd4 count
cd4_baseline <- function(data) {
  
  cd4_quantiles <- quantile(data$cd4_baseline, c(0.25, 0.5, 0.75), na.rm = TRUE)
  not_na_count <- sum(!is.na(data$cd4_baseline))
  min_cd4 <- min(data$cd4_baseline, na.rm = TRUE)
  max_cd4 <- max(data$cd4_baseline, na.rm = TRUE)
  
  list(data = data, quantiles = cd4_quantiles, not_na_count = not_na_count, min = min_cd4, max = max_cd4)
}

##baseline hiv-rna count
rna_baseline <- function(data) {
  
  rna_quantiles <- quantile(data$rna_baseline, c(0.25, 0.5, 0.75), na.rm = TRUE)
  not_na_count <- sum(!is.na(data$rna_baseline))
  min <- min(data$rna_baseline, na.rm = TRUE)
  max <- max(data$rna_baseline, na.rm = TRUE)
  list(data = data, quantiles = rna_quantiles, not_na_count = not_na_count, min = min, max = max)
}# i filter the ones out which are out of the timeframe and those with NA, this gives some more valus. For some, the most recent value may be a NA but there are other measures.

##ART regimens
regimens <- function(data) {
  data %>%
    group_by(regimen) %>%
    summarise(count = n()) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    arrange(desc(count)) %>%
    ungroup()
}

##drug combinations
drug_combinations <- function(data) {
  
  # Calculate the total number of entries
  total_count <- nrow(data)
  
  # Get top 5 treatments
  top_treatments <- data %>%
    group_by(treatment) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count)) %>%
    slice_head(n = 5)
  
  # Calculate count for "other" treatments
  other_count <- total_count - sum(top_treatments$count)
  
  # Aggregate "other" treatments
  other_treatments <- tibble(
    treatment = "other", 
    count = other_count
  )
  
  # Combine top 5 with "other" and then calculate percentages
  bind_rows(top_treatments, other_treatments) %>%
    mutate(percent = count / total_count * 100) %>%
    ungroup()
}

##plotting CD4/RNA groups
plot_group_counts <- function(data, plot_title) {
  
  base_groups <- data %>%
    group_by(cd4_group, rna_group) %>%
    summarise(n = n(), .groups = "drop")
  
  plot <- ggplot(base_groups, aes(x = cd4_group, y = rna_group, fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), vjust = 0.5, hjust = 0.5, color = "black") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(x = "Baseline CD4 Group",
         y = "Baseline RNA Group",
         fill = "Count",
         title = plot_title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  return(plot)
}


