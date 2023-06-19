#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above. Sometimes it is necessary to click "Open in browser" afterwards.
#

library(shiny)
library(plotly)
library(ggplot2)
library(dplyr)
library(here)

#### Date
lab_long <- readRDS("../../data_clean/lab_both_long.rds")
lab_long_tb <- readRDS("../../data_clean/lab_both_long_tb.rds")

suppression_cd4 <- sqrt(350)
suppression_rna <- log10(400)
#### User Interface
#### User Interface
ui <- fluidPage(
  selectInput("id_select", "Select ID", choices = unique(lab_long$id), multiple = TRUE),
  plotlyOutput("id_plot"),
  tableOutput("id_stats"),
  selectInput("id_select_tb", "Select ID for TB", choices = unique(lab_long_tb$id), multiple = TRUE),
  plotlyOutput("id_plot_tb"),
  tableOutput("id_stats_tb")
)

#### Server
server <- function(input, output) {
  color_palette <- reactive({
    c("blue", "red", "green", "purple", "orange")
  })
  
  xrange <- range(lab_long$time_diff, na.rm=TRUE)
  
  plot_data <- function(selected_data, add_vline=FALSE){
    id_colors <- setNames(color_palette()[1:length(unique(selected_data$id))], unique(selected_data$id))
    p <- plot_ly()
    for (id in unique(selected_data$id)) {
      id_data <- selected_data[selected_data$id == id, ]
      id_color <- id_colors[as.character(id)]
      p <- p %>%
        add_trace(data = id_data, x = ~time_diff, y = ~log10(rna_value), 
                  type = 'scatter', mode = 'lines', name = paste('RNA', id),
                  line = list(color = id_color, dash = 'dot', width = 2)) %>%
        add_trace(data = id_data, x = ~time_diff, y = ~sqrt(cd4_value), 
                  type = 'scatter', mode = 'lines', name = paste('CD4', id), 
                  line = list(color = id_color, width = 2), 
                  yaxis = 'y2') 
      
      if (add_vline) {
        p <- p %>%
          add_trace(x = c(id_data$time_diff_tb[1], id_data$time_diff_tb[1]), 
                    y = c(0, max(sqrt(id_data$cd4_value), na.rm = TRUE)),
                    mode = "lines", line = list(color = "red"), inherit = FALSE,
                    name = paste('Time TB', id)) #vertical line
      }
    }
    
    p <- p %>%
      add_lines(x = c(min(selected_data$time_diff, na.rm = TRUE), max(selected_data$time_diff, na.rm = TRUE)), 
                y = c(suppression_rna, suppression_rna), 
                line = list(color = "black", dash = 'dot'), inherit = FALSE,
                name = 'Suppression RNA') %>% 
      add_lines(x = c(min(selected_data$time_diff, na.rm = TRUE), max(selected_data$time_diff, na.rm = TRUE)), 
                y = c(suppression_cd4, suppression_cd4), 
                line = list(color = "black"), yaxis = 'y2', inherit = FALSE,
                name = 'Suppression CD4') %>% 
      layout(yaxis = list(title = 'RNA (log10 scale)'), 
             yaxis2 = list(title = 'CD4 (square root scale)', overlaying = 'y', side = 'right', showgrid = FALSE))
    
    p <- p %>%
      layout(
        yaxis = list(
          title = 'RNA (log10 scale)', 
          range = c(1, log10(max(selected_data$rna_value, na.rm = TRUE)))
        ), 
        yaxis2 = list(
          title = 'CD4 (square root scale)', 
          overlaying = 'y', 
          side = 'right', 
          showgrid = FALSE
        )
      )
    p
  }
  
  output$id_plot <- renderPlotly({
    selected_data <- lab_long[lab_long$id %in% input$id_select, ]
    plot_data(selected_data)
  })
  
  output$id_plot_tb <- renderPlotly({
    selected_data_tb <- lab_long_tb[lab_long_tb$id %in% input$id_select_tb, ]
    plot_data(selected_data_tb, add_vline=TRUE)
  })
  
  output$id_stats <- renderTable({
    selected_data <- lab_long[lab_long$id %in% input$id_select, ]
    stats <- selected_data %>%
      mutate(rna_value = log10(rna_value),
             cd4_value = sqrt(cd4_value)) %>%
      group_by(id) %>%
      summarise(min_rna = min(rna_value, na.rm = TRUE),
                max_rna = max(rna_value, na.rm = TRUE),
                min_cd4 = min(cd4_value, na.rm = TRUE),
                max_cd4 = max(cd4_value, na.rm = TRUE),
                median_rna = median(rna_value, na.rm = TRUE),
                median_cd4 = median(cd4_value, na.rm = TRUE),
                n_measurements_cd4 = sum(!is.na(cd4_value)),
                n_measurements_rna = sum(!is.na(rna_value))) %>%
      mutate(id = as.integer(id)) 
    
    stats
  })
  
  output$id_stats_tb <- renderTable({
    selected_data_tb <- lab_long_tb[lab_long_tb$id %in% input$id_select_tb, ]
    stats_tb <- selected_data_tb %>%
      mutate(rna_value = log10(rna_value),
             cd4_value = sqrt(cd4_value)) %>%
      group_by(id) %>%
      summarise(min_rna = min(rna_value, na.rm = TRUE),
                max_rna = max(rna_value, na.rm = TRUE),
                min_cd4 = min(cd4_value, na.rm = TRUE),
                max_cd4 = max(cd4_value, na.rm = TRUE),
                median_rna = median(rna_value, na.rm = TRUE),
                median_cd4 = median(cd4_value, na.rm = TRUE),
                n_measurements_cd4 = sum(!is.na(cd4_value)),
                n_measurements_rna = sum(!is.na(rna_value))) %>%
      mutate(id = as.integer(id)) 
    
    stats_tb
  })
}

shinyApp(ui = ui, server = server)

