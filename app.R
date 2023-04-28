library(shiny)
library(here)
library(plotly)

cor_mat <- read.csv(here('data/gene_ephys_corrs.csv'), row.names = 'gene') 
alphabetical_indices <- order(colnames(cor_mat))
cor_mat <- as.matrix(cor_mat[, alphabetical_indices])
gene_order <- read.csv(here('data/gene_order.csv'), row.names = 1)
gene_order <- gene_order$x
efeats_order <- read.csv(here('data/efeats_order.csv'), row.names = 1)
efeats_order <- efeats_order$x



forPlot <- read.csv(here('data/forPlot.csv'))


reference_celltype_colors <- c("C-PEP1"="#F8766D", "C-PEP2"="#E88526", "C-LTMR"="#D39200", "Abeta-HTMR"="#B79F00",
                               "C-PRURI"="#93AA00", "COOL"="#5EB300", "C-PEP3"="#00BA38", "Abeta-LTMR1a"="#00BF74",
                               "Adelta"="#00C19F", "C-COLD"="#00BFC4", "PROPRIO"="#00B9E3", "Abeta-LTMR2"="#00ADFA",
                               "C-NP"="#619CFF", "C-PEP4"="#AE87FF", "Adelta-LTMR"="#DB72FB", "Abeta-LTMR1b"="#F564E3",
                               "Adelta-HTMR"="#FF61C3")


create_gene_plot <- function(gene, efeature, data) {
  correlation <- cor(log1p(data[[gene]]), data[[efeature]], method = "spearman", use="pairwise.complete.obs")
  plot_title <- paste("Correlation:", round(correlation, 4))
  
  plot <- ggplot(data, aes(x = log1p(!!sym(gene)), y = !!sym(efeature))) +
    geom_point(aes(color = labels.p), size = 2) +
    geom_smooth(method = 'lm', se = FALSE, color = "black") +
    scale_color_manual(values = reference_celltype_colors) +
    ggtitle(plot_title)
  #stat_cor(method = "spearman", label.x = min(log1p(data[[gene]])), label.y = max(data[[efeature]], size = 4)) +
  #geom_text_repel(aes(label = labels.p), size = 3, force = 10) +
  #theme(legend.position = "none")
  
  return(plot)
}

generate_hover_text <- function(data) {
  hover_text <- matrix("", nrow = nrow(data), ncol = ncol(data))

  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      hover_text[i, j] <- paste("efeat:", colnames(data)[j], "<br>",
                                "gene:", rownames(data)[i], "<br>",
                                "correlation:", round(data[i, j], 4))
    }
  }

  return(hover_text)
}

ui <- fluidPage(
  fluidRow(
    column(12, actionButton("switch_button", "Switch clustered/alphabetical heatmap")),
    column(6, plotlyOutput("heat", height = "1200px")),
    column(6, plotlyOutput("scatterplot"))
  )
)

server <- function(input, output, session) {
  # clustered cormat
  reordered_cor_mat <- as.matrix(cor_mat[gene_order, efeats_order])
  
  # Heatmap data reactive variable
  heatmap_data <- reactiveVal(reordered_cor_mat)
  
  # Observe the switch button
  observeEvent(input$switch_button, {
    if (identical(heatmap_data(), reordered_cor_mat)) {
      cor_mat_alphabetical <- as.matrix(cor_mat)
      rownames(cor_mat_alphabetical) <- rev(rownames(cor_mat_alphabetical))
      heatmap_data(cor_mat_alphabetical)
    } else {
      heatmap_data(reordered_cor_mat)
    }
  })
  
  
  output$heat <- renderPlotly({
    hover_text <- generate_hover_text(heatmap_data())
    
    plot_ly(source = "heat_plot", colors = "RdBu", reversescale = T, ) %>%
      add_heatmap(
        x = colnames(heatmap_data()), 
        y = rownames(heatmap_data()), 
        z = heatmap_data(),
        text = hover_text, # Add the custom hover text
        hoverinfo = "text" # Display the custom hover text
      ) %>% 
      colorbar(limits = c(-1, 1)) %>% 
      layout(
        margin = list(l = 180), # Increase the left margin to make more space for row names
        yaxis = list(
          tickmode = "array",
          tickvals = seq(0, length(rownames(heatmap_data())) - 1, 1),
          ticktext = rownames(heatmap_data()),
          tickfont = list(size = 8) # Set the font size for row names
        ),
        xaxis = list(
          tickfont = list(size = 8) # Set the font size for column names
        )
      )
  })
  
  output$scatterplot <- renderPlotly({
    clickData <- event_data("plotly_click", source = "heat_plot")
    if (is.null(clickData)) return(NULL)
    
    gp <- create_gene_plot(gene=clickData[["y"]], efeature=clickData[["x"]], data = forPlot)
    ggplotly()
  })
  
}

shinyApp(ui, server)