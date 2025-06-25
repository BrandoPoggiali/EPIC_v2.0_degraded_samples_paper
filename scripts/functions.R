# Project: DNA degradation study
# File: functions.R
# Description: Contains reusable functions for plotting and data analysis.

#' Creates a faceted bar plot for QC metrics.
#'
#' @param qc_df_sorted A sorted data frame of QC metrics.
#' @param y The y-axis variable to plot (e.g., mean_intensity).
#' @param y_axis_title The title for the y-axis.
#' @return A ggplot object.
bar_plotter <- function(qc_df_sorted, y, y_axis_title = "") {
  fragment_lengths <- c("95 bp", "165 bp", "230 bp", "350 bp")
  frag_len_colour <- c("#F20001", "#39A1A3", "#F1D839", "#000000", "#A8A8A8")
  names(frag_len_colour) <- c(fragment_lengths, "Control")
  
  qc_df_sorted %>%
    ggplot(aes(DNA_input, {{y}})) +
    geom_bar(
      aes(fill = as.character(DNA_size), alpha = Duplicate),
      position = position_dodge2(), stat = 'identity'
    ) +
    scale_alpha_discrete(
      range = c(1, 0.6), labels = c("A", "B"), breaks = c("A", "B")
    ) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = frag_len_colour, guide = "none") +
    labs(x = "DNA input (ng)", alpha = "Replicates", y = y_axis_title) +
    facet_grid(
      ~ DNA_size, scales = "free_x", space = "free_x", drop = TRUE
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 35, vjust = 0.5, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      strip.text = element_text(size = 14)
    )
}


#' Creates a heatmap-style table plot for probe success rates or age prediction MAE.
#'
#' @param dataframe The data frame containing the data to plot.
#' @param x_variable The x-axis variable (e.g., DNA fragment size).
#' @param y_variable The y-axis variable (e.g., DNA input).
#' @param fill_variable The variable to map to the fill color (e.g., success rate).
#' @param label_variable The variable to display as text in the tiles.
#' @param x_axis_title Title for the x-axis.
#' @param y_axis_title Title for the y-axis.
#' @param fill_legend_title Title for the fill legend.
#' @param title_plot The main title for the plot.
#' @return A ggplot object.
plot_metric_table <- function(dataframe, x_variable, y_variable, fill_variable, label_variable,
                              x_axis_title, y_axis_title, fill_legend_title = "Value", title_plot = "") {
  
  graph <- ggplot(dataframe, aes(x = {{x_variable}}, y = {{y_variable}}, fill = {{fill_variable}})) +
    geom_tile(color = "black") +
    geom_text(aes(label = {{label_variable}}), color = "black", size = 4.0) +
    coord_fixed() +
    scale_fill_gradientn(
      colours = c("#FFFFFF", "#d0b5e3", "#792cae"),
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
      name = fill_legend_title
    ) +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      legend.key.height = grid::unit(1, "cm"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid = element_blank()
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(limits = rev, expand = c(0, 0)) +
    labs(x = x_axis_title, y = y_axis_title, title = title_plot)
  
  return(graph)
}