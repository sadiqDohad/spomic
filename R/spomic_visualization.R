usethis::use_package("ggplot2")

#' @export
plot_spomic <- function(spomic, point_size = 0.5, scale_dot_size = 2) {
  spomic@df$x_rescaled <- spomic@df$x - min(spomic@df$x)
  spomic@df$y_rescaled <- spomic@df$y - min(spomic@df$y)

  g <- ggplot(spomic@df, aes(x=x_rescaled, y=y_rescaled, group=cell_type)) +
    geom_point(aes(color = cell_type), size=point_size) +
    labs(title=spomic@details$sample) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size=8),
      panel.background = element_rect(fill = 'black'),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    coord_fixed(ratio = 1)
  return(g)
}

#' @export
plot_cell_proportions <- function(spomic, stack = TRUE) {
  proportions <- spomic@df |>
    dplyr::group_by(cell_type) |>
    dplyr::summarise(count = dplyr::n()) |>
    dplyr::mutate(proportion = count / sum(count),
                  label = paste0(round(proportion * 100, 1), "%"),
                  cell_type = forcats::fct_reorder(cell_type, proportion, .desc = TRUE))

  if(stack) {
    p <- ggplot(proportions, aes(x = "", y = proportion, fill = cell_type)) +
      geom_bar(stat = "identity", width = 0.25) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(title = "Proportions of Cell Types", x = NULL, y = "Proportion") +
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank())
  }

  if(!stack) {
    p <- ggplot(proportions, aes(x = reorder(cell_type, proportion), y = proportion)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(x = "Cell Type", y = "Proportion", title = "Cell Type Proportion") +
      theme_minimal()
  }

  return(p)
}

#' @export
plot_cell_pair <- function(spomic_obj, i, j, point_size = 0.75){
  df <- spomic_obj@df |>
    dplyr::mutate(cells = dplyr::case_when(cell_type == !!i | cell_type == !!j ~ as.character(cell_type),
                                           TRUE ~ "Other"))

  unique_cells <- unique(df$cells)
  default_colors <- c("red", "blue")  # Default colors for cellA and cellB
  other_color <- "gray50"

  # Constructing color mapping dynamically
  color_values <- c(setNames(default_colors, c(i, j)), Other = other_color)

  g <- ggplot(df, aes(x = x, y = y, color = cells)) +
    geom_point(size = point_size) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    scale_color_manual(values = color_values) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size=8),
      panel.background = element_rect(fill = 'black'),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    coord_fixed(ratio = 1)

  return(g)
}


#' @export
plotKDE <- function(spomic, cell) {
  require(viridis)
  full_data <- spomic@df
  # Compute convex hull indices
  hull_indices <- chull(full_data$x, full_data$y)
  # Extract the convex hull points in order
  hull_data <- full_data[hull_indices, ]

  # Filtered dataset for density estimation
  filtered_data <- spomic@df |>
    dplyr::filter(cell_type == cell)

  # Create the plot with convex hull outline
  p <- ggplot() +
    coord_fixed(ratio = 1) +
    geom_density_2d_filled(data = filtered_data, aes(x, y, fill = as.numeric(after_stat(level)))) +
    scale_fill_viridis_c(option = "magma") +  # Smooth gradient color
    geom_polygon(data = hull_data, aes(x, y), fill = NA, color = "white", linewidth = 1) +

    # Title & Theme
    ggtitle(paste("Kernel Density of", cell)) +
    # ggpubr::theme_pubclean() +
    labs(fill = "Kernel Density") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size=8),
      panel.background = element_rect(fill = 'black'),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    coord_fixed(ratio = 1)
  p
  return(p)
}

