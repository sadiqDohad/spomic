usethis::use_package("ggplot2")


#' @export
plotCellProportions <- function(spomic, stack = TRUE) {
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


#' #' @export
#' plotSpomic <- function(spomic_obj, point_size = 0.75){
#'   g <- ggplot(spomic_obj@df, aes(x = x, y = y, color = cell_type)) +
#'     geom_point(size = point_size) +
#'     coord_fixed(ratio = 1) +
#'     theme_minimal() +
#'     theme(plot.background = element_rect(fill = "black"),
#'           panel.grid.major = element_blank(),
#'           panel.grid.minor = element_blank(),
#'           axis.text = element_blank(),
#'           legend.title = element_blank(),
#'           legend.position = "bottom",
#'           legend.text = element_text(color = "white"),
#'           legend.spacing.y = unit(0.1, 'cm'),    # Reduce vertical spacing
#'           legend.key.size = unit(0.5, 'cm'),       # Increase color legend size
#'           legend.background = element_rect(fill = "black", color = NA),  # Black legend background
#'           legend.box.background = element_rect(fill = "black", color = NA)  # Black legend box background
#'     ) +
#'     guides(color = guide_legend(override.aes = list(size = 2.5)))  # Control color legend point size
#'
#'   return(g)
#' }
#' @export
plot_spomic <- function(spomic, point_size = 0.75){
  g <- ggplot(spomic@df, aes(x = x, y = y, color = cell_type)) +
    geom_point(size = point_size) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(color = "white"),
          legend.spacing.y = unit(0.1, 'cm'),    # Reduce vertical spacing
          legend.key.size = unit(0.5, 'cm'),       # Increase color legend size
          legend.background = element_rect(fill = "black", color = NA),  # Black legend background
          legend.box.background = element_rect(fill = "black", color = NA)  # Black legend box background
    ) +
    guides(color = guide_legend(override.aes = list(size = 2.5)))  # Control color legend point size

  return(g)
}

#' @export
plotCellPair <- function(spomic_obj, cellA, cellB, point_size = 0.75){
  df <- spomic_obj@df |>
    dplyr::mutate(cells = dplyr::case_when(cell_type == !!cellA | cell_type == !!cellB ~ as.character(cell_type),
                                           TRUE ~ "Other"))

  unique_cells <- unique(df$cells)
  default_colors <- c("red", "blue")  # Default colors for cellA and cellB
  other_color <- "gray50"

  # Constructing color mapping dynamically
  color_values <- c(setNames(default_colors, c(cellA, cellB)), Other = other_color)

  g <- ggplot(df, aes(x = x, y = y, color = cells)) +
    geom_point(size = point_size) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    scale_color_manual(values = color_values) +
    theme(plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(color = "white"),
          legend.spacing.y = unit(0.1, 'cm'),    # Reduce vertical spacing
          legend.key.size = unit(0.5, 'cm'),       # Increase color legend size
          legend.background = element_rect(fill = "black", color = NA),  # Black legend background
          legend.box.background = element_rect(fill = "black", color = NA)  # Black legend box background
    ) +
    guides(color = guide_legend(override.aes = list(size = 2.5)))
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
    ggtitle(paste(cell, "Kernel Density with Convex Hull Outline")) +
    ggpubr::theme_pubclean() +
    labs(fill = "Kernel Density")

  return(p)
}
