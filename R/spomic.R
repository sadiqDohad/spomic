#' @export
setClass(
  "Spomic",
  slots = c(
    details = "list", # This is flexible for you to add hyperparameters of various functions that you want to record
    df = "data.frame",
    pp = "ANY",
    results = "list" # This is flexible for you to add different resutls that you want to record
  ),
  prototype = list(
    pp = NULL  # Initialize as NULL to avoid undefined slot
  )
)

#' @export
validate_columns <- function(df, required_columns) {
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_columns, collapse = ", ")))
  }
}

#' @export
reorder_counterclockwise <- function(coords) {
  # Calculate centroid
  centroid_x <- mean(coords[, 1])
  centroid_y <- mean(coords[, 2])

  # Calculate angle from centroid
  angles <- atan2(coords[, 2] - centroid_y, coords[, 1] - centroid_x)

  # Reorder points by angle (counterclockwise)
  ordered_coords <- coords[order(angles), ]
  return(ordered_coords)
}

#' @export
spomic_to_pp <- function(df) {
  xrange <- range(df$x)
  yrange <- range(df$y)

  concave_hull <- concaveman::concaveman(cbind(df$x, df$y), concavity = 50)
  concave_hull_ccw <- reorder_counterclockwise(concave_hull)
  tight_window <- spatstat.geom::owin(poly = list(x = concave_hull_ccw[, 1],
                                   y = concave_hull_ccw[, 2]))

  pp <- spatstat.geom::ppp(df$x,
            df$y,
            window = tight_window,
            # window = spatstat.geom::owin(xrange, yrange),
            marks = factor(df$cell_type))
  return(pp)
}

#' @export
create_spomic <- function(p, drop_na = TRUE) {
  # Check input type and read data
  if (is.character(p)) {
    if (!file.exists(p)) stop("File does not exist: ", p)
    df <- read.csv(p)
  } else if (is.data.frame(p)) {
    df <- p
  } else {
    stop("Input must be a data.frame or path to a CSV file.")
  }

  if (drop_na) {
    df <- tidyr::drop_na(df)
  }

  if (nrow(df) == 0) stop("Data frame is empty.")

  # Validate required columns
  required_columns <- c("x", "y", "cell_type", "sample")
  validate_columns(df, required_columns)

  # Ensure correct column types
  if (!is.numeric(df$x) || !is.numeric(df$y)) {
    stop("Columns 'x' and 'y' must contain numeric values.")
  }

  # Ensure `sample` column has a single unique value
  unique_samples <- unique(df$sample)
  if (length(unique_samples) > 1) {
    stop("Data frame contains multiple unique sample IDs. Ensure only one sample per Spomic object.")
  }

  # Create Spomic object
  # object <- new("Spomic",
  #               sample = as.character(unique_samples),
  #               df = df,
  #               pp = spomic_to_pp(df),
  #               # null_envelope = data.frame(),
  #               hyperparameters = list(),
  #               results = list()
  # )
  object <- new("Spomic",
                details = list(),
                df = df,
                pp = spomic_to_pp(df),
                results = list()
  )

  object@details$sample <- unique(object@df$sample)
  object@details$hyperparameters <- list()

  return(object)
}

#' @export
set_spomic_hyperparameters <- function(spomic, r, colocalization_type = "Kcross.inhom", csr_nsim = 100) {
  spomic@details$hyperparameters$r <- r
  spomic@details$hyperparameters$colocalization_type <- colocalization_type
  spomic@details$hyperparameters$csr_nsim <- csr_nsim

  return(spomic)
}
