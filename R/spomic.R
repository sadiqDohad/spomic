#' @export
setClass(
  "Spomic",
  slots = c(
    sample = "character",
    df = "data.frame",
    pp = "ANY",  # Allow any object, including ppp
    # null_envelope = "data.frame",
    hyperparameters = "list",
    results = "list"
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
spomic_to_pp <- function(df) {
  xrange <- range(df$x)
  yrange <- range(df$y)
  pp <- spatstat.geom::ppp(df$x,
            df$y,
            window = spatstat.geom::owin(xrange, yrange),
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
  object <- new("Spomic",
                sample = as.character(unique_samples),
                df = df,
                pp = spomic_to_pp(df),
                # null_envelope = data.frame(),
                hyperparameters = list(),
                results = list()
  )

  return(object)
}

