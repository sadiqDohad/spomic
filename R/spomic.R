#' Spomic Class
#'
#' The `Spomic` S4 class is designed to store details, data, preprocessing steps,
#' and results for spatial analysis of single-cell spatial omics.  
#' This structure makes it easy to keep track of hyperparameters, 
#' input data, transformations, and output in one object.
#'
#' @slot details A `list` containing metadata or hyperparameters of various functions.
#' @slot df A `data.frame` containing the main dataset.
#' @slot pp An object of any type (`spp`) containing preprocessing steps.
#' @slot results A `list` containing analysis results, such as fitted models, metrics, or plots.
#'
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

#' Validate Required Columns in a Data Frame
#'
#' Checks whether a given data frame contains all required columns.
#' If any required column is missing, the function throws an error.
#'
#' @param df A `data.frame` to be validated.
#' @param required_columns A character vector of column names that must be present in `df`.
#'
#' @return This function is called for its side effect of throwing an error when required columns are missing.
#'         If all required columns are present, it returns `NULL` invisibly.
#'
#' @export
validate_columns <- function(df, required_columns) {
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_columns, collapse = ", ")))
  }
}

#' Reorder Coordinates Counterclockwise
#'
#' Reorders a set of 2D coordinates so that they are arranged in a
#' counterclockwise order around their centroid.
#'
#' @param coords A numeric matrix or data frame with exactly two columns:
#'        the first column for x-coordinates and the second for y-coordinates.
#'
#' @return A matrix or data frame of the same dimensions as `coords`,
#'         with rows reordered counterclockwise around the centroid.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Calculates the centroid of all given coordinates.
#'   \item Computes the polar angle of each point relative to the centroid.
#'   \item Sorts the points by angle to produce a counterclockwise ordering.
#' }
#'
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

#' Convert a Data Frame to a Point Pattern (ppp) Object
#'
#' Converts a data frame containing spatial coordinates and marks into a
#' `spatstat.geom::ppp` object. The observation window is defined by a
#' concave hull (counterclockwise ordered) around the points for a tight fit.
#'
#' @param df A `data.frame` with at least three columns:
#'   \describe{
#'     \item{x}{Numeric vector of x-coordinates.}
#'     \item{y}{Numeric vector of y-coordinates.}
#'     \item{cell_type}{Factor or character vector for point marks.}
#'   }
#'
#' @return A `ppp` object from the `spatstat.geom` package representing the
#'         spatial point pattern with the specified marks and tight observation window.
#'
#' @details
#' Steps performed:
#' \enumerate{
#'   \item Determine the concave hull of the points using \pkg{concaveman}.
#'   \item Reorder hull vertices counterclockwise for valid polygon construction.
#'   \item Create a \code{spatstat.geom::owin} window from the hull polygon.
#'   \item Construct the \code{spatstat.geom::ppp} object with the coordinates and marks.
#' }
#'
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

#' Create a Spomic Object from a Data Frame or CSV File
#'
#' Constructs a new \code{\linkS4class{Spomic}} object from a provided
#' data frame or path to a CSV file. The function validates required columns,
#' ensures proper data types, optionally removes missing values, and
#' generates a spatial point pattern (`ppp`) representation.
#'
#' @param p A `data.frame` or a character string giving the path to a CSV file.
#'          The data must contain the required columns: `"x"`, `"y"`,
#'          `"cell_type"`, and `"sample"`.
#' @param drop_na Logical, default `TRUE`. If `TRUE`, rows containing any `NA` values
#'        are removed using \code{\link[tidyr]{drop_na}}.
#'
#' @return A \code{\linkS4class{Spomic}} object with populated slots:
#' \describe{
#'   \item{details}{A list containing the sample ID and an empty hyperparameter list.}
#'   \item{df}{The cleaned and validated data frame.}
#'   \item{pp}{A `spatstat.geom::ppp` object representing the spatial point pattern.}
#'   \item{results}{An empty list for storing analysis results.}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Reads the data from a CSV file if `p` is a character string.
#'   \item Optionally drops rows with missing values.
#'   \item Validates the presence of required columns using \code{\link{validate_columns}}.
#'   \item Checks that `"x"` and `"y"` are numeric.
#'   \item Ensures only one unique `sample` value is present.
#'   \item Creates a \code{\linkS4class{Spomic}} object, including the point pattern via \code{\link{spomic_to_pp}}.
#' }
#'
#' @seealso
#' \code{\linkS4class{Spomic}},
#' \code{\link{validate_columns}},
#' \code{\link{spomic_to_pp}}
#'
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

#' Set Hyperparameters for a Spomic Object
#'
#' Updates the `details$hyperparameters` slot of a \code{\linkS4class{Spomic}} object
#' with analysis parameters, such as colocalization method, distance settings, and
#' simulation counts.
#'
#' @param spomic A \code{\linkS4class{Spomic}} object whose hyperparameters will be updated.
#' @param colocalization_type Character string specifying the colocalization analysis method.
#'        Default is `"Lcross"`.
#' @param fixed_distance Logical, default `FALSE`. If `TRUE`, uses a fixed distance instead of
#'        a range for analysis.
#' @param r Numeric or `NA`. The distance (or range of distances) for the analysis. Leave as `NA`
#'        if determined automatically.
#' @param csr_nsim Integer or `NA`. Number of simulations for complete spatial randomness (CSR)
#'        tests. Leave as `NA` if not applicable.
#'
#' @return The updated \code{\linkS4class{Spomic}} object with new hyperparameter values.
#'
#' @seealso \code{\linkS4class{Spomic}}, \code{\link{create_spomic}}
#' 
#' @export
set_spomic_hyperparameters <- function(spomic, colocalization_type = "Lcross", fixed_distance = FALSE, r = NA, csr_nsim = NA) {
  spomic@details$hyperparameters$colocalization_type <- colocalization_type
  spomic@details$hyperparameters$fixed_distance <- fixed_distance
  spomic@details$hyperparameters$r <- r
  spomic@details$hyperparameters$csr_nsim <- csr_nsim
  return(spomic)
}

#' Convert a Point Pattern (ppp) Object to a Data Frame
#'
#' Converts a \code{spatstat.geom::ppp} object into a \code{data.frame}
#' with columns for x and y coordinates, marks (cell types), and a sample identifier.
#'
#' @param pp A \code{spatstat.geom::ppp} object containing point coordinates and marks.
#' @param sample_name Character string giving the sample identifier to assign to
#'        the \code{sample} column in the output.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{x}{Numeric x-coordinates of points.}
#'   \item{y}{Numeric y-coordinates of points.}
#'   \item{cell_type}{Factor or character vector of point marks.}
#'   \item{sample}{Character string, identical for all rows, naming the sample.}
#' }
#'
#' @export
ppp_to_df <- function(pp, sample_name){
  df <- data.frame(x = pp$x,
                   y = pp$y,
                   cell_type = pp$marks,
                   sample = sample_name)
  return(df)
}
