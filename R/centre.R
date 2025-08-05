#' Centre Functional Time Series
#' 
#' Centers a functional time series using various methods.
#' 
#' @param x A functional time series object.
#' @param type Type of centering: "mean", "var", "median", or "trimmed".
#' 
#' @return A centered functional time series.
#' 
#' @details
#' The available centering types are:
#' \itemize{
#'   \item "mean": Center by mean function
#'   \item "var": Center by variance function
#'   \item "median": Center by median function
#'   \item "trimmed": Center by trimmed mean function
#' }
#' 
#' @examples
#' \dontrun{
#' # Load example data
#' data(pm_10_GR)
#' 
#' # Center by mean
#' centered_data <- centre(pm_10_GR, type = "mean")
#' }
#' 
#' @seealso \code{\link{mean.fts}}, \code{\link{median.fts}}
#' 
#' @export
centre <- function(x, type)
{
    switch(type, mean = func.mean(t(x)), var = func.var(t(x)),
                 median = depth.FM_fun(t(x))$median,
                 trimmed = depth.FM_fun(t(x))$mtrim)
}
