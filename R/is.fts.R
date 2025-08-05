#' Check if Object is a Functional Time Series
#' 
#' Tests whether an object is a functional time series.
#' 
#' @param x An object to be tested.
#' 
#' @return Logical. TRUE if the object is a functional time series, FALSE otherwise.
#' 
#' @examples
#' \dontrun{
#' # Load example data
#' data(pm_10_GR)
#' 
#' # Check if it's a functional time series
#' is.fts(pm_10_GR)
#' }
#' 
#' @seealso \code{\link{ftsm}}
#' 
#' @export
is.fts <- function(x) 
{
    inherits(x, "fts") & length(x$x) > 0
}
