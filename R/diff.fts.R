#' Difference Functional Time Series
#' 
#' Computes differences of a functional time series.
#' 
#' @param x A functional time series object of class \code{fts} or \code{sfts}.
#' @param lag Lag for differencing (default: 1).
#' @param differences Order of differencing (default: 1).
#' @param ... Additional arguments passed to \code{diff()}.
#' 
#' @return A functional time series object with differenced values.
#' 
#' @examples
#' \dontrun{
#' # Load example data
#' data(pm_10_GR)
#' 
#' # Compute first differences
#' diff_data <- diff(pm_10_GR, lag = 1, differences = 1)
#' 
#' # Plot differenced data
#' plot(diff_data)
#' }
#' 
#' @seealso \code{\link{ftsm}}, \code{\link{is.fts}}
#' 
#' @export
diff.fts <- function (x, lag = 1, differences = 1, ...) 
{
    if (class(x)[1] == "fts"|class(x)[1] == "sfts"){
        x$y <- t(diff(t(x$y), lag, differences, ...))
        return(x)
    }
    else{
        stop("object is not a functional time series or a sliced functional time series")
    }
}

