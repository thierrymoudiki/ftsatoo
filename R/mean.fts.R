#' Mean Function for Functional Time Series
#' 
#' Computes the mean function of a functional time series using various methods.
#' 
#' @param x A functional time series object of class \code{fts}, \code{fds}, or \code{sfts}.
#' @param method Method for computing the mean: "coordinate" (default), "FM", "mode", "RP", "RPD", or "radius".
#' @param na.rm Logical. If TRUE, remove missing values (default: TRUE).
#' @param alpha Parameter for radius depth method.
#' @param beta Parameter for radius depth method.
#' @param weight Weight parameter for radius depth method.
#' @param ... Additional arguments passed to depth functions.
#' 
#' @return A list containing:
#' \item{x}{Grid points}
#' \item{y}{Mean function values}
#' 
#' @details
#' The available methods are:
#' \itemize{
#'   \item "coordinate": Coordinate-wise mean
#'   \item "FM": Fraiman-Muniz depth-based mean
#'   \item "mode": Mode depth-based mean
#'   \item "RP": Random projection depth-based mean
#'   \item "RPD": Random projection depth-based mean
#'   \item "radius": Radius depth-based mean
#' }
#' 
#' @examples
#' \dontrun{
#' # Load example data
#' data(pm_10_GR)
#' 
#' # Compute mean function
#' mean_func <- mean(pm_10_GR, method = "coordinate")
#' 
#' # Plot mean function
#' plot(mean_func$x, mean_func$y, type = "l", 
#'      xlab = "Time", ylab = "Mean function")
#' }
#' 
#' @seealso \code{\link{median.fts}}, \code{\link{var.fts}}, \code{\link{sd.fts}}
#' 
#' @references
#' Fraiman, R., & Muniz, G. (2001). Trimmed means for functional data.
#' Test, 10(2), 419-440.
#' 
#' @export
mean.fts <- function (x, method = c("coordinate", "FM", "mode", "RP", "RPD", "radius"), 
                         na.rm = TRUE, alpha, beta, weight, ...) 
{
   if (class(x)[1] == "fts"|class(x)[1] == "fds"|class(x)[1] == "sfts"){
       method = match.arg(method)
       if (method == "coordinate"){
          loc <- rowMeans(x$y, na.rm = na.rm)
       }
       if (method == "FM"){
          loc <- depth.FM(x)$mtrim
       }
       if (method == "mode"){
          loc <- depth.mode(x)$mtrim
       }
       if (method == "RP"){
          loc <- depth.RP(x)$mtrim
       }
       if (method == "RPD"){
          loc <- depth.RPD(x)$mtrim    
       }
       if (method == "radius"){
       	  loc <- depth.radius(x, alpha, beta, weight)$mtrim
       }
       if (class(x)[1] == "fds"){
           warning("Object is not a functional time series.")
       }
       return(list(x = x$x, y = loc))
   }
   else {
        stop("Not a functional object.")
   }
}
