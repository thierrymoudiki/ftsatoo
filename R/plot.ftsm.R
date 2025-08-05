#' Plot Functional Time Series Model
#' 
#' Creates diagnostic plots for a functional time series model.
#' 
#' @param x An object of class \code{ftsm} returned by \code{ftsm()}.
#' @param components Number of components to plot (default: all components).
#' @param components.start Starting component number (default: 0).
#' @param xlab1 Label for x-axis of basis function plots.
#' @param ylab1 Label for y-axis of basis function plots.
#' @param xlab2 Label for x-axis of coefficient plots.
#' @param ylab2 Label for y-axis of coefficient plots.
#' @param mean.lab Label for mean function plot.
#' @param level.lab Label for level component plot.
#' @param main.title Main title for the plot.
#' @param interaction.title Title for interaction plots.
#' @param basiscol Color for basis function lines.
#' @param coeffcol Color for coefficient lines.
#' @param outlier.col Color for outlier points.
#' @param outlier.pch Plotting character for outlier points.
#' @param outlier.cex Size of outlier points.
#' @param ... Additional arguments passed to plotting functions.
#' 
#' @return A plot showing basis functions and their corresponding coefficients.
#' 
#' @examples
#' \dontrun{
#' # Fit functional time series model
#' fit <- ftsm(pm_10_GR, order = 3)
#' 
#' # Plot the model
#' plot(fit)
#' 
#' # Plot only first 2 components
#' plot(fit, components = 2)
#' }
#' 
#' @seealso \code{\link{ftsm}}, \code{\link{forecast.ftsm}}
#' 
#' @export
plot.ftsm <- function (x, components, components.start = 0, xlab1 = x$y$xname, ylab1 = "Basis function",
xlab2 = "Time", ylab2 = "Coefficient", mean.lab = "Mean",
level.lab = "Level", main.title = "Main effects", interaction.title = "Interaction",
basiscol = 1, coeffcol = 1, outlier.col = 2, outlier.pch = 19,
outlier.cex = 0.5, ...)
{
    oldpar <- par(no.readonly = TRUE)
    mean <- is.element("mean", colnames(x$basis))
    level <- is.element("level", colnames(x$basis))
    if(!(mean | level))
    components.start = max(1, components.start)
    m <- mean + level
    order <- ncol(x$basis) - m
    if (missing(components))
    components <- order
    n <- components -  components.start + 1
    par(mfcol = c(2, n))
    if(components.start == 0){
        if (mean)
        plot(x$y$x, x$basis[, "mean"], type = "l", lty = 1, xlab = xlab1,
        ylab = mean.lab, main = main.title, col = basiscol,
        ...)
        if (level)
        plot.ts(x$coeff[, "level"], xlab = xlab2, ylab = level.lab,
        main = ifelse(mean, "", main.title), col = coeffcol,
        ...)
        if (m == 1)
        plot(0, 0, type = "n", xaxt = "n", yaxt = "n", bty = "n",
        xlab = "", ylab = "")
    }
    for (i in max(1, components.start):components) {
        yl1 <- ifelse(n > 1, paste(ylab1, i), ylab1)
        yl2 <- ifelse(n > 1, paste(ylab2, i), ylab2)
        plot(x$y$x, x$basis[, m + i], type = "l", lty = 1,
        xlab = xlab1, ylab = yl1, col = basiscol, ...)
        if (i == 1)
        title(interaction.title)
        plot.ts(x$coeff[, m + i], xlab = xlab2, ylab = yl2,
        col = coeffcol, ...)
        if (!is.null(x$wt)) {
            if (sum(x$wt < 0.1)/length(x$wt) < 0.2) {
                points(time(x$coeff)[x$wt < 0.1], x$coeff[x$wt <
                0.1, i + m], pch = outlier.pch, col = outlier.col,
                cex = outlier.cex)
            }
        }
    }
    par(oldpar)
}

