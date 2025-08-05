#' Functional Partial Least Squares Regression
#' 
#' Performs functional partial least squares regression for functional time series.
#' 
#' @param data A functional time series object of class \code{fts}.
#' @param order Number of components to include in the model (default: 6).
#' @param type Type of PLS algorithm: "simpls" (default) or "nipals".
#' @param unit.weights Logical. If TRUE, use unit weights (default: TRUE).
#' @param weight Logical. If TRUE, use weighted PLS (default: FALSE).
#' @param beta Weight parameter for exponential weighting (default: 0.1).
#' @param interval Logical. If TRUE, compute prediction intervals (default: FALSE).
#' @param method Method for prediction intervals: "delta" (default) or "boota".
#' @param alpha Significance level for prediction intervals (default: 0.05).
#' @param B Number of bootstrap replications (default: 100).
#' @param adjust Logical. If TRUE, adjust for bias (default: FALSE).
#' @param backh Number of steps back for validation (default: 10).
#' 
#' @return An object of class \code{fm} containing:
#' \item{x1}{Time points}
#' \item{y1}{Grid points}
#' \item{ypred}{Predicted functional time series}
#' \item{y}{Original functional time series}
#' \item{Ypred}{Predicted values}
#' \item{B}{Regression coefficients}
#' \item{P}{X loadings}
#' \item{Q}{Y loadings}
#' \item{T}{X scores}
#' \item{R}{Weights}
#' \item{fitted}{Fitted values}
#' \item{residuals}{Residuals}
#' \item{meanX}{Mean of X}
#' \item{meanY}{Mean of Y}
#' \item{call}{Function call}
#' 
#' @examples
#' \dontrun{
#' # Load example data
#' data(pm_10_GR)
#' 
#' # Fit functional PLS regression
#' fit <- fplsr(pm_10_GR, order = 3)
#' 
#' # Plot results
#' plot(fit)
#' }
#' 
#' @seealso \code{\link{forecastfplsr}}, \code{\link{plotfplsr}}
#' 
#' @references
#' Hyndman, R.J., & Shang, H.L. (2009). Forecasting functional time series.
#' Journal of the Korean Statistical Society, 38(3), 199-221.
#' 
#' @export
fplsr <- function (data, order = 6, type = c("simpls", "nipals"), unit.weights = TRUE, 
    weight = FALSE, beta = 0.1, interval = FALSE, method = c("delta", 
        "boota"), alpha = 0.05, B = 100, adjust = FALSE, backh = 10) 
{
    type = match.arg(type)
    rawdata = t(data$y)
    n = dim(rawdata)[1]
    Xtrain = rawdata[1:(n - 1), ]
    Ytrain = rawdata[2:n, ]
    Xtest = as.numeric(rawdata[n, ])
    if (interval == FALSE) {
        if (type == "simpls") {
            if (unit.weights == TRUE) 
            {
                output = unitsimpls(Xtrain, Ytrain, Xtest, order, weight = weight, beta = beta)
                fitted = t(output$T %*% t(output$Q)) + colMeans(Ytrain)
                colnames(fitted) = rownames(Xtrain)
                residuals = t(Ytrain) - fitted
                
                Ypred_mat = as.matrix(output$Ypred)
                colnames(Ypred_mat) = 1:ncol(Ypred_mat)
                
                Xtrain_mat = as.matrix(colMeans(Xtrain))
                colnames(Xtrain_mat) = 1:ncol(Xtrain_mat)
                
                Ytrain_mat = as.matrix(colMeans(Ytrain))
                colnames(Ytrain_mat) = 1:ncol(Ytrain_mat)
                out = list(x1 = as.numeric(rownames(Xtrain)), 
                           y1 = as.numeric(colnames(Xtrain)), 
			                     ypred = fts(1:dim(Xtrain)[2], t(Xtrain), xname = data$xname, yname = data$yname),
			                     y = fts(1:dim(Ytrain)[2], t(Ytrain), xname = data$xname, yname = data$yname), 
                           Ypred = fts(1:dim(Ytrain)[2], Ypred_mat, xname = data$xname, yname = data$yname), 
			                     B = output$B, P = output$P, Q = output$Q, T = output$T, R = output$R, 
			                     fitted = fts(1:dim(Xtrain)[2], fitted, xname = data$xname, yname = "Fitted values"), 
                           residuals = fts(1:dim(Xtrain)[2], residuals, xname = data$xname, yname = "Residual"), 
                           meanX = fts(1:dim(Xtrain)[2], Xtrain_mat, xname = data$xname, yname = data$yname), 
                           meanY = fts(1:dim(Ytrain)[2], Ytrain_mat , xname = data$xname, yname = data$yname), 
                           call = match.call())
                return(structure(out, class = "fm"))
            }
            else {
                output = simpls(Xtrain, Ytrain, Xtest, order, weight = weight, beta = beta)
                fitted = t(output$T %*% t(output$Q)) + colMeans(Ytrain)
                colnames(fitted) = rownames(Xtrain)
                residuals = t(Ytrain) - fitted
                
                Ypred_mat = as.matrix(output$Ypred)
                colnames(Ypred_mat) = 1:ncol(Ypred_mat)
                
                Xtrain_mat = as.matrix(colMeans(Xtrain))
                colnames(Xtrain_mat) = 1:ncol(Xtrain_mat)
                
                Ytrain_mat = as.matrix(colMeans(Ytrain))
                colnames(Ytrain_mat) = 1:ncol(Ytrain_mat)
                out = list(x1 = as.numeric(rownames(Xtrain)), 
                           y1 = as.numeric(colnames(Xtrain)), 
			                     ypred = fts(1:dim(Xtrain)[2], t(Xtrain), xname = data$xname, yname = data$yname),
			                     y = fts(1:dim(Ytrain)[2], t(Ytrain), xname = data$xname, yname = data$yname), 
                           Ypred = fts(1:dim(Ytrain)[2], Ypred_mat, xname = data$xname, yname = data$yname), 
			                     B = output$B, P = output$P, Q = output$Q, T = output$T, R = output$R, 
			                     fitted = fts(1:dim(Xtrain)[2], fitted, xname = data$xname, yname = "Fitted values"), 
                           residuals = fts(1:dim(Xtrain)[2], residuals, xname = data$xname, yname = "Residual"), 
                           meanX = fts(1:dim(Xtrain)[2], Xtrain_mat, xname = data$xname, yname = data$yname), 
                           meanY = fts(1:dim(Ytrain)[2], Ytrain_mat, xname = data$xname, yname = data$yname), 
                           call = match.call())
                return(structure(out, class = "fm"))
            }
        }
        else {
            output = nipals(Xtrain, Ytrain, Xtest, order, weight = weight, beta = beta)
            
            Ypred_mat = matrix(output$Ypred, dim(Ytrain)[2], )
            colnames(Ypred_mat) = 1:ncol(Ypred_mat)

            Xtrain_mat = as.matrix(colMeans(Xtrain))            
            colnames(Xtrain_mat) = 1:ncol(Xtrain_mat)
            
            Ytrain_mat = as.matrix(colMeans(Ytrain))
            colnames(Ytrain_mat) = 1:ncol(Ytrain_mat)
            
            fitted_mat_value = t(output$fitted.values[, , order])
            colnames(fitted_mat_value) = rownames(Xtrain)
            
            residual_mat_value = t(output$residuals[, , order])
            colnames(residual_mat_value) = rownames(Xtrain)
            
            out = list(x1 = as.numeric(rownames(Xtrain)), y1 = as.numeric(colnames(Xtrain)), 
      				         ypred = fts(1:dim(Xtrain)[2], t(Xtrain), xname = data$xname, yname = data$yname),
                       y = fts(1:dim(Ytrain)[2], t(Ytrain), xname = data$xname, yname = data$yname), 
      				         Ypred = fts(1:dim(Ytrain)[2], Ypred_mat, xname = data$xname, yname = data$yname), 
      				         P = output$P, Q = output$Q, B = output$B, T = output$T, R = output$R, 
      				         meanX = fts(1:dim(Xtrain)[2], Xtrain_mat, xname = data$xname, yname = data$yname), 
      				         meanY = fts(1:dim(Ytrain)[2], Ytrain_mat, xname = data$xname, yname = data$yname), 
      				         Yscores = output$Yscores, projection = output$projection, 
      				         fitted = fts(1:dim(Xtrain)[2], fitted_mat_value, xname = data$xname, yname = "Fitted values"), 
      				         residuals = fts(1:dim(Xtrain)[2], residual_mat_value, xname = data$xname, yname = "Residual"), 
      				         Xvar = output$Xvar, Xtotvar = output$Xtotvar, 
                       call = match.call())
            return(structure(out, class = "fm"))
        }
    }
    else {
        fplsrPI(t(Xtrain), t(Ytrain), Xtest, order, method = method, 
            alpha = alpha, B = B, weight = weight, beta = beta, 
            adjust = adjust, backh = backh)
    }
}
