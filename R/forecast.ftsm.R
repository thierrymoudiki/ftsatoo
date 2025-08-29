#' Forecast Functional Time Series Model
#' 
#' Forecasts a functional time series model using various univariate time series methods.
#' 
#' @param object An object of class \code{ftsm} returned by \code{ftsm()}.
#' @param h Forecast horizon (number of periods to forecast).
#' @param method Forecasting method: "ets" (default), "arima", "ar", "ets.na", "rwdrift", "rw", "struct", or "arfima".
#' @param level Confidence level for prediction intervals (default: 80).
#' @param jumpchoice Method for handling the jump-off point: "fit" (default) or "actual".
#' @param pimethod Method for prediction intervals: "parametric" (default) or "nonparametric".
#' @param B Number of bootstrap replications for nonparametric prediction intervals (default: 100).
#' @param usedata Number of observations to use for fitting (default: all available).
#' @param adjust Logical. If TRUE, adjust forecasts for bias (default: TRUE).
#' @param model ETS model specification (optional).
#' @param damped Logical. If TRUE, use damped trend (optional).
#' @param stationary Logical. If TRUE, force stationarity (default: FALSE).
#' @param ... Additional arguments passed to the forecasting method.
#' 
#' @return An object of class \code{ftsf} containing:
#' \item{method}{Forecasting method used}
#' \item{x}{Time points}
#' \item{y}{Original functional time series}
#' \item{fitted}{Fitted values}
#' \item{residuals}{Residuals}
#' \item{mean}{Point forecasts}
#' \item{lower}{Lower prediction intervals}
#' \item{upper}{Upper prediction intervals}
#' \item{level}{Confidence level}
#' \item{xname}{Name of x variable}
#' \item{yname}{Name of y variable}
#' 
#' @examples
#' \dontrun{
#' # Fit functional time series model
#' fit <- ftsm(pm_10_GR, order = 3)
#' # Forecast 12 periods ahead
#' fc <- forecast(fit, h = 12, method = "ets")
#' # Plot forecasts
#' plot(fc)
#' 
#' fit <- ftsm(pm_10_GR, order = 3, mean=FALSE); fc <- forecast(fit, h = 12, method = "ridge2f")
#' # Forecast 12 periods ahead
#' # Plot forecasts
#' plot(fc)
#' 
#' }
#' 
#' @seealso \code{\link{ftsm}}, \code{\link{plot.ftsf}}
#' 
#' @references
#' Hyndman, R.J., & Shang, H.L. (2009). Forecasting functional time series.
#' Journal of the Korean Statistical Society, 38(3), 199-221.
#' 
#' @export
forecast.ftsm <- function (object, h = 10, method = c("ets", "arima", "ar", "ets.na",
    "rwdrift", "rw", "struct", "arfima", "ridge2f"), level = 80, jumpchoice = c("fit", "actual"), 
    pimethod = c("parametric", "nonparametric"), B = 100, usedata = nrow(object$coeff), 
    adjust = TRUE, model = NULL, damped = NULL, stationary = FALSE, 
    ...)
{
    method <- match.arg(method)
    jumpchoice <- match.arg(jumpchoice)
    pimethod <- match.arg(pimethod)
    if (jumpchoice == "actual") {
        var.col <- apply(object$coeff, 2, var)
        idx <- order(var.col)[1]
        if (var.col[idx] > 1e-08) 
            stop("No mean function fitted. So jumpchoice cannot be 'actual'")
        trueval <- object$fitted$y + object$residuals$y
        n <- ncol(trueval)
        object$basis[, idx] <- trueval[, n]
        for (i in 1:ncol(object$basis)) {
            if (i != idx) 
                object$coeff[, i] <- object$coeff[, i] - object$coeff[n, i]
        }
    }
    else if (jumpchoice != "fit") 
        stop("Unknown jump choice")
    nb <- ncol(object$basis)
    l <- nrow(object$coeff)
    
    if(method=="ar" | method=="arfima")
        stationary <- TRUE
    if(any(stationary)==TRUE & !(method == "arima" | method == "ar" | method == "arfima"))
        stop("Choose a stationary method")
  
    meanfcast <- varfcast <- matrix(NA, nrow = h, ncol = nb)
    obs <- fitted <- matrix(NA, nrow = l, ncol = nb)
    qconf <- qnorm(0.5 + level/200)
    usedata <- min(usedata, l)
    ytsp <- tsp(object$coeff)
    x <- xx <- ts(as.matrix(object$coeff[(l - usedata + 1):l, 
        ]), start = ytsp[1] + l - usedata, frequency = ytsp[3])
    xx[object$weights[(l - usedata + 1):l] < 0.1, ] <- NA
    fmodels <- list()
    if (method == "ets") {
        if (is.null(model)) 
            model <- c("ANN", rep("ZZZ", nb - 1))
        else if (length(model) == 1) 
            model <- c("ANN", rep(model, nb - 1))
        else if (length(model) == nb - 1) 
            model <- c("ANN", model)
        else stop("Length of model does not match number of coefficients")
        if (!is.null(damped)) {
            if (length(damped) == 1) 
                damped <- c(FALSE, rep(damped, nb))
            else if (length(damped) == nb - 1) 
                damped <- c(FALSE, damped)
            else stop("Length of damped does not match number of coefficients")
        }
        for (i in 1:nb) {
            if (!is.null(damped)) 
                fmodels[[i]] <- ets(x[, i], model = model[i], 
                  damped = damped[i], ...)
            else fmodels[[i]] <- ets(x[, i], model = model[i], 
                ...)
            pegelsfit <- forecast(fmodels[[i]], h = h, level = level)
            meanfcast[, i] <- pegelsfit$mean
            varfcast[, i] <- ((pegelsfit$upper[, 1] - pegelsfit$lower[, 
                1])/(2 * qconf[1]))^2
            fitted[, i] <- pegelsfit$fitted
        }
    }
    else if (method == "ets.na") {      
      if (is.null(model)) 
        model <- c("ANN", rep("AZN", nb - 1))
      else if (length(model) == 1) 
        model <- c("ANN", rep(model, nb -1))
      else if (length(model) == nb - 1) 
        model <- c("ANN", model)
      else stop("Length of model does not match number of coefficients")      
      for (i in 1:nb) {
        barima <-  pegelsna(xx[, i], model = model[i])
        fitted[,i] <- fitted(barima)
        pred <- forecast(barima,h=h,level=level)
        fmodels[[i]] <- pred
        meanfcast[,i] <- pred$mean
        varfcast[,i] <- ((pred$upper[,1]-pred$lower[,1])/(2*qconf[1]))^2
      }
    }
    else if (method == "arima") {
        if (length(stationary) == 1) 
            stationary <- c(TRUE, rep(stationary, nb))
        else if (length(stationary) == nb - 1) 
            stationary <- c(TRUE, stationary)
        else stop("Length of stationary does not match number of coefficients")
        for(i in 1:nb)
        {
            if(var(xx[,i],na.rm=TRUE) < 1e-8)
            {
                cc <- mean(xx[,i],na.rm=TRUE)
                fmodels[[i]] <- list("Constant",cc)
                meanfcast[,i] <- rep(cc,h)
                varfcast[,i] <- rep(0,h)
                fitted[,i] <- rep(cc,length(xx[,i]))
            }
            else
            {
                barima <- auto.arima(xx[,i],stationary=stationary[i],...)
                fit_barima <- fitted(barima)
                if(length(fitted(barima)) != length(xx[,i]))
                {
                    if(any(is.na(fit_barima)))
                    {
                        fitted[which(is.na(xx[,i])),i] <- rep(NA, length(which(is.na(xx[,i]))))
                        fitted[which(!is.na(xx[,i])), i] <- fit_barima[-which(is.na(fit_barima))]
                    }
                    else
                    {
                        fitted[which(is.na(xx[,i])),i] <- rep(NA, length(which(is.na(xx[,i]))))
                        fitted[which(!is.na(xx[,i])), i] <- fit_barima
                    }
                }
                else
                {
                    fitted[,i] <- fit_barima
                }
                pred <- forecast(barima,h=h,level=level)
                fmodels[[i]] <- pred
                meanfcast[,i] <- pred$mean
                varfcast[,i] <- ((pred$upper[,1]-pred$lower[,1])/(2*qconf[1]))^2
            }
        }
    }
    else if (method == "rwdrift") {
      for (i in 1:nb) {
        if(var(xx[,i],na.rm=TRUE) < 1e-8)
        {
          cc <- mean(xx[,i],na.rm=TRUE)
          fmodels[[i]] <- list("Constant",cc)
          meanfcast[,i] <- rep(cc,h)
          varfcast[,i] <- rep(0,h)
          fitted[,i] <- rep(cc,length(xx[,i]))
        }
        else
        {
        barima <-  Arima(xx[,i], order = c(0,1,0), include.drift = TRUE)
        fitted[,i] <- fitted(barima)
        pred <- forecast(barima,h=h,level=level)
        fmodels[[i]] <- pred
        meanfcast[,i] <- pred$mean
        varfcast[,i] <- ((pred$upper[,1]-pred$lower[,1])/(2*qconf[1]))^2
        }
      }
    }
    else if (method == "rw") {
      for (i in 1:nb) {
        if(var(xx[,i],na.rm=TRUE) < 1e-8)
        {
          cc <- mean(xx[,i],na.rm=TRUE)
          fmodels[[i]] <- list("Constant",cc)
          meanfcast[,i] <- rep(cc,h)
          varfcast[,i] <- rep(0,h)
          fitted[,i] <- rep(cc,length(xx[,i]))
        }
        else
        {
        barima <-  Arima(xx[,i], order = c(0,1,0), include.drift = FALSE)
        fitted[,i] <- fitted(barima)
        pred <- forecast(barima,h=h,level=level)
        fmodels[[i]] <- pred
        meanfcast[,i] <- pred$mean
        varfcast[,i] <- ((pred$upper[,1]-pred$lower[,1])/(2*qconf[1]))^2
        }
      }
    }
    else if(method=="struct")
    {
         for (i in 1:nb)
         {
            if(var(xx[,i],na.rm=TRUE) < 1e-8)
            {
                cc <- mean(xx[,i],na.rm=TRUE)
                meanfcast[,i] <- rep(cc,h)
                varfcast[,i] <- rep(0,h)
                fitted[,i] <- rep(cc,length(xx[,i]))
            }
            else
            {
                    fitStruct <- struct.forecast(xx[,i],h=h,level=level,...)
                    meanfcast[,i] <- fitStruct$mean
                    varfcast[,i] <- fitStruct$var
                    fitted[,i] <- fitStruct$fitted
            }
        }
    }
    else if(method=="arfima")
    {
        for(i in 1:nb)
        {
            if(var(xx[,i],na.rm=TRUE) < 1e-8)
            {
                cc <- mean(xx[,i],na.rm=TRUE)
                fmodels[[i]] <- list("Constant",cc)
                meanfcast[,i] <- rep(cc,h)
                varfcast[,i] <- rep(0,h)
                fitted[,i] <- rep(cc,length(xx[,i]))
            }
            else
            {
                barfima <- arfima(xx[,i],...)
                fitted[,i] <- fitted(barfima)
                pred <- forecast(barfima,h=h,level=level)
                fmodels[[i]] <- pred
                meanfcast[,i] <- pred$mean
                varfcast[,i] <- ((pred$upper[,1]-pred$lower[,1])/(2*qconf[1]))^2
            }
        }
    }
    else if(method=="ridge2f")
    {
        # Load required library
        if (!requireNamespace("ahead", quietly = TRUE)) {
            stop("Package 'ahead' is required for ridge2f method but not installed.")
        }
        
        # Use ridge2f for multivariate forecasting of all coefficients
        misc::debug_print(x)
        ridge_fit <- ahead::ridge2f(
            y = x,  # x is the coefficient matrix (time series of all components)
            h = h,
            ...     # Pass through additional arguments
        )

        # Extract forecasts and prediction intervals
        meanfcast <- ridge_fit$mean
        
        # Handle prediction intervals
        if ("upper" %in% names(ridge_fit) && "lower" %in% names(ridge_fit)) {
            # If ridge2f provides prediction intervals directly
            varfcast <- ((ridge_fit$upper - ridge_fit$lower) / (2 * qconf[1]))^2
        } else {
            # If ridge2f doesn't provide intervals, use a default approach
            # This could be enhanced with bootstrap or other methods
            varfcast[] <- var(ridge_fit$residuals, na.rm = TRUE)  # Simple fallback
        }
        
        # Fitted values (if available from ridge2f)
        if ("fitted" %in% names(ridge_fit)) {
            fitted <- ridge_fit$fitted
        } else {
            # Calculate fitted values manually if needed
            fitted <- x - ridge_fit$residuals
        }
        
        # Store model for later use
        fmodels <- ridge_fit  # Store the ridge model
        misc::debug_print(ridge_fit)
    }
    else stop("Unknown method")
    ytsp <- try(tsp(object$fitted$time), silent=TRUE)
    if (inherits(ytsp, 'try-error'))
        ytsp <- tsp(object$fitted)
    error <- ts(object$coeff - fitted, start = ytsp[1], frequency = ytsp[3])
    misc::debug_print(object$coeff)
    misc::debug_print(fitted)
    ferror <- onestepfcast <- object$y
    onestepfcast$y <- object$basis %*% t(fitted)
    onestepfcast$yname <- "One step forecasts"
    misc::debug_print(onestepfcast)
    misc::debug_print(onestepfcast$y)
    misc::debug_print(object$y$y)
    #colnames(onestepfcast$y) = colnames(object$y$y)
    #ferror$y <- object$y$y - onestepfcast$y
    #ferror$yname <- "One step errors"
    basis_obj_fore = object$basis %*% t(meanfcast)
    misc::debug_print(basis_obj_fore)
    colnames(basis_obj_fore) = 1:h
    fmean <- fts(object$y$x, basis_obj_fore, start = ytsp[2] + 
        1/ytsp[3], frequency = ytsp[3], xname = object$y$xname, yname = "Forecasts")
    colnames(fmean$y) = seq(ytsp[2]+1/ytsp[3], ytsp[2]+h/ytsp[3], by=1/ytsp[3])    
    res <- object$residuals
    res$y <- res$y^2
    vx <- rowMeans(res$y, na.rm = TRUE)
    modelvar <- object$basis^2 %*% t(varfcast)
    totalvar <- sweep(modelvar, 1, vx + object$mean.se^2, "+")
    # if (adjust & nb > 1) {
    #     #adj.factor <- rowMeans(ferror$y^2, na.rm = TRUE)/totalvar[, 
    #     #    1]
    #     totalvar <- sweep(totalvar, 1, adj.factor, "*")
    # }
    # else adj.factor <- 1
    adj.factor <- 1
    if (length(qconf) > 1) 
        stop("Multiple confidence levels not yet implemented")
    if (pimethod == "parametric") {
        # SPECIAL CASE of "ridge2f"
      misc::debug_print(pimethod)
        if (method == "ridge2f" && "upper" %in% names(fmodels) && "lower" %in% names(fmodels)) {
            # Special case for ridge2f: use native prediction intervals if available
            ridge_fit <- fmodels  # The stored ridge2f model        
            # Transform ridge2f intervals back to functional space
            flower_coeff <- ridge_fit$lower %*% t(object$basis)
            fupper_coeff <- ridge_fit$upper %*% t(object$basis)        
            flower <- fts(object$y$x, t(flower_coeff), start = ytsp[2] + 
                1/ytsp[3], frequency = ytsp[3], xname = object$y$xname, yname = "Forecast lower limit")
            fupper <- fts(object$y$x, t(fupper_coeff), start = ytsp[2] + 
                1/ytsp[3], frequency = ytsp[3], xname = object$y$xname, yname = "Forecast upper limit")
            colnames(flower$y) = colnames(fupper$y) = seq(ytsp[2]+1/ytsp[3], ytsp[2]+h/ytsp[3], by=1/ytsp[3])  
            misc::debug_print(flower)
            misc::debug_print(fupper)
            misc::debug_print(meanfcast)
            misc::debug_print(ytsp)
            # Create coefficient forecasts using ridge2f intervals directly
            coeff <- list()
            for (i in 1:nb) {
                coeff[[i]] <- structure(list(
                    mean = ts(meanfcast[, i], start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3]), 
                    lower = ts(ridge_fit$lower[, i], start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3]), 
                    upper = ts(ridge_fit$upper[, i], start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3]), 
                    level = level, x = x[, i], method = method, 
                    model = fmodels
                ), class = "forecast")
            }            
            names(coeff) <- paste("Basis", 1:nb)
            misc::debug_print(coeff)
        return(structure(list(mean = fmean, lower = flower, upper = fupper, 
            fitted = onestepfcast, error = ferror, coeff = coeff, 
            coeff.error = error, var = list(model = modelvar, 
                error = vx, mean = object$mean.se^2, total = totalvar, 
                coeff = varfcast, adj.factor = adj.factor
            ), model = object), 
            class = "ftsf"))
            names(coeff) <- paste("Basis", 1:nb)
            return(structure(list(mean = fmean, lower = flower, upper = fupper, 
            fitted = onestepfcast, error = ferror, coeff = coeff, 
            coeff.error = error, var = list(model = modelvar, 
                error = vx, mean = object$mean.se^2, total = totalvar, 
                coeff = varfcast, adj.factor = adj.factor
            ), model = object), 
            class = "ftsf"))
        } else {
        tmp <- qconf * sqrt(totalvar)
        flower <- fts(object$y$x, fmean$y - tmp, start = ytsp[2] + 
            1/ytsp[3], frequency = ytsp[3], xname = object$y$xname, yname = "Forecast lower limit")
        fupper <- fts(object$y$x, fmean$y + tmp, start = ytsp[2] + 
            1/ytsp[3], frequency = ytsp[3], xname = object$y$xname, yname = "Forecast upper limit")
        colnames(flower$y) = colnames(fupper$y) = seq(ytsp[2]+1/ytsp[3], ytsp[2]+h/ytsp[3], by=1/ytsp[3])        
        coeff <- list()
        for (i in 1:nb) {
            coeff[[i]] <- structure(list(mean = ts(meanfcast[, 
                i], start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3]), lower = ts(meanfcast[,
                i] - qconf * sqrt(varfcast[, i]), start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3]), upper = ts(meanfcast[, i] + qconf * sqrt(varfcast[, i]), start = ytsp[2] + 1/ytsp[3],
                frequency = ytsp[3]), level = level, x = x[, i], method = method, 
                model = fmodels[[i]]), class = "forecast")
        }
        names(coeff) <- paste("Basis", 1:nb)
        return(structure(list(mean = fmean, lower = flower, upper = fupper, 
            fitted = onestepfcast, error = ferror, coeff = coeff, 
            coeff.error = error, var = list(model = modelvar, 
                error = vx, mean = object$mean.se^2, total = totalvar, 
                coeff = varfcast, adj.factor = adj.factor
            ), model = object), 
            class = "ftsf"))}
    }
    else { # nonparametric
        if (method == "ridge2f")
        {
            junk = ftsmPI(object, B = B, level = level, h = h, fmethod = method)
        } else {
            junk = ftsmPI2(object, B = B, level = level, h = h, ...)
        }        
        colnames(junk$lb) = colnames(junk$ub) = seq(ytsp[2]+1/ytsp[3], ytsp[2]+h/ytsp[3], by=1/ytsp[3])
        lb = fts(object$y$x, junk$lb, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3],
            xname = object$y$xname, yname = "Forecast lower limit")
        ub = fts(object$y$x, junk$ub, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3], 
            xname = object$y$xname, yname = "Forecast upper limit")
        return(structure(list(mean = fmean, bootsamp = junk$bootsamp, 
            lower = lb, upper = ub, model = object), class = "ftsf"))
    }
}
