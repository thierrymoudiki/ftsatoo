#' Forecast a functional time series model
#'
#' @description Produces forecasts from a functional time series model using various
#' univariate time series forecasting methods applied to the coefficients.
#'
#' @param object An object of class \code{ftsm} (functional time series model).
#' @param h Number of periods for forecasting.
#' @param method Forecasting method to use for coefficients. Options include:
#'   \code{"ets"} (exponential smoothing), \code{"arima"} (ARIMA model),
#'   \code{"ar"} (autoregressive), \code{"ets.na"} (ETS with NA handling),
#'   \code{"rwdrift"} (random walk with drift), \code{"rw"} (random walk),
#'   \code{"struct"} (structural model), \code{"arfima"} (ARFIMA model).
#'   For more options, use argument \code{FUN}. 
#' @param level Confidence level for prediction intervals.
#' @param jumpchoice Method for handling the jump choice: \code{"fit"} (use fitted values)
#'   or \code{"actual"} (use actual values).
#' @param pimethod Method for prediction intervals: \code{"parametric"} (parametric intervals)
#'   or \code{"nonparametric"} (nonparametric/bootstrap intervals).
#' @param B Number of bootstrap replications for nonparametric prediction intervals.
#' @param usedata Number of observations to use for fitting the coefficient models.
#' @param adjust Logical indicating whether to adjust the variance using one-step errors.
#' @param model Model specification for ETS method. See details.
#' @param damped Logical or vector indicating whether to use damped trend for ETS models.
#' @param stationary Logical or vector indicating whether to force stationary models
#'   for ARIMA/AR/ARFIMA methods.
#' @param FUN Custom forecasting function to use instead of built-in methods.
#' @param ... Additional arguments passed to the forecasting functions.
#'
#' @details
#' This function forecasts a functional time series by applying univariate time series
#' forecasting methods to the coefficients obtained from a functional principal
#' component analysis. The forecasts are then transformed back to functional space.
#'
#' For the \code{method = "ets"} option, the \code{model} parameter can be specified as:
#' \itemize{
#'   \item{A single character string: Applied to all coefficients except the first}
#'   \item{A vector of length \code{nb-1}: Applied to coefficients 2 through nb}
#'   \item{A vector of length \code{nb}: Applied to all coefficients}
#' }
#' The first coefficient (mean function) always uses \code{"ANN"} model unless specified otherwise.
#'
#' @return An object of class \code{ftsf} containing:
#' \item{mean}{Point forecasts as a functional time series object}
#' \item{lower}{Lower prediction bounds (parametric method)}
#' \item{upper}{Upper prediction bounds (parametric method)}
#' \item{fitted}{One-step fitted values}
#' \item{error}{One-step forecast errors}
#' \item{coeff}{Forecasts and prediction intervals for each coefficient}
#' \item{coeff.error}{Errors in coefficient estimation}
#' \item{var}{Variance components}
#' \item{model}{The original \code{ftsm} object}
#' \item{bootsamp}{Bootstrap samples (nonparametric method)}
#'
#' @seealso \code{\link{ftsm}}, \code{\link{plot.ftsf}}
#'
#' @examples
#' \dontrun{
#' # Fit functional time series model
#' fmodel <- ftsm(y = fts(data))
#' 
#' # Forecast using ETS method
#' forecast_ets <- forecast(fmodel, h = 10, method = "ets")
#' 
#' # Forecast using ARIMA method
#' forecast_arima <- forecast(fmodel, h = 10, method = "arima")
#' 
#' # Forecast with custom function
#' custom_forecast <- forecast(fmodel, h = 10, FUN = forecast::thetaf)
#' 
#' # Plot forecasts
#' plot(forecast_ets)
#' }
#'
#' @export
#' @importFrom forecast ets auto.arima Arima arfima forecast
#' @importFrom stats tsp ts qnorm
#' @importFrom ahead genericforecast
forecast.ftsm <- function (object, h = 10, method = c("ets", "arima", "ar", "ets.na",
                                                      "rwdrift", "rw", "struct", "arfima"), level = 80, jumpchoice = c("fit", "actual"), 
                           pimethod = c("parametric", "nonparametric"), B = 100, usedata = nrow(object$coeff), 
                           adjust = TRUE, model = NULL, damped = NULL, stationary = FALSE, 
                           FUN = NULL, 
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
  
  if(!is.null(FUN))
  {
    # multivariate forecasting case
    fit_obj <- try(FUN(xx, h=h, level=level, ...), 
                   silent=TRUE)
    if (!inherits(fit_obj, "try-error"))
    { # multivariate forecasting case
      fitted <- try(fit_obj$fitted, silent=TRUE)
      if (inherits(fitted, "try-error"))
        fitted <- NULL 
      pred <- fit_obj 
      fmodels <- pred 
      meanfcast <- pred$mean
      varfcast <- ((pred$upper-pred$lower)/(2*qconf[1]))^2 
    }
    else {
      # univariate forecasting case 
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
          fit_obj <- try(ahead::genericforecast(FUN=FUN, y=xx[,i], 
                                                h=h, level=level, 
                                                ...), silent=TRUE)
          if (!inherits(fit_obj, "try-error"))
          {
            fitted[,i] <- fitted(fit_obj)
            pred <- fit_obj
            fmodels[[i]] <- pred
            meanfcast[,i] <- pred$mean
            varfcast[,i] <- ((pred$upper-pred$lower)/(2*qconf[1]))^2 
          } else {
            stop(paste("ahead::genericforecast failed to fit series", i))
          }
        }
      }
    }

  } else {
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
    else stop("Unknown method")
  }
  ytsp <- tsp(object$fitted$time)
  error <- ts(object$coeff - fitted, start = ytsp[1], frequency = ytsp[3])
  ferror <- onestepfcast <- object$y
  onestepfcast$y <- object$basis %*% t(fitted)
  onestepfcast$yname <- "One step forecasts"
  colnames(onestepfcast$y) = colnames(object$y$y)
  ferror$y <- object$y$y - onestepfcast$y
  ferror$yname <- "One step errors"
  basis_obj_fore = object$basis %*% t(meanfcast)
  colnames(basis_obj_fore) = 1:h
  fmean <- fts(object$y$x, basis_obj_fore, start = ytsp[2] + 
                 1/ytsp[3], frequency = ytsp[3], xname = object$y$xname, yname = "Forecasts")
  colnames(fmean$y) = seq(ytsp[2]+1/ytsp[3], ytsp[2]+h/ytsp[3], by=1/ytsp[3])    
  res <- object$residuals
  res$y <- res$y^2
  vx <- rowMeans(res$y, na.rm = TRUE)
  modelvar <- object$basis^2 %*% t(varfcast)
  totalvar <- sweep(modelvar, 1, vx + object$mean.se^2, "+")
  if (adjust & nb > 1) {
    adj.factor <- rowMeans(ferror$y^2, na.rm = TRUE)/totalvar[, 
                                                              1]
    totalvar <- sweep(totalvar, 1, adj.factor, "*")
  }
  else adj.factor <- 1
  if (length(qconf) > 1) 
    stop("Multiple confidence levels not yet implemented")
  if (pimethod == "parametric") {
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
                                                          coeff = varfcast, adj.factor = adj.factor), model = object), 
                     class = "ftsf"))
  }
  else {
    if (is.null(FUN)) # default forecasting models 
    {
      junk = ftsmPI(object, B = B, level = level, h = h, fmethod = method) 
    } else { # custom and multivariate 
      junk = try(ftsmPI(object, B = B, level = level, h = h, 
                        fmethod = "other", 
                    FUN=FUN, ...), silent=TRUE)
      if (inherits(junk, "try-error")) # multivariate case 
      {
        junk = try(ftsmPI2(object, B = B, level = level, h = h,
                          FUN=FUN, ...), silent=TRUE)
      }
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