ftsmPI <- function (object, B, level, h, fmethod = c("ets", "arima", "other"), 
                    FUN=NULL, ...) 
{
  data = object$y$y
  p = nrow(data)
  n = ncol(data)
  ncomp = dim(object$basis)[2] - 1
  if((n - ncomp - h + 1) <= 0)
  {
      ncomp = 1
  }
  mdata = apply(data, 1, mean)
  mdata2 = array(rep(as.matrix(mdata), B * h), dim = c(p, B, h))
  sdata = scale(t(data), scale = FALSE)
  load = as.matrix(svd(sdata)$v[, 1:ncomp])
  sco = sdata %*% load
  olivia = matrix(, ncomp, h)
  
  if (fmethod == "other")
  {
    stopifnot(!is.null(FUN))
    fit_obj <- lapply(ahead::genericforecast(FUN = FUN, 
                                             y = sco[, i]), 
                      h = h, ...)
    for (i in 1:ncomp) {
      olivia[i, ] = fit_obj[[i]]$mean
    }
    forerr = matrix(NA, (n - ncomp - h + 1), ncomp)
    for (i in h:(n - ncomp)) {
      k = i + (ncomp - h)
      fore = matrix(, 1, ncomp)
      for (j in 1:ncomp) {
        fore[, j] = ahead::genericforecast(FUN = FUN, 
                                           y = sco[1:k, j], 
                                           h = h, ...)$mean[h]
      }
      forerr[i - h + 1, ] = sco[k + h, ] - fore
    }
  }
  
  if (fmethod == "ets") {
    for (i in 1:ncomp) {
      olivia[i, ] = forecast(ets(sco[, i]), h = h)$mean
    }
  }
  
  if (fmethod == "arima") {
    for (i in 1:ncomp) {
      olivia[i, ] = forecast(auto.arima(sco[, i]), h = h)$mean
    }
  }
  
  if (fmethod %in% c("ets", "arima")){
    forerr = matrix(NA, (n - ncomp - h + 1), ncomp)
    for (i in h:(n - ncomp)) {
      k = i + (ncomp - h)
      fore = matrix(, 1, ncomp)
      if (fmethod == "ets") {
        for (j in 1:ncomp) {
          fore[, j] = forecast(ets(sco[1:k, j]), h = h)$mean[h]
        }
      }
      if (fmethod == "arima") {
        for (j in 1:ncomp) {
          fore[, j] = forecast(auto.arima(sco[1:k, j]), 
                               h = h)$mean[h]
        }
      }
      forerr[i - h + 1, ] = sco[k + h, ] - fore
    }
  }
  
  resi = t(sdata) - load %*% t(sco)
  q = array(NA, dim = c(p, B, h))
  for (j in 1:h) {
    for (i in 1:p) {
      q[i, , j] = sample(resi[i, ], size = B, replace = TRUE)
    }
  }
  ny = array(NA, dim = c(ncomp, B, h))
  for (j in 1:h) {
    for (i in 1:ncomp) {
      ny[i, , j] = sample(forerr[, i], size = B, replace = TRUE)
    }
  }
  oli = array(rep(olivia, B * h), dim = c(ncomp, B, h))
  fo = array(NA, dim = c(ncomp, B, h))
  for (j in 1:h) {
    for (i in 1:B) {
      fo[, i, j] = oli[, i, j] + ny[, i, j]
    }
  }
  pred = array(NA, dim = c(p, B, h))
  for (j in 1:h) {
    for (i in 1:B) {
      pred[, i, j] = load %*% fo[, i, j] + mdata2[, i, 
                                                  j] + q[, i, j]
    }
  }
  k1 = k2 = matrix(NA, p, h)
  for (j in 1:h) {
    for (i in 1:p) {
      k1[i, j] = quantile(pred[i, , j], (100 - level)/200, 
                          na.rm = TRUE)
      k2[i, j] = quantile(pred[i, , j], 1 - (100 - level)/200, 
                          na.rm = TRUE)
    }
  }
  return(list(bootsamp = pred, lb = k1, ub = k2))
}

# for ridge2
ftsmPI2 <- function(object, B, level, h, ...) {
  # Load required library
  if (!requireNamespace("ahead", quietly = TRUE)) {
    stop("Package 'ahead' is required but not installed.")
  }
  
  data = object$y$y
  p = nrow(data)
  n = ncol(data)
  ncomp = dim(object$basis)[2] - 1
  
  if((n - ncomp - h + 1) <= 0) {
    ncomp = 1
  }
  
  mdata = apply(data, 1, mean)
  mdata2 = array(rep(as.matrix(mdata), B * h), dim = c(p, B, h))
  sdata = scale(t(data), scale = FALSE)
  load = as.matrix(svd(sdata)$v[, 1:ncomp])
  sco = sdata %*% load
  
  # Generate forecasts using ridge2f for all components at once
  # ridge2f can handle multivariate series directly
  ridge_fit = ahead::ridge2f(
    y = sco,  # Pass all components as multivariate series
    h = h,
    ...  # Pass through any additional arguments
  )
  
  olivia = t(ridge_fit$mean)  # Transpose to get ncomp x h matrix
  
  # Calculate forecast errors for bootstrap
  forerr = matrix(NA, (n - ncomp - h + 1), ncomp)
  
  for (i in h:(n - ncomp)) {
    k = i + (ncomp - h)
    
    # Use ridge2f for out-of-sample forecasting on all components
    ridge_fit = ahead::ridge2f(
      y = sco[1:k, ],  # Pass subset of all components
      h = h,
      ...  # Pass through any additional arguments
    )
    
    # Extract the h-step ahead forecast for all components
    fore = ridge_fit$mean[h, ]  # h-th row contains h-step forecasts
    
    forerr[i - h + 1, ] = sco[k + h, ] - fore
  }
  
  # Calculate residuals
  resi = t(sdata) - load %*% t(sco)
  
  # Bootstrap sampling for residuals
  q = array(NA, dim = c(p, B, h))
  for (j in 1:h) {
    for (i in 1:p) {
      q[i, , j] = sample(resi[i, ], size = B, replace = TRUE)
    }
  }
  
  # Bootstrap sampling for forecast errors
  ny = array(NA, dim = c(ncomp, B, h))
  for (j in 1:h) {
    for (i in 1:ncomp) {
      ny[i, , j] = sample(forerr[, i], size = B, replace = TRUE)
    }
  }
  
  # Create bootstrap forecasts
  oli = array(rep(olivia, B * h), dim = c(ncomp, B, h))
  fo = array(NA, dim = c(ncomp, B, h))
  
  for (j in 1:h) {
    for (i in 1:B) {
      fo[, i, j] = oli[, i, j] + ny[, i, j]
    }
  }
  
  # Generate final predictions
  pred = array(NA, dim = c(p, B, h))
  for (j in 1:h) {
    for (i in 1:B) {
      pred[, i, j] = load %*% fo[, i, j] + mdata2[, i, j] + q[, i, j]
    }
  }
  
  # Calculate prediction intervals
  k1 = k2 = matrix(NA, p, h)
  for (j in 1:h) {
    for (i in 1:p) {
      k1[i, j] = quantile(pred[i, , j], (100 - level)/200, na.rm = TRUE)
      k2[i, j] = quantile(pred[i, , j], 1 - (100 - level)/200, na.rm = TRUE)
    }
  }
  
  return(list(bootsamp = pred, lb = k1, ub = k2))
}