#' Functional Autocorrelation Function
#' 
#' Computes the functional autocorrelation function for a functional time series.
#' 
#' @param fun_data A matrix or data frame containing functional time series data.
#' @param lag_value_range Vector of lag values to compute autocorrelation for (default: seq(0, 20, by = 1)).
#' 
#' @return A vector of autocorrelation values corresponding to the lag values.
#' 
#' @details
#' The functional autocorrelation function measures the correlation between 
#' functional observations at different time lags. For lag 0, the value is set to NA.
#' 
#' @examples
#' \dontrun{
#' # Load example data
#' data(pm_10_GR)
#' 
#' # Compute functional autocorrelation
#' acf_values <- facf(pm_10_GR$y)
#' 
#' # Plot autocorrelation function
#' plot(seq(0, 20, by = 1), acf_values, type = "l", 
#'      xlab = "Lag", ylab = "Autocorrelation")
#' }
#' 
#' @seealso \code{\link{ftsm}}, \code{\link{forecast.ftsm}}
#' 
#' @export
facf <-
function(fun_data, lag_value_range = seq(0, 20, by = 1))
{
    center_dat = scale(fun_data, center = TRUE, scale = FALSE)
    T = nrow(center_dat)
    gamma_l <- function(lag, T)
    {
        gamma_lag_sum = 0
        if(lag >= 0)
        {
            for(ij in 1:(T-lag))
            {
                gamma_lag_sum = gamma_lag_sum + as.matrix(center_dat[ij,]) %*% t(as.matrix(center_dat[(ij+lag),]))
            }
        }
        else
        {
            for(ij in 1:(T+lag))
            {
                gamma_lag_sum = gamma_lag_sum + as.matrix(center_dat[ij-lag,]) %*% t(as.matrix(center_dat[ij,]))
            }
        }
        return(gamma_lag_sum/T)
    }
    
    rho_val = vector("numeric", length(lag_value_range))
    for(ik in 1:length(lag_value_range))
    {
        lag_value = lag_value_range[ik] 
        if(lag_value == 0)
        {
            rho_val[ik] = NA
        }
        else
        {
            rho_val[ik] = sqrt(sum((gamma_l(lag = lag_value, T = T))^2))/sum(diag(gamma_l(lag = 0, T = T)))
        }
    }
    return(rho_val)
}
