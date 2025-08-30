# ftsatoo

[![Documentation](https://img.shields.io/badge/documentation-is_here-green)](https://techtonique.github.io/ftsatoo/index.html)


Derived from R package ftsa. 

Functions for visualizing, modeling, (generic) forecasting and hypothesis testing of functional time series.

```R
library(ftsatoo)
library(ahead)

fit_obj <- ftsm(ElNino_ERSST_region_1and2, order=3, method = "M")
summary(fit_obj)
plot(fit_obj)

# Here's the novelty
(fcast <- forecast(object = fit_obj, h = 12, level=95, 
                    pimethod = "nonparametric",
                    FUN=forecast::thetaf))   
plot(fcast)                    
(fcast <- forecast(object = fit_obj, h = 12, level=95, 
                    pimethod = "nonparametric",
                    FUN=ahead::dynrmf))     
plot(fcast)                                      
```