# simulations

(Added 10-28-2021) For sample size n=1000 three types of simulations were done for the single-censored model and the statistics mean (first moment), median and 95th percentile. The estimates, derived from the nonparamtric MLE, were compared with the estimates, based on the Weibull and log-normal model, when the incubation time data are generated from the Weibull distribution. The log-normal estimates are clearly off for all these statistics, as shown in the simulations. Conversely, if we would have used the log-normal distribution to generate the incubation times, the Weibull estimates would have been off. In fact, the parametric estimates are inconsistent if the underlying distribution is of a different type.

The first moment is estimated by the mean of the nonparatric MLE (no smoothing is needed in this case). The estimates of the median and 95th percentile are based on the SMLE. R scripts are given in the directory "simulations". The parametric estimates are computed using the R package nloptr.


