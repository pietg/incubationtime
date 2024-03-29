# incubationtime

As a start we analyze the data of travelers from Wuhan, as analyzed by the Dutch Centre for Infectious Disease Control (RIVM): https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.5.2000062 A sample of 88 people was used in the estimation of the distribution of the incubation time for COVID-19. The RIVM used (low dimensional) parametric maximum likelihood to analyze these data using, e.g., the Weibull distribution as a possible distribution for the incubation time. We use nonparametric maximum likelihood instead.

The original data file is original_data_Wuhan.tsv. This was transformed into a data file transformed_data_Wuhan.txt. This, in turn, was transformed into the input file inputdata_Wuhan.txt, where the time, spent in Wuhan, was shifted making the left point equal to zero. For traveler number 67, who apparently had a connecting flight, the duration of stay in Wuhan was changed from 0 to 1 day.

Using the Weibull maximum likelihood method, the estimate was computed by two methods. One is a very simple method using Weibull.cpp where one does not have to compute the derivatives of the log likelihood and can use the Hooke-Jeeves pattern search method. The other one can be found in R_Weibull_Estimation.R, where we use the R package lbfgs, and where the gradient (derivatives of the log likelihood) has to be provided. The results are remarkably similar.

The parametric and nonparametric method are both present in analysis_EM.R and analysis_ICM.R, where we use the R package Rcpp, to connect to the files NPMLE_EM.cpp, NPMLE_ICM.cpp and Weibull.cpp, in which the nonparametric and parametric maximum likelihood estimators are computed. For instructions on how to use the R package Rcpp on Mac, Windows or Linux, see, e.g., https://teuder.github.io/rcpp4everyone_en/020_install.html. At the end of the scripts analysis_EM.R and analysis_ICM.R a picture of the difference in the estimation of the density by the two methods is shown. The estimates of the first moments are similar, though.

The nonparametric maximum likelihood estimator is computed by the EM algorithm and by the much more efficient iterative convex minorant algorithm. For the EM method, one has to run analysis_EM.R and for the iterative convex minorant algorithm one has to use analysis_ICM.R. The two methods give exactly the same result for the MLE, but the iterative convex minorant (ICM) algorithm needs much less iterations to get the result, with the necessary and sufficient (Fenchel) duality conditions satisfied at level 10^{-10}.

This repository complements the paper http://www.nieuwarchief.nl/serie5/pdf/naw5-2020-21-3-181.pdf in "Nieuw Archief voor Wiskunde" ("New Archive for Mathematics", a journal for the Dutch mathematicians), discussing the different ways of estimating the distribution of the incubation time. More details about the methods are given in incubation.pdf: https://onlinelibrary.wiley.com/doi/abs/10.1111/stan.12231

Note that the computation of the MLE does not need a choice of kernel or bandwidth. The bandwidths for the SMLE and density estimate can be chosen by a data-adaptive method, described in incubation.pdf, using the smoothed bootstrap. We demonstrate this method for the density estimate in analysis_ICM.R, using 1000 bootstrap samples for each point on a grid. The R script for this is: bandwidth_choice.R.

Usually the choice of the kernel does not make much difference, as long as the kernel is sufficiently smooth at the boundary (for this reason we have a preference for the triweight above the Epanechnikov kernel).

The directory Lauer_analysis contains an R script for the analysis of a data file of 181 subjects corresponding to the paper "The Incubation Period of Coronavirus Disease 2019 (COVID-19) From Publicly Reported Confirmed Cases: Estimation and Application", Lauer et al., Annals of Internal Medicine, 2020, 577--582. These data are *doubly* censored (the time of becoming symptomatic is only known to belong to an interval). 

The R script analysis_IQM_Lauer.R gives the MLE, SMLE and density estimate of the incubation time distribution for this data set and compares this with the results, using the Weibull and log-normal distributions, coming from the R package coarseDataTools. Moreover, the R script analysis_IQM_testdata.R gives the analysis for a test data set of size n=500, where the incubation times are generated from a mixture of normal distributions.

The directory "simulations" contains (for sample size n=1000, but the sample size can be varied) three types of simulations for the single-censored model and the statistics mean (first moment), median and 95th percentile. The estimates, derived from the nonparamtric MLE, were compared with the estimates, based on the Weibull and log-normal model, when the incubation time data are generated from the Weibull distribution. The log-normal estimates are clearly off for all these statistics, as shown in the simulations. Conversely, if we would have used the log-normal distribution to generate the incubation times, the Weibull estimates would have been off. In fact, the parametric estimates are inconsistent if the underlying distribution is of a different type.

The first moment is estimated by the mean of the nonparametric MLE (no smoothing is needed in this case). The estimates of the median and 95th percentile are based on the SMLE. R scripts are given in the directory "simulations". The parametric estimates are computed using the R package nloptr.

The directory "bootstrap" contains R scripts for the construction of bootstrap 95% confidence intervals for the distribution function and density of the Wuhan data. They are all derived from the nonparametric MLE. The smoothed bootstrap is used; the ordinary bootstrap is inconsistent for the MLE in this case.

The directory "support reduction" contains scripts, reproducing the pictures of the paper single_double_int_censoring.pdf (https://arxiv.org/abs/2310.04225). Presently it contains the R script reproducing Figures 5 and 7.


