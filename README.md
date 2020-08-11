# incubationtime

As a start we analyze the data of travelers of Wuhan, as analyzed by the Dutch Centre for Infectious Disease Control (RIVM): https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.5.2000062
A sample of 88 people was used in the estimation of the distribution of the incubation time for Covid 19. The RIVM used parametric maximum likelihood to analyze these data using, e.g., the Weibull distribution as a possible distribution for the incubation time. We use nonparametric maximum likelihood instead.

The original data file is original_data_Wuhan.tsv. This was transformed into a data file ``transformed_data_Wuhan''. This, in turn, was transformed into the input file inputdata_Wuhan.txt, where the time, spent in Wuhan, was shifted making the left point equal to zero. For traveler number 67, who apparently had a connecting flight, the duration of stay in Wuhan was changed from 0 to 1 day.

Using the Weibull maximum likelihood method, the estimate was computed by two methods. One is a very simple method using Weibull.cpp where one does not have to compute the derivatives of the log likelihood and can use the Hooke-Jeeves pattern search method. The other one can be found in R_Weibull_Estimation.R, where we use the R package lbfgs, and where the gradient (derivatives of the log likelihood) have to be provided. The results are remarkably similar.
The parametric and nonparametric method are both present in analysis.R, where we use the R package Rcpp, to connect to the files NPMLE.cpp and Weibul.cpp, in which the nonparametric and parametric maximum likelihood estimators are computed. To use the R package Rcpp, one either has to use Rtools (on Windows) or a terminal version of Xcode (on a Mac). At the end of the script analysis.R a picture of the striking difference in the estimation of the density by the two methods is shown. The estimates of the first moments is similar, though.

This repository complements a paper in "Nieuw Archief voor Wiskunde" ("New archive of Mathematics", a journal for the Dutch mathematicians), discussing the different ways of estimating the distribution of the incubation time.
