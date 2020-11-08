# incubationtime

As a start we analyze the data of travelers from Wuhan, as analyzed by the Dutch Centre for Infectious Disease Control (RIVM): https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.5.2000062 A sample of 88 people was used in the estimation of the distribution of the incubation time for COVID-19. The RIVM used (low dimensional) parametric maximum likelihood to analyze these data using, e.g., the Weibull distribution as a possible distribution for the incubation time. We use nonparametric maximum likelihood instead.

The original data file is original_data_Wuhan.tsv. This was transformed into a data file transformed_data_Wuhan.txt. This, in turn, was transformed into the input file inputdata_Wuhan.txt, where the time, spent in Wuhan, was shifted making the left point equal to zero. For traveler number 67, who apparently had a connecting flight, the duration of stay in Wuhan was changed from 0 to 1 day.

Using the Weibull maximum likelihood method, the estimate was computed by two methods. One is a very simple method using Weibull.cpp where one does not have to compute the derivatives of the log likelihood and can use the Hooke-Jeeves pattern search method. The other one can be found in R_Weibull_Estimation.R, where we use the R package lbfgs, and where the gradient (derivatives of the log likelihood) has to be provided. The results are remarkably similar.

The parametric and nonparametric method are both present in analysis_EM.R and analysis_ICM.R, where we use the R package Rcpp, to connect to the files NPMLE_EM.cpp, NPMLE_ICM.cpp and Weibull.cpp, in which the nonparametric and parametric maximum likelihood estimators are computed. For instructions on how to use the R package Rcpp on Mac, Windows or Linux, see, e.g., https://teuder.github.io/rcpp4everyone_en/020_install.html. At the end of the scripts analysis_EM.R and analysis_ICM.R a picture of the difference in the estimation of the density by the two methods is shown. The estimates of the first moments are similar, though.

The nonparametric maximum likelihood estimator is computed by the EM algorithm and by the much more efficient iterative convex minorant algorithm. For the EM method, one has to run analysis_EM.R and for the iterative convex minorant algorithm one has to use analysis_ICM.R. The two methods give exactly the same result for the MLE, but the iterative convex minorant (ICM) algorithm needs much less iterations to get the result, with the necessary and sufficient (Fenchel) duality conditions satisfied at level 10^{-10}.

This repository complements the paper http://www.nieuwarchief.nl/serie5/pdf/naw5-2020-21-3-181.pdf in "Nieuw Archief voor Wiskunde" ("New Archive for Mathematics", a journal for the Dutch mathematicians), discussing the different ways of estimating the distribution of the incubation time. More details about the methods are given in incubation.pdf.

Note that the computation of the MLE does not need a choice of kernel or bandwidth.The bandwidths for the SMLE and density estimate can be chosen by a data-adaptive method, described in incubation.pdf, using the smoothed bootstrap. We demonstrate this method for the density estimate in analysis_ICM.R, using 1000 bootstrap samples for each point on a grid on [3,7]. Usually the choice of the kernel does not make much difference, as long as the kernel is sufficiently smooth at the boundary (for this reason we have a preference for the triweight above the Epanechnikov kernel).

The directory incubation_continuous contains two files: a C++ file main.cpp and a makefile, which can be run in Linux or on a Mac and produces an executable called incubation_continuous.exe. Underlying theory is given in Section 4 and the appendix of incubation.pdf in the main directory of this repository.

The executable can be run from a terminal window, typing ./incubation_continuous.exe. It produces a file with 1000 lines of 10 variances of the density estimate in the continuous incubation time model. In this model the data are not discretized to days, and the simulations produce the nonparametric estimates of the density of the incubation time at 2,3,... 11 days. The executable also produces the example of an input file for the computation of the density estimate produced by the simulation (data_exp.txt) and the MLE and its masses (MLE.txt and ICM_masses, resp.).

The theory predicts the values of the variances of the estimates of the density, which can be found by running inteq.R in the main directory, which is presently set to evaluate the theoretic value at the point (number of days) 6, but can also be used to evaluate at another number of days between 2 and 11 or a number between these bounds that is not an integer. The executable was tested on a Mac and using Linux.
