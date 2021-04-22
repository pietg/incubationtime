# Analysis Lauer data and a test for doubly censored data

This directory contains an R script for the analysis of a data file of 181 subjects corresponding to the paper "The Incubation Period of Coronavirus Disease 2019 (COVID-19) From Publicly Reported Confirmed Cases: Estimation and Application", Lauer et al., Annals of Internal Medicine, 2020, 577--582. These data are *doubly* censored (the time of becoming symptomatic is only known to belong to an interval). 
The R script analysis_IQM_Lauer.R gives the MLE, SMLE and density estimate of the incubation time distribution for this data set and compares this with the results, using the Weibull and log-normal distributions, coming from the R package coarseDataTools. Moreover, the R script analysis_IQM_testdata.R gives the analysis for a test data set of size n=500, where the incubation times are generated from a mixture of normal distributions. The R script for generating the dtata from the mixture of normals was written by Slavik Koval.


