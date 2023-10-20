# Confidence intervals

The R script analysis_SR_single.R produces Figure 5 in single_double_int_censoring.pdf, which gives confidence intervals and coverages of 1000 samples of size 1000 for the singly interval censoring model. If one runs the script, one will see the numbers of the samples drawn go by and at the end two pictures, one of the confidence interval for the last sample and one for the coverage of the confidence intervals of the real parameters for the points 3 to 10 over the 1000 samples.

The R script analysis_SR_doubly.R produces a figure analagous to Figure 7 in single_double_int_censoring.pdf, which gives confidence intervals and coverages of 1000 samples of size 1000 for the doubly interval censoring model. In this case the means of the diagonals of the Fisher information matrices over 1000 samples were used as estimates of the variances, to illustrate the validity of the asymptotic theory. In practice on would use bootstrap samples to this purpose.

All the confidence intervals of the paper were constructed by a version of the present support reduction algorithm. The manuscript can also be found on arXiv: https://arxiv.org/abs/2310.04225

For instructions on how to use the R package Rcpp on Mac, Windows or Linux, see, e.g., https://teuder.github.io/rcpp4everyone_en/020_install.html.


