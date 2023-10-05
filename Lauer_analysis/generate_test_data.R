set.seed(43)

## Generate bimodal distribution
## as a mixture of two normal distributions
## using rejection sampling method

n <- 1000 # sample size
M <- 20  # right bound of distribution

## predefined density function
double_normal <- function(x) {
  res <- (0.8 * dnorm(x, 7, 2) + 0.2 * dnorm(x, 13, 1))/0.999813896704
}

## test plot 1
x <- seq(0, 20, by = 0.1)
plot(x, double_normal(x), type = "l")

## compute maximum of the density function
y <- double_normal(x)
ymax <- max(y)

## initialize the list for output
incubation_times <- c()

## Do rejection sampling
while (length(incubation_times) < n) {
  # ste1: generate random uniform from [0, M]
  U <- runif(1, max = M)
  
  # step2: sample a standard Uniform[0,1]
  V <- runif(1)
  
  # step3: check if we accept it
  if (V < double_normal(U) / (M * ymax)) {
    incubation_times <- c(incubation_times, U)
  }
}

incubation_times <- round(incubation_times)

V1 <- sample(1:30, n, replace=T)
V <- vector(,n)
for (i in 1: n){
	V[i] <- sample(1:V1[i],1)
	}
 

# left interval
left_ <- sample(1:3, n, replace=T)
right_ <- sample(1:3, n, replace=T)

tmp <- V + incubation_times

V2 <- tmp - left_
V3 <- tmp + right_

df <- data.frame(V1, V2, V3)

write.csv(df, "test_data.csv", row.names=FALSE)
