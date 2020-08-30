
[![Build Status](https://travis-ci.com/AnthonyChristidis/splitSelect.svg?branch=master)](https://travis-ci.com/AnthonyChristidis/splitSelect) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/splitSelect)](https://cran.r-project.org/package=splitSelect) [![Downloads](http://cranlogs.r-pkg.org/badges/splitSelect)](https://cran.r-project.org/package=splitSelect)

splitSelect
============

This package provides functions for generating all possible splits of variables into groups, and computing the best split selection regression estimator for low-dimensional data.

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=splitSelect).

``` r
install.packages("splitSelect", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/splitSelect).

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/splitSelect")
```

### Usage

Here is some code to generate all possible splits of variables into groups.

``` r
# Loading library
library(splitSelect)

# Setting number of variables and groups
p <- 8
G <- 4
use.all <- TRUE

# Generate the number of partitions
my.partitions <- generate_partitions(p, G, use.all=use.all)
my.partitions

# Generate the number of splits
nsplit(p, G, use.all=use.all)
# Generate the number of splits (fixed partition)
nsplit(p, G, use.all=use.all,
       fix.partition=matrix(c(2,2,2,2), nrow=1))

# Generate the splits
all.splits <- generate_splits(p, G, use.all=use.all)
head(all.splits)
nrow(all.splits)
# Generate the splits (fixed partition)
all.splits <- generate_splits(p, G, use.all=use.all,
                              fix.partition=matrix(c(2,2,2,2), nrow=1))
head(all.splits)
nrow(all.splits)

# Generate samples of splits
sample.splits <- rsplit(10000, p, G, fix.partition=matrix(c(2,2,2,2), nrow=1))
sample.splits
```

Here is some code to apply to compute the best split selection estimator for simulated data with spurious correlation in the training set.

``` r
# Download the packages
install.packages("simTargetCov")
install.packages("glmnet")
install.packages("SplitReg")

# Setting the parameters
p <- 6
n <- 30
n.test <- 5000
group.beta <- 5
beta <- c(rep(1, 2), rep(group.beta, p-2))
rho <- 0.1
r <- 0.9
SNR <- 3
# Creating the target matrix with "kernel" set to rho
target_cor <- function(r, p){
  Gamma <- diag(p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      Gamma[i,j] <- Gamma[j,i] <- r^(abs(i-j))
    }
  }
  return(Gamma)
}
# AR Correlation Structure
Sigma.r <- target_cor(r, p)
Sigma.rho <- target_cor(rho, p)
sigma.epsilon <- as.numeric(sqrt((t(beta) %*% Sigma.rho %*% beta)/SNR))

# Number of cores
nb.cores <- parallel::detectCores()-1
# Registering the clusters
cl <- parallel::makeCluster(nb.cores)
doParallel::registerDoParallel(cl)

# Set the seed
set.seed(0)

# Simulate some data
x.train <- simTargetCov::simTargetCov(n=n, p=p, target=Sigma.r)
y.train <- 1 + x.train %*% beta + rnorm(n=n, mean=0, sd=sigma.epsilon)
x.test <- mvnfast::rmvn(n.test, mu=rep(0,p), sigma=Sigma.rho)
y.test <- 1 + x.test %*% beta + rnorm(n.test, sd=sigma.epsilon)

# Best Split Selection for Regression
system.time(
  split.out <- cv.splitSelect(x.train, y.train, G=2, use.all=TRUE,
                              fix.partition=list(matrix(c(2,4,
                                                          3,3), ncol=2, byrow=TRUE)), fix.split=NULL,
                              intercept=TRUE, group.model="glmnet", alpha=0, nfolds=10,
                              parallel=TRUE, cores=nb.cores)
)
split.predictions <- predict(split.out, newx=x.test)
mean((split.predictions-y.test)^2)/sigma.epsilon^2

# Ending the cluster
parallel::stopCluster(cl)

# Ridge Regression
cv.ridge <- glmnet::cv.glmnet(x.train, y.train, alpha=0)
ridge <- glmnet::glmnet(x.train, y.train, alpha=0, lambda=cv.ridge$lambda.min) 
ridge.predictions <- predict(ridge, newx=x.test)
mean((ridge.predictions-y.test)^2)/sigma.epsilon^2

# Lasso
cv.lasso <- glmnet::cv.glmnet(x.train, y.train, alpha=1)
lasso <- glmnet::glmnet(x.train, y.train, alpha=1, lambda=cv.lasso$lambda.min)
lasso.predictions <- predict(lasso, newx=x.test)
mean((lasso.predictions-y.test)^2)/sigma.epsilon^2

# Elastic Net
cv.elastic <- glmnet::cv.glmnet(x.train, y.train, alpha=3/4)
elastic <- glmnet::glmnet(x.train, y.train, alpha=3/4, lambda=cv.elastic$lambda.min)
elastic.predictions <- predict(elastic, newx=x.test)
mean((elastic.predictions-y.test)^2)/sigma.epsilon^2

# SplitReg
cv.splitreg <- SplitReg::cv.SplitReg(x.train, y.train, num_models=3, alpha=1e-2)
splitreg.predictions <- predict(cv.splitreg, newx=x.test)
mean((splitreg.predictions-y.test)^2)/sigma.epsilon^2

# Looking at the MSPEs for all the possible splits (out-of-sample)
split.mspes <-
  sapply(1:nrow(split.out$splits), function(x, n.test, x.test, split.out, y.test) 
  {mean((y.test-cbind(rep(1, n.test), x.test) %*% split.out$betas[,x])^2)}, 
  n.test, x.test, split.out, y.test)/sigma.epsilon^2
# Minimum MSPE for the splits (out-of-sample)
min(split.mspes)
split.mspes[split.out$optimal.split]
# Optimal splits comparison (out-of-sample)
split.out$splits[which.min(split.mspes),]
split.out$optimal.split.var
# Optimal betas comparison (out-of-sample)
split.out$betas[,which.min(split.mspes), drop=FALSE]
coef(split.out)
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
