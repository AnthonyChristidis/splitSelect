% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.splitSelect.R
\name{coef.splitSelect}
\alias{coef.splitSelect}
\title{Coefficients for splitSelect object}
\usage{
\method{coef}{splitSelect}(object, ...)
}
\arguments{
\item{object}{An object of class splitSelect.}

\item{...}{Additional arguments for compatibility.}
}
\value{
A matrix with the coefficients of the \code{splitSelect} object.
}
\description{
\code{coef.splitSelect} returns the coefficients for a splitSelect object.
}
\examples{
# Setting the parameters
p <- 4
n <- 30
n.test <- 5000
beta <- rep(5,4)
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
sigma.epsilon <- as.numeric(sqrt((t(beta) \%*\% Sigma.rho \%*\% beta)/SNR))
# Simulate some data
x.train <- mvnfast::rmvn(30, mu=rep(0,p), sigma=Sigma.r)
y.train <- 1 + x.train \%*\% beta + rnorm(n=n, mean=0, sd=sigma.epsilon)

# Generating the coefficients for a fixed partition of the variables
\donttest{
split.out <- splitSelect(x.train, y.train, G=2, use.all=TRUE,
                         fix.partition=list(matrix(c(2,2), ncol=2, byrow=TRUE)), fix.split=NULL,
                         intercept=TRUE, group.model="glmnet", alphas=0)     
coef(split.out)
}

}
\seealso{
\code{\link{splitSelect}}
}
\author{
Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
}
