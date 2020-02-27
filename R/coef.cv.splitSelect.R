#' 
#' @title Coefficients for splitSelect object
#'
#' @description \code{coef.cv.splitSelect} returns the coefficients for a cv.splitSelect for new data.
#' 
#' @param object An object of class cv.splitSelect.
#' @param optimal.only A boolean variable (TRUE default) to indicate if only the coefficient of the optimal split are returned.
#' @param ... Additional arguments for compatibility.
#' 
#' @return A matrix with the coefficients of the \code{cv.splitSelect} object.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Setting the parameters
#' p <- 4
#' n <- 30
#' n.test <- 5000
#' beta <- rep(5,4)
#' rho <- 0.1
#' r <- 0.9
#' SNR <- 3
#' # Creating the target matrix with "kernel" set to rho
#' target_cor <- function(r, p){
#'   Gamma <- diag(p)
#'   for(i in 1:(p-1)){
#'     for(j in (i+1):p){
#'       Gamma[i,j] <- Gamma[j,i] <- r^(abs(i-j))
#'     }
#'   }
#'   return(Gamma)
#' }
#' # AR Correlation Structure
#' Sigma.r <- target_cor(r, p)
#' Sigma.rho <- target_cor(rho, p)
#' sigma.epsilon <- as.numeric(sqrt((t(beta) %*% Sigma.rho %*% beta)/SNR))
#' # Simulate some data
#' x.train <- mvnfast::rmvn(30, mu=rep(0,p), sigma=Sigma.r)
#' y.train <- 1 + x.train %*% beta + rnorm(n=n, mean=0, sd=sigma.epsilon)
#' 
#' # Generating the coefficients for a fixed split
#' \donttest{
#' split.out <- cv.splitSelect(x.train, y.train, G=2, use.all=TRUE,
#'                             fix.partition=list(matrix(c(2,2), 
#'                                                ncol=2, byrow=TRUE)), 
#'                             fix.split=NULL,
#'                             intercept=TRUE, group.model="glmnet", alphas=0, nfolds=10)
#' coef(split.out)
#' }
#' 
#' @seealso \code{\link{splitSelect}}
#' 
coef.cv.splitSelect <- function(object, optimal.only = TRUE, ...){
  
  # Check input data
  if(!any(class(object) %in% "cv.splitSelect"))
    stop("The object should be of class \"cv.splitSelect\".")
  
  # Returning the coefficients
  if(optimal.only)
    return(object$betas[,object$optimal.split, drop=FALSE]) else
      return(object$betas)
}
