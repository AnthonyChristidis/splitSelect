#' 
#' @title Predictions for splitSelect object
#'
#' @description \code{predict.splitSelect} returns the prediction for splitSelect for new data.
#' 
#' @param object An object of class splitSelect.
#' @param newx A matrix with the new data.
#' @param ... Additional arguments for compatibility.
#' 
#' @return A matrix with the predictions of the \code{splitSelect} object.
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
#' x.test <- mvnfast::rmvn(n.test, mu=rep(0,p), sigma=Sigma.rho)
#' y.test <- 1 + x.test %*% beta + rnorm(n.test, sd=sigma.epsilon)
#' 
#' # Generating the coefficients for a fixed split
#' \donttest{
#' split.out <- splitSelect(x.train, y.train, G=2, use.all=TRUE,
#'                          fix.partition=list(matrix(c(2,2), ncol=2, byrow=TRUE)), fix.split=NULL,
#'                          intercept=TRUE, group.model="glmnet", alphas=0)
#' predict(split.out, newx=x.test)
#' }
#' 
#' @seealso \code{\link{splitSelect}}
#' 
predict.splitSelect <- function(object, newx, ...){
  
  # Check input data
  if(!any(class(object) %in% "splitSelect"))
    stop("The object should be of class \"splitSelect\"")
  if(is.matrix(newx)){
    if(ncol(newx)!=ncol(object$splits))
      stop("The dimension of newx is invalid.")
  } else if(length(newx)!=ncol(object$splits))
    stop("The number of variables for newx is invalid.")
  
  # Matrix to store the predictions
  predictions <- matrix(nrow=nrow(newx), ncol=nrow(object$splits))
  
  # Regression Case
  if(object$family=="gaussian"){
    
    # Removing the intercepts
    if(!is.null(object$intercepts))
      object$betas <- object$betas[-1,, drop=FALSE]
    
    # Computing the predictions
    for(newx.id in 1:nrow(newx)){
      for(split.id in 1:nrow(object$splits)){
        predictions[newx.id, split.id] <- newx[newx.id,] %*% object$betas[,split.id, drop=FALSE]
      }
    }
    
    # Adding the intercepts
    if(!is.null(object$intercepts))
      for(split.id in 1:nrow(object$splits))
        predictions[,split.id] <- predictions[,split.id] + object$intercepts[split.id]
      
  } else if(object$family=="binomial"){
    
    # Computing the predictions
    for(newx.id in 1:nrow(newx)){
      for(split.id in 1:nrow(object$splits)){
        predictions[newx.id, split.id] <- exp(newx[newx.id,] %*% object$betas[,split.id, drop=FALSE])/
          (1+exp(newx[newx.id,] %*% object$betas[,split.id, drop=FALSE]))
      }
    }
    # Getting the binary prediction
    predictions <- round(predictions, 0)
  }

  # Returning the coefficients
  return(predictions)
}








