#' 
#' @importFrom stats glm
#' 
#' @title Split Selection for Regression - Coefficients Generation
#'
#' @description \code{splitSelect_coef} generates the coefficients for a particular split of variables into groups.
#' 
#' @param x Design matrix.
#' @param y Response vector.
#' @param variables.split A vector with the split of the variables into groups as values.
#' @param intercept Boolean variable to determine if there is intercept (default is TRUE) or not.
#' @param family Description of the error distribution and link function to be used for the model. Must be one of "gaussian" or "binomial".
#' @param group.model Model used for the groups. Must be one of "glmnet" or "LS".
#' @param lambdas The shinkrage parameters for the "glmnet" regularization. If NULL (default), optimal values are chosen.
#' @param alphas Elastic net mixing parameter. Should be between 0 (default) and 1.
#' 
#' @return A vector with the regression coefficients for the split.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Setting the parameters
#' p <- 6
#' n <- 30
#' n.test <- 5000
#' group.beta <- -3
#' beta <- c(rep(1, 2), rep(group.beta, p-2))
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
#' splitSelect_coef(x.train, y.train, variables.split=matrix(c(1,2,1,2,1,2), nrow=1))
#' 
splitSelect_coef <- function(x, y, variables.split, 
                             intercept=TRUE, 
                             family=c("gaussian", "binomial")[1],
                             group.model=c("glmnet", "LS", "Logistic")[1], 
                             lambdas=NULL, alphas=0){
  
  # Check input for the family (link)
  if(!(family %in% c("gaussian", "binomial"))){
    stop("group.model should be one of \"gaussian\" or \"binomial\".")
  }
  # Check input for lambdas and alphas
  if(!is.null(lambdas)){
    if(length(lambdas)!=length(unique(variables.split)))
      lambdas <- rep(lambdas[1], length(unique(variables.split)))
  }
  if(!is.null(alphas)){
    if(length(alphas)!=length(unique(variables.split)))
      alphas <- rep(alphas[1], length(unique(variables.split)))
  }
  
  # Storing the number of samples
  n <- nrow(x)
  # Storing the number of variables
  p <- ncol(x)

  # Storage of the beta coefficients
  final.beta <- numeric(p)
  
  # Storing the number of groups
  G <- unique(variables.split)
  
  # Centering the design matrix
  x.c <- scale(x, center=TRUE, scale=FALSE)
  # Centering the response
  y.c <- scale(y, center=TRUE, scale=FALSE)
  
  
  if(group.model=="LS"){
    
    if(intercept){
      
      for(g in G){
        x.g <- x.c[,variables.split==g, drop=FALSE]
        beta.g <- solve(t(x.g)%*%x.g)%*%t(x.g)%*%y.c
        final.beta[variables.split==g] <- beta.g
      }
      # Adjusting for the number of groups
      final.beta <- final.beta/length(G)
      # Computing the intercept
      split.intercept <- as.numeric(mean(y)) - apply(x, 2, mean)%*%final.beta
      final.beta <- c(split.intercept, final.beta)
    } else{
      
      for(g in G){
        x.g <- x[,variables.split==g, drop=FALSE]
        beta.g <- solve(t(x.g)%*%x.g)%*%t(x.g)%*%y
        final.beta[variables.split==g] <- beta.g
      }
      # Adjusting for the number of groups
      final.beta <- final.beta/length(G)
    }
    
    # Returning the adaptive SPLIT estimate
    return(final.beta)
    
  } else if(group.model=="Logistic"){
    
    if(intercept){
    
      # for(g in G){
      #   x.g <- x.c[,variables.split==g, drop=FALSE]
      #   beta.g <- solve(t(x.g)%*%x.g)%*%t(x.g)%*%y.c
      #   final.beta[variables.split==g] <- beta.g
      # }
      # # Adjusting for the number of groups
      # final.beta <- final.beta/length(G)
      # # Computing the intercept
      # split.intercept <- as.numeric(mean(y)) - apply(x, 2, mean)%*%final.beta
      # final.beta <- c(split.intercept, final.beta)
    } else{
      
      for(g in G){
        x.g <- x[,variables.split==g, drop=FALSE]
        beta.g <- glm(y ~ x.g, family=family)$coef
        final.beta[variables.split==g] <- beta.g
      }
      # Adjusting for the number of groups
      final.beta <- final.beta/length(G)
    }
    
    # Returning the adaptive SPLIT estimate
    return(final.beta)
    
  } else if(group.model=="glmnet"){

    if(is.null(lambdas)){

      # Computing the glmnet estimates for the variables
      for(g in G){
        if(sum(variables.split==g)==1){
          x.g <- cbind(0, x[, variables.split==g, drop=FALSE])
          beta.g <- glmnet::cv.glmnet(x.g, y, alpha=alphas[g], 
                                      intercept=intercept, grouped=FALSE, 
                                      family=family)
          beta.g <- as.numeric(coef(beta.g, s="lambda.min"))[-(1:2)]
        } else{
          x.g <- x[, variables.split==g, drop=FALSE]
          beta.g <- glmnet::cv.glmnet(x.g, y, alpha=alphas[g], 
                                      intercept=intercept, grouped=FALSE, 
                                      family=family)
          beta.g <- as.numeric(coef(beta.g, s="lambda.min"))[-1]
        }
        final.beta[variables.split==g] <- beta.g
      }
    } else{
      # Computing the glmnet estimates for the variables
      for(g in G){
        if(sum(variables.split==g)==1){
          x.g <- cbind(0, x[, variables.split==g, drop=FALSE])
          beta.g <- glmnet::glmnet(x.g, y, alpha=alphas[g], lambda=lambdas[g], 
                                   intercept=intercept, grouped=FALSE, 
                                   family=family)
          beta.g <- as.numeric(coef(beta.g))[-(1:2)]
        } else{
          x.g <- x[, variables.split==g, drop=FALSE]
          beta.g <- glmnet::glmnet(x.g, y.c, alpha=alphas[g], lambda=lambdas[g], 
                                   intercept=intercept, grouped=FALSE, 
                                   family=family)
          beta.g <- as.numeric(coef(beta.g))[-1]
        }
        final.beta[variables.split==g] <- beta.g
      }
    }
    
    # Adjusting for the number of groups
    final.beta <- final.beta/length(G)
    # Compute the intercept
    if(intercept){
      split.intercept <- as.numeric(mean(y) - apply(x, 2, mean)%*%final.beta)
      final.beta <- c(split.intercept, final.beta)
    }
    
    # Returning the adaptive SPLIT estimate
    return(final.beta)
  }
}
















