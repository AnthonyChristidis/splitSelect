#' 
#' @import foreach
#' 
#' @title SPLIT Regression Modeling for Low-Dimensional Data - Cross-Validation
#'
#' @description \code{cv.SPLIT} performs the best split selection algorithm with cross-validation
#' 
#' @param x Design matrix.
#' @param y Response vector.
#' @param intercept Boolean variable to determine if there is intercept (default is TRUE) or not.
#' @param G Number of groups into which the variables are split. Can have more than one value.
#' @param alpha Elastic net mixing parameter. Should be between 0 (default) and 1.
#' @param group.model Model used for the groups. Must be one of "glmnet" or "LS".
#' @param nsample Number of sample splits for each value of G. If NULL, then all splits will be considered (unless there is overflow).
#' @param use.all Boolean variable to determine if all variables must be used (default is TRUE).
#' @param fix.partition Optional list with G elements indicating the partitions (in each row) to be considered for the splits.
#' @param fix.split Optional matrix with p columns indicating the groups (in each row) to be considered for the splits.
#' @param nfolds Number of folds for the cross-validation procedure.
#' @param parallel Boolean variable to determine if parallelization of the function. Default is FALSE.
#' @param cores Number of cores for the parallelization for the function.
#' 
#' @return An object of class cv.SPLIT.
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
#' # split.out <- cv.SPLIT(x.train, y.train, G=2, use.all=TRUE,
#' #                      fix.partition=list(matrix(c(2,2), 
#' #                                                ncol=2, byrow=TRUE)), 
#' #                      fix.split=NULL,
#' #                      intercept=TRUE, group.model="glmnet", alpha=0, nfolds=10)
#'
cv.SPLIT <- function(x, y, intercept = TRUE,
                     G, use.all = TRUE,
                     group.model=c("glmnet", "LS")[1], alpha = 0,
                     nsample = NULL, fix.partition = NULL, fix.split = NULL,
                     nfolds = 10,
                     parallel=FALSE, cores=getOption('mc.cores', 2L)){
  
  # Check input data
  if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
    stop("x should belong to one of the following classes: matrix, data.frame.")
  } else if (all(!inherits(y, "matrix"), all(!inherits(y, "numeric")))) {
    stop("y should belong to one of the following classes: matrix, numeric.")
  } else if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
    stop("x should not have missing, infinite or nan values.")
  } else if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
    stop("y should not have missing, infinite or nan values.")
  } else {
    if(inherits(y, "matrix")) {
      if (ncol(y)>1){
        stop("y should be a vector")
      }
      y <- as.numeric(y)
    }
    len_y <- length(y)
    if (len_y != nrow(x)) {
      stop("y and x should have the same number of rows.")
    }
  }
  if(!is.null(alpha)){
    if (!inherits(alpha, "numeric")) {
      stop("alpha should be numeric")
    } else if (any(alpha < 0, alpha > 1)) {
      stop("alpha should be a numeric value between 0 and 1.")
    }
  }
  if(!(group.model %in% c("glmnet", "LS"))){
    stop("group.model should be one of \"glmnet\" or \"LS\".")
  }
  p <- ncol(x) # Storing the number of variables
  
  # Getting the full adaptive SPLIT estimate
  out.split <- SPLIT(x, y, intercept=intercept,
                     G=G, use.all=use.all,
                     group.model=group.model, alpha=alpha,
                     nsample=nsample, fix.partition=fix.partition, fix.split=fix.split,
                     parallel=parallel, cores=cores)
  
  # Creating the folds
  folds <- caret::createFolds(1:nrow(x), nfolds)

  if(parallel){
      
    # Registering cores if not already declared
    if (!foreach::getDoParRegistered()){
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)
    }
    
    # Determining the splits to store for each core
    parallel.id <- split(1:nrow(out.split$splits), factor(sort(rank(1:nrow(out.split$splits))%%cores)))
    
    # Hack initialization
    subset.ind <- NULL
    # Parallel computation for the subsets
    splits.mspes <- foreach::foreach(subset.ind=1:min(cores, length(parallel.id)), 
                                     .packages=c("SPLIT"), 
                                     .export=c("SPLIT", "SPLIT_generate_coefficients",
                                               "construct.SPLIT")) %dopar% { 
      
      # Data to store the mspes for the core
      core.splits <- parallel.id[[subset.ind]]
      core.splits.mspes <- numeric(length(core.splits))
                                  
      # Looping over the different splits
      for(split.id in 1:length(core.splits)){
        
        # Looping over the folds
        for(fold.id in 1:nfolds){
          x.train <- x[-folds[[fold.id]],]; x.test <- x[folds[[fold.id]],, drop=FALSE]
          y.train <- y[-folds[[fold.id]]]; y.test <- y[folds[[fold.id]]]
          fold.split <- SPLIT(x=x.train, y=y.train, intercept=intercept,
                              G=G, group.model=group.model, alpha = alpha,
                              fix.split = matrix(out.split$splits[core.splits[split.id],], nrow=1))
          if(intercept)
            split.pred <- cbind(rep(1, length(y.test)), x.test) %*% fold.split$betas else
              split.pred <- x.test %*% fold.split$betas
          core.splits.mspes[split.id] <- core.splits.mspes[split.id] + mean((y.test - split.pred)^2)
        }
        # Adjust mspes for number of folds
        core.splits.mspes[split.id] <- core.splits.mspes[split.id]/nfolds
      }
      # Return the core.splits.mspes
      return(core.splits.mspes)
    }
    # Unlisting the MSPEs for all the splits
    splits.mspes <- unlist(splits.mspes)
    
  } else{
    
    # Variable to store the MSPEs
    splits.mspes <- numeric(nrow(out.split$splits))
    
    # Looping over the different splits
    for(split.id in 1:nrow(out.split$splits)){
      
      # Looping over the folds
      for(fold.id in 1:nfolds){
        x.train <- x[-folds[[fold.id]],]; x.test <- x[folds[[fold.id]],]
        y.train <- y[-folds[[fold.id]]]; y.test <- y[folds[[fold.id]]]
        fold.split <- SPLIT(x=x.train, y=y.train, intercept=intercept,
                            G=G, group.model=group.model, alpha = alpha,
                            fix.split = matrix(out.split$splits[split.id,], nrow=1))
        if(intercept)
          split.pred <- cbind(rep(1, length(y.test)), x.test) %*% fold.split$betas else
            split.pred <- x.test %*% fold.split$betas
        splits.mspes[split.id] <- splits.mspes[split.id] + mean((y.test - split.pred)^2)
      }
      # Adjust mspes for number of folds
      splits.mspes[split.id] <- splits.mspes[split.id]/nfolds
    }
  }
  
  # Construct the output
  cv.SPLIT.out <- list(betas=out.split$betas, splits=out.split$splits, 
                       optimal.split=which.min(splits.mspes), splits.mspes=splits.mspes,
                       optimal.split.var=out.split$splits[which.min(splits.mspes),])
  fn.call <- match.call()
  cv.SPLIT.out <- construct.cv.SPLIT(object=cv.SPLIT.out, fn_call=fn.call, x=x, y=y, intercept=intercept)
  
  # Returning the splits and the coefficients
  return(cv.SPLIT.out)
}








