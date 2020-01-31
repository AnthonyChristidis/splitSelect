# -----------------------------------------------------------------------
# Object Construction for SPLIT object
# 
# object: the SPLIT object
# fn_call: the function call
# x: the design matrix
# y: the response vector
construct.SPLIT <- function(object, fn_call, x, y, intercept=TRUE){
  class(object) <- append("SPLIT", class(object))
  num_splits <- nrow(object$splits)
  mux_train <- apply(x, 2, mean)
  muy_train <- mean(y)
  if(intercept)
    object$intercepts <- object$betas[1,]
  object$call <- fn_call
  return(object)
}

# -----------------------------------------------------------------------
# Object Construction for cv.SPLIT object
# 
# object: the cv.SPLIT object
# fn_call: the function call
# x: the design matrix
# y: the response vector
construct.cv.SPLIT <- function(object, fn_call, x, y, intercept=TRUE){
  class(object) <- append("cv.SPLIT", class(object))
  num_splits <- nrow(object$splits)
  mux_train <- apply(x, 2, mean)
  muy_train <- mean(y)
  if(intercept)
    object$intercepts <- object$betas[1,]
  object$call <- fn_call
  return(object)
}
































