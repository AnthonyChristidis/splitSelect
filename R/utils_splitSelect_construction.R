# -----------------------------------------------------------------------
# Object Construction for splitSelect object
# 
# object: the splitSelect object
# fn_call: the function call
# x: the design matrix
# y: the response vector
construct.splitSelect <- function(object, fn_call, x, y, intercept=TRUE){
  class(object) <- append("splitSelect", class(object))
  num_splits <- nrow(object$splits)
  mux_train <- apply(x, 2, mean)
  muy_train <- mean(y)
  if(intercept)
    object$intercepts <- object$betas[1,]
  object$call <- fn_call
  return(object)
}

# -----------------------------------------------------------------------
# Object Construction for cv.splitSelect object
# 
# object: the cv.splitSelect object
# fn_call: the function call
# x: the design matrix
# y: the response vector
construct.cv.splitSelect <- function(object, fn_call, x, y, intercept=TRUE){
  class(object) <- append("cv.splitSelect", class(object))
  num_splits <- nrow(object$splits)
  mux_train <- apply(x, 2, mean)
  muy_train <- mean(y)
  if(intercept)
    object$intercepts <- object$betas[1,]
  object$call <- fn_call
  return(object)
}
































