# -----------------------------------------------------------------------
# Object Construction for splitSelect object
# 
# object: the splitSelect object
# fn_call: the function call
# x: the design matrix
# y: the response vector
# intercept: boolean for intercept or not
# family: the family for the errors and link function
construct.splitSelect <- function(object, fn_call, x, y, intercept=TRUE, family=family){
  class(object) <- append("splitSelect", class(object))
  object$num_splits <- nrow(object$splits)
  if(intercept)
    object$intercepts <- object$betas[1,]
  object$call <- fn_call
  object$family <- family
  return(object)
}

# -----------------------------------------------------------------------
# Object Construction for cv.splitSelect object
# 
# object: the cv.splitSelect object
# fn_call: the function call
# x: the design matrix
# y: the response vector
# intercept: boolean for intercept or not
# family: the family for the errors and link function
construct.cv.splitSelect <- function(object, fn_call, x, y, intercept=TRUE, family=family){
  class(object) <- append("cv.splitSelect", class(object))
  object$num_splits <- nrow(object$splits)
  if(intercept)
    object$intercepts <- object$betas[1,]
  object$call <- fn_call
  object$family <- family
  return(object)
}
































