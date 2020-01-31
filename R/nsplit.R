#' 
#' @title Compute Total Number of Possible Splits
#'
#' @description \code{nsplits} returns the total number of possible splits of variables into groups.
#' 
#' @param p Number of variables or objects to split.
#' @param G Number of groups into which the variables are split. 
#' @param use.all Boolean variable to determine if all variables must be used (default is TRUE).
#' @param fix.partition Optional matrix with G columns (or list if more than one value of G) indicating the partitions (in each row) to be considered for the splits.
#' 
#' @return A numeric vector with the total number of possible splits.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Compute the total number of possible splits of 6 variables into 3 groups
#' # We use all the variables
#' out.n.splits.all <- nsplit(p=6, G=3, use.all=TRUE)
#' out.n.splits.all
#' # We don't enforce using all the variables
#' out.n.splits <- nsplit(p=6, G=3, use.all=FALSE)
#' out.n.splits
#' 
nsplit <- function(p, G, use.all = TRUE,
                   fix.partition = NULL){
  
  # Case where G has more than one element
  if(length(G)>1){
    n.splits <- numeric(length(G))
    
    for(g.ind in 1:length(G)){
      n.splits[g.ind] <- nsplit(p, G[g.ind], use.all=use.all, fix.partition=fix.partition[[g.ind]])
    }
    return(n.splits)
  }
  
  # Check input data 
  if(any(c(!is.numeric(p), length(p)!=1, p<1, !(p%%1==0))))
    stop("p should be a positive interger.")
  if(any(c(!is.numeric(G), any(G<1), !any((G%%1==0)), any(G>p))))
    stop("G should be a positive interger less or equal to p.")
  if(!is.null(fix.partition)){
    if(any(!is.matrix(fix.partition), ncol(fix.partition)!=G))
      stop("fix.partition should be a matrix with G columns")
    fix.sum <- rowSums(fix.partition)
    if(use.all)
      if(any(fix.sum!=p))
        stop("The number of variables used does not sum to p. Set use.all to FALSE if needed.") else
          if(any(fix.sum>p))
            stop("Some fixed partitions have an invalid number of total variables")
  }

  # Case where G is equal to p
  if(G==p)
    return(1)
  
  # Vector to store the number of splits
  n.splits <- numeric(1)
  
  # Looping over the number of groups G
  for(g in G){
    # Generate the partitions
    if(is.null(fix.partition))
      partitions <- generate_partitions(p, G, use.all) else
        partitions <- fix.partition
    
    # Number of splits
    for(n.partition in 1:nrow(partitions)){
      n.splits <- n.splits + multicool::multinom(c(partitions[n.partition,],p-sum(partitions[n.partition,])), counts=TRUE)/
        prod(factorial(table(partitions[n.partition,])))
    }
  }

  # Return the total number of possible splits
  return(n.splits)
}




