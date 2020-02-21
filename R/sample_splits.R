#' 
#' @title Generate Samples of Splits Possibilities
#'
#' @description \code{rsplit} returns a matrix with random splits of the variables in groups.
#' 
#' @param n Number of sample splits.
#' @param p Number of variables or objects to split.
#' @param G Number of groups into which the variables are split.
#' @param use.all Boolean variable to determine if all variables must be used (default is TRUE).
#' @param fix.partition Optional matrix with G columns indicating the partitions (in each row) to be considered for the splits.
#' @param verbose Boolean variable to determine if console output for cross-validation progress is printed (default is TRUE).
#' 
#' @return A matrix or list with the number of possible objects in each group using splits.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Generating sample splits of 6 variables in 3 groups
#' # Using all the variables
#' random.splits <- rsplit(100, 6, 3)
#' # Using fixed partitions
#' random.splits.fixed <- rsplit(100, 6, 3, fix.partition=matrix(c(2,2,2), nrow=1))
#' 
rsplit <- function(n, p, G, use.all = TRUE,
                   fix.partition = NULL,
                   verbose=TRUE){
  
  # Case where G has more than one element
  if(length(G)>1){
    all.splits <- matrix(ncol=p, nrow=0)
    
    for(g in 1:length(G)){
      current.splits <- rsplit(n, p, G[g], use.all=use.all, fix.partition=fix.partition[[g]])
      all.splits <- rbind(all.splits, current.splits)
    }
    return(all.splits)
  }
  
  # Check input data 
  if(any(c(!is.numeric(n), length(n)!=1, n<1, !(n%%1==0))))
    stop("n should be a positive interger.")
  if(any(c(!is.numeric(p), length(p)!=1, p<1, !(p%%1==0))))
    stop("p should be a positive interger.")
  if(any(c(!is.numeric(G), length(G)!=1, any(G<1), !any(G%%1==0), any(G>=p))))
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
  
  # Case where there is only one group (G=1)
  if(G==1){
    
    # Case where there is only one group
    if(verbose)
      cat("Only a single split is possible for G=1, and is the only value returned.\n")
    return(matrix(rep(1,p), ncol=p))
    
    # Create matrix to store splits
    final.splits <- matrix(ncol=p, nrow=n)
    
    # Determining the number of variables selected
    if(is.null(fix.partition))
      n.vars <- sample(1:p, n, replace=TRUE) else
        n.vars <- sample(c(fix.partition), n, replace=TRUE)
    # Selecting the variables
    for(i in 1:n){
      final.splits[i,sample(1:p, n.vars[i], replace=FALSE)] <- 1
    }
    
    # Returning the final splits
    return(final.splits)
  }
  
  # Generate the partitions
  if(is.null(fix.partition))
    partitions <- generate_partitions(p, G, use.all) else
      partitions <- fix.partition
    
  # Create matrix to store splits
  final.splits <- matrix(ncol=p, nrow=n)
    
  # Generating the random splits
  for(j in 1:n){
    
    # Determine the partition
    current.partition <- partitions[sample(1:nrow(partitions), 1),]
    vars.randomization <- sample(1:p) 
    # Assigning the splits
    for(g in 1:G){
      if(g==1)
        g.start <- 1 else
          g.start <- sum(current.partition[1:(g-1)]) + 1
      g.end <- sum(current.partition[1:g])
      final.splits[j,vars.randomization[g.start:g.end]] <- g
    }
  }
  
  # Returning the final splits
  return(final.splits)
}








