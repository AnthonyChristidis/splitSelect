#' 
#' @importFrom stats coef rbinom
#' @importFrom utils combn
#' 
#' @title Generate Splits Partitions Possibilities
#'
#' @description \code{generate_partitions} returns a matrix with the number of possible objects in each group using splits.
#' 
#' @param p Number of variables or objects to split.
#' @param G Number of groups into which the variables are split.
#' @param use.all Boolean variable to determine if all variables must be used (default is TRUE).
#' 
#' @return A matrix or list with the number of possible objects in each group using splits.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Generating the possible split partitions of 6 variables in 3 groups
#' # Using all the variables
#' split.3groups.all <- generate_partitions(6, 3)
#' split.3groups.all
#' # Without using all the variables
#' split.3groups <- generate_partitions(6, 3, use.all=FALSE)
#' split.3groups
#' 
generate_partitions <- function(p, G, use.all = TRUE){
  
  # Case where G has more than one element
  if(length(G)>1){
    group.possibilities <- list()
    
    for(g in 1:length(G)){
      group.possibilities[[g]] <- generate_partitions(p, G[g], 
                                                      use.all=use.all)
    }
    return(group.possibilities)
  }
  
  # Check input data 
  if(any(c(!is.numeric(p), length(p)!=1, p<1, !(p%%1==0))))
    stop("p should be a positive interger.")
  if(any(c(!is.numeric(G), length(p)!=1, any(G<1), !any(G%%1==0), any(G>p))))
    stop("G should be a positive interger less or equal to p.")
  
  # Case where there is only one group (G=1)
  if(G==1){
    if(use.all)
      return(matrix(p, ncol=1)) else
        return(matrix(1:(p-1), ncol=1))
  }
  
  # Modifying G for extra group if not all variables must be used
  if(!use.all)
    G <- G + 1
  
  # Generating markers for the groups
  group.marker <- t(combn(1:(p+(G-1)), G-1))
  if((G-1)!=1)
    group.marker <- group.marker[c(group.marker[,1]!=1),] else
      group.marker <- matrix(group.marker[c(group.marker[,1]!=1),], ncol=1) # Case where this could be a vector
  
  # Generating the possibilities for the number of variables in each group
  group.possibilities <- as.matrix(group.marker[,1]-1, ncol=1)
  if(G>2){ # Case where G>2 for the subsequent columns
    for(g in 2:(G-1)){
      group.possibilities <- cbind(group.possibilities, 
                                   group.marker[,g]-(group.marker[,g-1]+1))
    }
  }
  rm(group.marker)
  # Adding the last group with the number of variables
  group.possibilities <- cbind(group.possibilities, p - rowSums(group.possibilities))
  # Removing cases where there is an empty group
  if(!use.all)
    group.possibilities <- group.possibilities[!(apply(group.possibilities[,-G], 1, function(t) return(any(t==0)))),][,-G] else
      group.possibilities <- group.possibilities[!(apply(group.possibilities, 1, function(t) return(any(t==0)))),]
  # Removing duplicated cases
  group.possibilities <- unique(t(apply(group.possibilities, 1, function(t) return(sort(t)))))
  
  # Returning the possible number of variables per group
  return(group.possibilities)
}








