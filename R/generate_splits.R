#' 
#' @title Generate Splits Possibilities
#'
#' @description \code{generate_splits} returns a matrix with the different splits of the variables in reach row.
#' 
#' @param p Number of variables or objects to split.
#' @param G Number of groups into which the variables are split.
#' @param use.all Boolean variable to determine if all variables must be used (default is TRUE).
#' @param fix.partition Optional matrix with G columns (or list if more than one value of G) indicating the partitions (in each row) to be considered for the splits.
#' @param verbose Boolean variable to determine if console output for cross-validation progress is printed (default is TRUE).
#' 
#' @return A matrix with the different splits of the variables in the groups.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Generating the possible splits of 6 variables in 3 groups
#' # Using all the variables
#' split.3groups.all <- generate_splits(6, 3)
#' split.3groups.all
#' # Without using all the variables
#' split.3groups <- generate_splits(6, 3, use.all=FALSE)
#' split.3groups
#' 
generate_splits <- function(p, G, use.all = TRUE,
                            fix.partition = NULL,
                            verbose=TRUE){
  
  # Case where G has more than one element
  if(length(G)>1){
    all.splits <- matrix(ncol=p, nrow=0)
    
    for(g in 1:length(G)){
      current.splits <- generate_splits(p, G[g], use.all=use.all, fix.partition=fix.partition[[g]])
      all.splits <- rbind(all.splits, current.splits)
    }
    return(all.splits)
  }
  
  # Check input data 
  if(any(c(!is.numeric(p), length(p)!=1, p<1, !(p%%1==0))))
    stop("p should be a positive interger.")
  if(any(c(!is.numeric(G), any(G<1), length(G)!=1,!any(G%%1==0), any(G>p))))
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
  
  # Determining if the computation should be done
  total.splits <- nsplit(p, G, use.all = use.all,
                         fix.partition = fix.partition)
  if(total.splits > 1e5)
    stop("There are over 100,000 total number of possible splits. The number of splits is over the limit.") else if(verbose)
      cat("Generating", total.splits, "total splits...\n")
  
  # Case where there is only one group (G=1)
  if(G==1){
    
    if(use.all)
      return(matrix(rep(1, p), ncol=p)) else{
        # Create matrix to store splits
        final.splits <- matrix(ncol=p, nrow=0)
        
        # Variables selection for the group
        variables.selection <- sapply(1:p, function(x) t(combn(p,x)))
        
        for(j in 1:length(variables.selection)){
          # Creating the matrix for the variables selection
          variables.selection.matrix <- matrix(NA, nrow(variables.selection[[j]]), ncol=p)
          # Assigning variables to group
          for(i in 1:nrow(variables.selection[[j]]))
            variables.selection.matrix[i,variables.selection[[j]][i,]] <- 1
          
          # Adding matrix to final splits matrix
          final.splits <- rbind(final.splits, variables.selection.matrix)
        }
        
        # Returning the final splits
        return(final.splits)
      }
  }
  
  # Generate the partitions
  if(is.null(fix.partition))
    partitions <- generate_partitions(p, G, use.all) else
      partitions <- fix.partition
  
  # Looping over the different partitions
  for(n.partition in 1:nrow(partitions)){
    
    # Storing the current partition 
    current.partition <- partitions[n.partition,]
    
    # Storing the different number of group sizes
    partition.new <- unique(current.partition)
    
    # Looping over the different number of group sizes
    for(partition.iter in 1:(length(partition.new))){
      split.combn <- t(combn(p,partition.new[partition.iter]))
    subsplits <- create.subSplits(split.combn, sum(current.partition==partition.new[partition.iter]))
    if(partition.iter==1)
      partition.splits <- subsplits else
        partition.splits <- merge.combinations(partition.splits, subsplits)
    }
    
    # Setting the groups
    partition.groups <- matrix(nrow=nrow(partition.splits), ncol=p)
    for(g in 1:G){
      if(g==1)
        g.start <- 1 else
          g.start <- sum(current.partition[1:(g-1)]) + 1
      g.end <- sum(current.partition[1:g])
      for(group.row in 1:nrow(partition.splits))
        partition.groups[group.row,partition.splits[group.row, g.start:g.end]] <- g
    }
    
    # Setting the final splits
    if(n.partition==1)
      final.splits <- partition.groups else
        final.splits <- rbind(final.splits, partition.groups)
  }
  
  # Returning the final splits
  return(final.splits)
}








