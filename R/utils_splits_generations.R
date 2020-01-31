# -----------------------------------------------------------------------
# Function for merging combinations of splits
# 
# first.data: first argument for matrix (of any type) for the splits
# second.data: second argument for matrix (of any type) for the splits
merge.combinations <- function(first.data, second.data){
  
  # Creating the matrix by merging the combinations from each group
  final.data <- first.data[rep(1:(nrow(first.data)), each=nrow(second.data)),]
  final.data <- cbind(final.data, second.data[rep(1:nrow(second.data), times=nrow(first.data)),])
  
  # Removing split groupings where a variable is duplicated (i.e. in more than one group)
  final.data <- final.data[apply(final.data, 1, function(x){return(!any(duplicated(x)))}),]
  
  # Returning the split groupings
  final.data <- unname(final.data)
  return(final.data)
}

# Function to check if object is interger(0)
is.integer0 <- function(x){
  return(is.integer(x) && length(x) == 0L)
}

# -----------------------------------------------------------------------
# Function for creating splits for a partition
# 
# comb.data: data
# second.data: second argument for matrix (of any type) for the splits
create.subSplits <- function(comb.data, n.times){
  
  # Storing the number of variables
  p <- ncol(comb.data)
  
  # Indicator for the groups
  p.ind <- c(1, (1:(n.times-1))*(p) + 1)
  
  # Initiation the merging
  new.comb.data <- comb.data
  
  # Looping for the number of times the group is repeated
  if(n.times>1){
    for(n.iter in 2:n.times){
      new.comb.data <- merge.combinations(new.comb.data, comb.data)
      new.comb.data <- new.comb.data[(new.comb.data[,p.ind[n.iter]] > new.comb.data[,p.ind[n.iter-1]]),]
    }
    return(new.comb.data)
  } else
    return(new.comb.data)
}






















