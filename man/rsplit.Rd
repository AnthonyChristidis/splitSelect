% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_splits.R
\name{rsplit}
\alias{rsplit}
\title{Generate Samples of Splits Possibilities}
\usage{
rsplit(n, p, G, use.all = TRUE, fix.partition = NULL, verbose = TRUE)
}
\arguments{
\item{n}{Number of sample splits.}

\item{p}{Number of variables or objects to split.}

\item{G}{Number of groups into which the variables are split.}

\item{use.all}{Boolean variable to determine if all variables must be used (default is TRUE).}

\item{fix.partition}{Optional matrix with G columns indicating the partitions (in each row) to be considered for the splits.}

\item{verbose}{Boolean variable to determine if console output for cross-validation progress is printed (default is TRUE).}
}
\value{
A matrix or list with the number of possible objects in each group using splits.
}
\description{
\code{rsplit} returns a matrix with random splits of the variables in groups.
}
\examples{
# Generating sample splits of 6 variables in 3 groups
# Using all the variables
random.splits <- rsplit(100, 6, 3)
# Using fixed partitions
random.splits.fixed <- rsplit(100, 6, 3, fix.partition=matrix(c(2,2,2), nrow=1))

}
\author{
Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
}
