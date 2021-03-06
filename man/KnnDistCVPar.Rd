% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ParallelKNNCrossValidation_functions.R
\name{KnnDistCVPar}
\alias{KnnDistCVPar}
\title{K-Nearest Neighbour correct cross-validation with distance input using parallel processing}
\usage{
KnnDistCVPar(
  DistMat,
  GroupMembership,
  K,
  EqualIter = 100,
  SampleSize = NA,
  TieBreaker = c("Random", "Remove", "Report"),
  Verbose = FALSE
)
}
\arguments{
\item{DistMat}{is a square matrix of pairwise distances among all reference specimens.}

\item{GroupMembership}{a character or factor vector in the same order as the distance data to denote group membership.}

\item{K}{is the number of nearest neighbours that the method will use for assigning group classification.}

\item{EqualIter}{sets the number of iterations resampling to equal sample size will be carried out.}

\item{SampleSize}{is the sample number that groups will be subsampled to if \code{Equal} is set to TRUE. The default is set to NA and will therefore use the smallest sample size of the groups provided.}

\item{TieBreaker}{is the method used to break ties if there is no majority resulting from K. Three methods are available('Random', 'Remove' and 'Report'): Random randomly returns one of tied classifications; Remove returns 'UnIDed' for the classification; Report returns a the multiple classifications as a single character string with tied classifications separated by '_'. NOTE: for correct cross-validation proceedures the results of both Report will be considered an incorrect identification even if one of the multiple reported classifications is correct.}

\item{Verbose}{determines whether the cross-validation results for each reference specimen is returned. Note that if this is set to TRUE and Equal is set to TRUE the funtion will return a list with the results of each iteration which will slow the process dramatically and take a lot of local memory.}
}
\value{
Returns a matrix of the leave-one-out classifications for all the specimens along with their known classification.
}
\description{
This function takes a square matrix of distances among specimens of known group membership
and returns the results of a leave-one-out correct cross validation identification for each
specimen to provide a correct cross-validation percentage.
}
\details{
The function is primarily for use with resampling unequal groups to equal sample size a set
number of times. This process is carried out with parrallel processing.

This function applies both a weighted approach and an unweighted approach and returns both results.

Note that this function is faster when datasets are large and/or when greater numbers of resampling
iterations are used. For small samples and few resampling iterations the function is unlikely to be
much faster, this is because in addition to the time it takes to carry out calculations the parallel
processing will need to compile the results at the end. This process adds additional time to the
process.
}
\author{
Ardern Hulme-Beaman
}
