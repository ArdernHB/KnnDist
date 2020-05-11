



#' Function for resampling groups to equal sample size
#'
#' This function takes a character of factor vector and returns a new vector of equal sample size
#' along with a vector of index positions to sort whatever other data might be associated with the
#' factor vector into the new equal sample size order. The user can also provide a user defined minimum
#' sample size if the user wishes the resampling procedure to sample to a smaller number of specimens
#' than the smallest sample size. In this way the function can also be used in a leave-p-out cross-
#' validation exercise.
#'
#'
#' @inheritParams KnnDistCV
#' @param GroupSize the sample size you wish the function to resample the data to. If NA supplied then the smallest sample size of the groups supplied will be used.
#' @return Returns a list of 2 vectors. The first vector is the indexed location of the new equally size sample
#'
#'
#' @author Ardern Hulme-Beaman



BalancedGrps <- function(GroupMembership, GroupSize=NA){
  #Data=PairwiseShapeDistMat; GroupMembership=chr(Groups[GrpPos])

  if (is.na(GroupSize)){
    minSampSize <- min(table(as.character(GroupMembership)))
  } else {
    minSampSize <- GroupSize
  }

  sampleindex <- 1:length(GroupMembership)
  RandSamp <- stats::aggregate(x = sampleindex, by=list(factor(GroupMembership)), sample, size=minSampSize)
  Index <- c(t(RandSamp$x))
  GroupMem <- c(t(stats::aggregate(RandSamp$Group.1, list(RandSamp$Group.1), rep, minSampSize)$x))
  return(list(IndexedLocations=Index, Newfactors=GroupMem))

}



#' Function for breaking tied votes
#'
#' This function is primarily for use in other k-NN functions. It takes a vector and applies
#' the chosen tiebreaker method. This includes either randomly selecting one of the tied groups,
#' returning all groups collapsed into a single character string, or returning 'UnIDed'.
#'
#' @param X a vector of the tied classifications to be considered.
#' @inheritParams KnnDistCV
#' @return Returns a character string of the identification based on the chosen tiebreaker method.
#'
#'
#' @author Ardern Hulme-Beaman


TieBreakerFunction <- function(X, TieBreaker){
  if (TieBreaker=='Random' && length(TieBreaker)==1){
    ReturnedVote <- sample(X, size = 1)
  } else if (TieBreaker=='Remove' && length(TieBreaker)==1){
    ReturnedVote <- 'UnIDed'
  } else if (TieBreaker=='Report' && length(TieBreaker)==1){
    ReturnedVote <- paste(X, collapse = '_')
  } else {
    warning('Warning: TieBreaker not set or not set to recognisible method. Function will revert to Report for TieBreaker argument. If you have supplied a TieBreaker argument please check capitalisation and spelling.')
    ReturnedVote <- paste(X, collapse = '_')
  }
  return(ReturnedVote)
}




#' Function for calculating the majority vote
#'
#' This function carried out k-NN analyses on the provided data.
#'
#' @param X a vector or matrix of the nearest neighbours. If \code{Weighting=TRUE} then a matrix is expected with the first column being the membership names of nearest neighbours and the second column should be the distances to each neighbour; if \code{Weighting=FALSE} then the just a vector of the membership names of nearest neighbours is required.
#' @inheritParams KnnDistCV
#' @param Weighting is a logical TRUE or FALSE determining whether classification will be based on a weighted K. In this method a simple weighting system is used where weights are 1/distance, where distances is the distance of the considered neighbour from the unknown.
#' @return Returns the classification of the unknown individual to be identified.
#'
#'
#' @author Ardern Hulme-Beaman



KVote <- function(X, K, Weighting=FALSE, TieBreaker=c('Random', 'Remove', 'Report')){
  #VotingData=KIDmat[,1]
  #VotingData=c(rep('East', 5), rep('West', 5))
  #Weighting=Weighted
  #VotingData=cbind(KIDmat[,1], Kweightmat[,1])

  VotingData <- X

  if (Weighting==FALSE){
    if (is.vector(VotingData)){
      MajorityVote <- names(which(table(VotingData[1:K])==max(table(VotingData[1:K]))))
    } else {
      stop('Error: VotingData is not a vector, KVote function expects VotingData to be a vecotr when Weighting=FALSE')
    }

  } else {
    if (dim(VotingData)[2]<2){
      stop('Error: VotingData must be a matrix of 2 columns: the first column must be the nearest neighbour classifications, the second column must be the weighting values')
    }
    WeightingScores <- stats::aggregate(x = 1/chr2nu(VotingData[1:K,2]), by = list(as.factor(VotingData[1:K,1])), FUN = sum)
    MajorityVote <- as.character(WeightingScores$Group.1[which(WeightingScores$x==max(WeightingScores$x))])
  }

  if (length(MajorityVote)>1){
    ReturnedVote <- TieBreakerFunction(MajorityVote, TieBreaker)
  } else {
    ReturnedVote <- MajorityVote
  }

  return(ReturnedVote)
}

