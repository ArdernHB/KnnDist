
#' K-Nearest Neighbour single specimen identification
#'
#' This function takes a vector of distances to an unknown specimen that is to be identified and returns
#' an identification for the specimen. This function applies both a weighted approach and an unweighted
#' appraoch and returns both results.
#'
#' @param X a numeric vector of distances from the unknown specimen to the all the reference specimens.
#' @inheritParams KnnDistCV
#' @details The function calculates the classification based on a majority vote from k nearest neighbours.
#' The method also calculates the classification using a weighted approach. In the weighted approach each
#' of the k nearest neighbours are given a weight calculated as 1/distance. The weights are summed by group
#' for the k nearest neighbours and the unknown specimen is assigned to the group with the highest summed weight.
#'
#' @return Returns a list of two objects, the classification result from a weighted approach and the classification result from an unweighted approach.
#'
#'
#' @author Ardern Hulme-Beaman
#' @export


KnnIDingSingleInd <- function(X, K, GroupMembership, TieBreaker = c('Random', 'Remove', 'Report')){
  #Dist2Unknown=ProcDTableRes[1,-1]
  #Weighted=TRUE
  #Membership=Groups

  Dist2Unknown <- X
  VotingWeights <- sort(chr2nu(Dist2Unknown), index.return=TRUE)
  VotingData <- GroupMembership[VotingWeights$ix]

  UnweightedMajorityVote <- names(which(table(VotingData[1:K])==max(table(VotingData[1:K]))))

  WeightingScores <- stats::aggregate(x = 1/chr2nu(VotingWeights$x[1:K]), by = list(as.factor(VotingData[1:K])), FUN = sum)
  WeightedMajorityVote <- as.character(WeightingScores$Group.1[which(WeightingScores$x==max(WeightingScores$x))])

  if (length(WeightedMajorityVote)>1){
    ReturnedWeightedVote <- TieBreakerFunction(WeightedMajorityVote, TieBreaker)
  } else {
    ReturnedWeightedVote <- WeightedMajorityVote
  }

  if (length(UnweightedMajorityVote)>1){
    ReturnedUnweightedVote <- TieBreakerFunction(UnweightedMajorityVote, TieBreaker)
  } else {
    ReturnedUnweightedVote <- UnweightedMajorityVote
  }

  return(list(Weighted.Result=ReturnedWeightedVote, Unweighted.Result=ReturnedUnweightedVote))

}




#' k-Nearest Neighbour multiple specimen identification
#'
#' This function is for an unbalanced kNN identification design applied to multiple unknown specimens.
#'
#' @param UnknownIdentifier the name used in the \code{GroupMembership} argument to denote specimens to be identified. Only one name can be supplied as Unknown; default is set to 'Unknown'.
#' @inheritParams KnnDistCV
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classificaiton.
#'
#'
#' @author Ardern Hulme-Beaman
#' @export



KnnDistIDingGroup <- function(DistMat, GroupMembership, UnknownIdentifier = 'Unknown', K, TieBreaker){
  #Rpraetor$Lat.Long
  #DistMat=ProcDTableRes; GroupMembership=UnkGroups; UnknownIdentifier = 'Unknown'; UnknownSpecimenIDs=rownames(Rpraetor$Lat.Long); K=1; Weighted=TRUE; TieBreaker = 'Random'
  UnkPos <- which(GroupMembership%in%UnknownIdentifier)
  UnkDists <- DistMat[UnkPos,-UnkPos]
  colnames(UnkDists) <- GroupMembership[-UnkPos]

  RawIDRes <- apply(X = UnkDists, MARGIN = 1, FUN = KnnIDingSingleInd, K=K, Membership=colnames(UnkDists), TieBreaker = TieBreaker)

  Results <- do.call('rbind', RawIDRes)

  return(Results)

}


#' K-Nearest Neighbour multiple specimen identification with resampling to equal sample size
#'
#' This function is for a balanced KNN identification design applied to multiple unknown specimens.
#' Groups are resampled iteratively to equal sample size; note this is done by downsampleing the groups
#' to a sample size set by the user or if left to default to the sample size of the smallest groups.
#' Bootstrap resampling is not provided as an option because this can be problematic for nearest neighbour
#' approaches because of duplication of neighbours.
#'
#'
#' @inheritParams KnnDistIDingGroup
#' @inheritParams KnnDistCV
#' @param SpecimenIDs should be a list of the specimens unique identifiers for all specimens in the \code{DistMat} object and ensuring they are in the same order as the \code{DistMat} object.
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classificaiton.
#'
#'
#' @author Ardern Hulme-Beaman
#' @export



KnnDistIDingBal <- function(DistMat, GroupMembership, UnknownIdentifier = 'Unknown', SpecimenIDs, K, EqualIter=100, TieBreaker, SampleSize=NA){

  #DistMat=PairwiseShapeDistMat; GroupMembership=Groups[GrpPos]; K=10; EqualIter=1000
  #DistMat=SalmonProcDist; GroupMembership= CompDatasort$info$species; UnknownIdentifier="Nunalleq"; UnknownSpecimenIDs=CompDatasort$info$ID[ArchPos]; K=1; EqualIter=1000
  #DistMat = SalmonProcDist; GroupMembership = CompDatasort$info$species; K = 10; equal = TRUE; EqualIter = 100
  #UnknownIdentifier="Nunalleq"
  #UnknownSpecimenIDs <- CompDatasort$info$ID[UnkPos]


  UnkPos <- which(GroupMembership%in%UnknownIdentifier)

  UnkDists <- DistMat[UnkPos,-UnkPos]

  KResultsIterUnweighted <- KResultsIterWeighted <- array(0, dim = c(EqualIter, length(unique(GroupMembership[-UnkPos])), length(UnkPos)), dimnames = list(1:EqualIter, sort(unique(GroupMembership[-UnkPos])), rownames(UnkDists)))

  for (Eq in 1:EqualIter){
    #Eq <- 1

    Balancing <- BalancedGrps(GroupMembership = as.character(GroupMembership[-UnkPos]), GroupSize = SampleSize)

    BalGroupMembership <- c(Balancing$Newfactors, as.character(GroupMembership[UnkPos]))
    BalGroupIndex <- c(Balancing$IndexedLocations, UnkPos)

    IterRes <- KnnDistIDingGroup(DistMat = DistMat[BalGroupIndex,BalGroupIndex], GroupMembership=BalGroupMembership, UnknownIdentifier, K=K, TieBreaker = TieBreaker)

    for (i in 1:length(UnkPos)){
      #i <- 1
      #IterRes[i,1] <- 'Unk'
      if (IterRes[i,1]%in%dimnames(KResultsIterUnweighted)[[2]]){
        KResultsIterWeighted[Eq,which(names(KResultsIterWeighted[Eq,,i])%in%IterRes[i,1]),i] <- 1
        KResultsIterUnweighted[Eq,which(names(KResultsIterUnweighted[Eq,,i])%in%IterRes[i,2]),i] <- 1
      }
    }

  }

  apply2 <- function(X, mar, fun){
    apply(X, MARGIN=mar, FUN=fun)
  }

  Results <- list(Unweighted.IDs = t(apply(KResultsIterUnweighted, MARGIN = 3, FUN = apply2, mar=2, fun=sum))/EqualIter,
                  Weighted.IDs = t(apply(KResultsIterWeighted, MARGIN = 3, FUN = apply2, mar=2, fun=sum))/EqualIter)


  return(Results)
}


