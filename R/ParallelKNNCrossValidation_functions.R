

#' K-Nearest Neighbour correct cross-validation with distance input using parallel processing
#'
#' This function description NEEDS COMPLETION
#'
#'
#'
#' @inheritParams KnnDistCV
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classificaiton.
#'
#'
#' @author Ardern Hulme-Beaman
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export


KnnDistCVPar <- function(DistMat, GroupMembership, K, EqualIter=100, SampleSize=NA, TieBreaker=c('Random', 'Remove', 'Report'), Verbose=FALSE){
  #K=2
  #DistMat=ProcDTableRes; GroupMembership=Groups; K=10; Equal=TRUE; EqualIter=100
  #Weighted=TRUE; TieBreaker='Report'



  chr2nu <- function(X){
    as.numeric(as.character(X))
  }

  ParOutput <- function(PreviousResults, ResultList){

    NewResults <- PreviousResults
    for (i in 1:length(ResultList)){
      NewResults[[i]] <- rbind(PreviousResults[[i]], ResultList[[i]])
    }

    return(NewResults)
  }


  BalancedGrps <- function(GroupMembership, GroupSize=SampleSize){
    #Data=PairwiseShapeDistMat; GroupMembership=chr(Groups[GrpPos])

    if (is.na(GroupSize)){
      minSampSize <- min(table(GroupMembership))
    } else {
      minSampSize <- GroupSize
    }

    sampleindex <- 1:length(GroupMembership)
    RandSamp <- stats::aggregate(x = sampleindex, by=list(factor(GroupMembership)), sample, size=minSampSize)
    Index <- c(t(RandSamp$x))
    GroupMem <- c(t(stats::aggregate(RandSamp$Group.1, list(RandSamp$Group.1), rep, minSampSize)$x))
    return(list(IndexedLocations=Index, Newfactors=GroupMem))

  }


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
      if (TieBreaker=='Random' && length(TieBreaker)==1){
        ReturnedVote <- sample(MajorityVote, size = 1)
      } else if (TieBreaker=='Remove' && length(TieBreaker)==1){
        ReturnedVote <- 'UnIDed'
      } else if (TieBreaker=='Report' && length(TieBreaker)==1){
        ReturnedVote <- paste(MajorityVote, collapse = '_')
      } else {
        warning('Warning: TieBreaker not set or not set to recognisible method. Function will revert to Report for TieBreaker argument. If you have supplied a TieBreaker argument please check capitalisation and spelling.')
        ReturnedVote <- paste(MajorityVote, collapse = '_')
      }
    } else {
      ReturnedVote <- MajorityVote
    }

    return(ReturnedVote)
  }

  ParEqualIter <- function(DistData, GrpMem, ParK=K, ParTieBreaker=TieBreaker, ParVerbose=Verbose){
    #DistData=DistMat; GrpMem=GroupMembership; ParK=K; ParTieBreaker='Report'; ParVerbose=FALSE

    BalancingGrps <- BalancedGrps(GrpMem, SampleSize)

    BalancedDistMat <- DistData[BalancingGrps$IndexedLocations,BalancingGrps$IndexedLocations]

    #this sorts each row of the distance matrix individually into an order of smallest to greatest distance
    #column order is maintained
    AdjSortedDistMat <- SortedDistMat <- BalancedDistMat
    rownames(AdjSortedDistMat) <- rownames(SortedDistMat) <- 1:dim(SortedDistMat)[1]

    for (i in 1:dim(BalancedDistMat)[1]){
      #i <- 2

      AdjSortedDistMat[,i] <- sort(BalancedDistMat[,i])
      SortedDistMat[,i] <- as.character(BalancingGrps$Newfactors[sort(SortedDistMat[,i], index.return=TRUE)$ix])
    }

    #Removing the first row bbecause this will represent the distance 0 i.e. the distance from the column specimen to itself
    KIDmat <- SortedDistMat[-1,]
    Kweightmat <- AdjSortedDistMat[-1,]

    KArray <- array(data = NA, dim = c(dim(KIDmat), 2))
    KArray[,,1] <- KIDmat
    KArray[,,2] <- Kweightmat

    dimnames(KArray) <- list(dimnames(Kweightmat)[[1]], dimnames(Kweightmat)[[2]], c('Call', 'Weight'))

    WeightedRes <- apply(X = KArray, MARGIN = 2, FUN = KVote, K=ParK, Weighting=TRUE, TieBreaker=ParTieBreaker)
    UnweightedRes <- apply(X = KIDmat, MARGIN = 2, FUN = KVote, K=ParK, Weighting=FALSE, TieBreaker=ParTieBreaker)

    WeightedCCVPercent <- sum(BalancingGrps$Newfactors==WeightedRes)/length(BalancingGrps$Newfactors)
    UnweightedCCVPercent <- sum(BalancingGrps$Newfactors==UnweightedRes)/length(BalancingGrps$Newfactors)

    ResultsSummaryTable <- c(WeightedCCVPercent, UnweightedCCVPercent)

    if (ParVerbose==TRUE){
      ResTable <- cbind(rownames(BalancedDistMat), BalancingGrps$Newfactors, WeightedRes, UnweightedRes)
      colnames(ResTable) <- c('ID', 'True.Classification', 'Weighted.Classification', 'Unweighted.Classification')
      return(list(ResultsSummaryTable, ResTable))
    } else {
      return(list(ResultsSummaryTable, NA))
    }
  }


  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  a <- 1
  ParResults <- foreach::foreach(a = 1:EqualIter, .combine = ParOutput) %dopar% {
    ParEqualIter(DistData = DistMat, GrpMem = GroupMembership, ParK = K, ParTieBreaker = 'Report', ParVerbose = Verbose)
  }


  parallel::stopCluster(clust)

  colnames(ParResults[[1]]) <- c('Weighted', 'Unweighted')

  if (Verbose==TRUE){
    names(ParResults) <- c('Iteration.Summaries', 'Verbose.Output')
    return(ParResults)
  } else {
    return(ParResults[[1]])
  }

}


#' Stepwise K-Nearest Neighbour correct cross-validation with distance input using parallel processing
#'
#' This function NEEDS COMPLETION
#'
#'
#' @inheritParams KnnDistCVStepwise
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classificaiton.
#' @details When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
#'
#' @section Citations:
#'
#' Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
#' https://CRAN.R-project.org/package=shapes
#'
#'
#' @keywords shape distances
#' @keywords Geometric morphometric distances
#' @author Ardern Hulme-Beaman
#' @import shapes
#' @import svMisc
#' @export


KnnDistCVStepwisePar <- function(DistMat, GroupMembership, Kmax, EqualIter=100, SampleSize=NA, TieBreaker=c('Random', 'Remove', 'Report'), PlotResults=TRUE){

  #Kmax=20
  #DistMat=ProcDTableRes; GroupMembership=Groups; K=10; Equal=TRUE; EqualIter=100
  #Weighted=TRUE; TieBreaker='Report'
  #PrintProg=TRUE
  #Verbose=TRUE; Equal=TRUE

  MinSamp <- min(table(as.character(GroupMembership)))

  if (Kmax>MinSamp){
    warning('Kmax is set to higher than the smallest sample size.')
  }



  chr2nu <- function(X){
    as.numeric(as.character(X))
  }

  ParOutput <- function(PreviousResults, ResultList){

    NewResults <- PreviousResults
    for (i in 1:length(ResultList)){
      NewResults[[i]] <- rbind(PreviousResults[[i]], ResultList[[i]])
    }

    return(NewResults)
  }


  BalancedGrps <- function(GroupMembership, GroupSize=SampleSize){
    #Data=PairwiseShapeDistMat; GroupMembership=chr(Groups[GrpPos])

    if (is.na(GroupSize)){
      minSampSize <- min(table(GroupMembership))
    } else {
      minSampSize <- GroupSize
    }

    sampleindex <- 1:length(GroupMembership)
    RandSamp <- stats::aggregate(x = sampleindex, by=list(factor(GroupMembership)), sample, size=minSampSize)
    Index <- c(t(RandSamp$x))
    GroupMem <- c(t(stats::aggregate(RandSamp$Group.1, list(RandSamp$Group.1), rep, minSampSize)$x))
    return(list(IndexedLocations=Index, Newfactors=GroupMem))

  }


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
      if (TieBreaker=='Random' && length(TieBreaker)==1){
        ReturnedVote <- sample(MajorityVote, size = 1)
      } else if (TieBreaker=='Remove' && length(TieBreaker)==1){
        ReturnedVote <- 'UnIDed'
      } else if (TieBreaker=='Report' && length(TieBreaker)==1){
        ReturnedVote <- paste(MajorityVote, collapse = '_')
      } else {
        warning('Warning: TieBreaker not set or not set to recognisible method. Function will revert to Report for TieBreaker argument. If you have supplied a TieBreaker argument please check capitalisation and spelling.')
        ReturnedVote <- paste(MajorityVote, collapse = '_')
      }
    } else {
      ReturnedVote <- MajorityVote
    }

    return(ReturnedVote)
  }

  ParEqualIterStepwise <- function(DistData, GrpMem, ParKmax=Kmax, ParTieBreaker=TieBreaker, ParSampleSize=SampleSize){
    #DistData=DistMat; GrpMem=GroupMembership; ParKmax=K; ParTieBreaker='Report'; ParVerbose=FALSE; ParSampleSize=NA

    BalancingGrps <- BalancedGrps(GrpMem, ParSampleSize)

    BalancedDistMat <- DistData[BalancingGrps$IndexedLocations,BalancingGrps$IndexedLocations]

    #this sorts each row of the distance matrix individually into an order of smallest to greatest distance
    #column order is maintained
    AdjSortedDistMat <- SortedDistMat <- BalancedDistMat
    rownames(AdjSortedDistMat) <- rownames(SortedDistMat) <- 1:dim(SortedDistMat)[1]

    for (i in 1:dim(BalancedDistMat)[1]){
      #i <- 2

      AdjSortedDistMat[,i] <- sort(BalancedDistMat[,i])
      SortedDistMat[,i] <- as.character(BalancingGrps$Newfactors[sort(SortedDistMat[,i], index.return=TRUE)$ix])
    }

    #Removing the first row bbecause this will represent the distance 0 i.e. the distance from the column specimen to itself
    KArray <- array(data = NA, dim = c(dim(SortedDistMat[-1,]), 2))
    KArray[,,1] <- as.matrix(SortedDistMat[-1,])
    KArray[,,2] <- as.matrix(AdjSortedDistMat[-1,])

    dimnames(KArray) <- list(dimnames(AdjSortedDistMat[-1,])[[1]], dimnames(AdjSortedDistMat[-1,])[[2]], c('Call', 'Weight'))

    ResultsTable <- list(Unweighted.Results=matrix(NA, nrow = 1, ncol = ParKmax), Weighted.Results=matrix(NA, nrow = 1, ncol = ParKmax))

    for (K in 1:ParKmax){

      WeightedRes <- apply(X = KArray, MARGIN = 2, FUN = KVote, K=K, Weighting=TRUE, TieBreaker=TieBreaker)
      UnweightedRes <- apply(X = KArray[,,1], MARGIN = 2, FUN = KVote, K=K, Weighting=FALSE, TieBreaker=TieBreaker)

      WeightedCCVPercent <- sum(BalancingGrps$Newfactors==WeightedRes)/length(BalancingGrps$Newfactors)
      UnweightedCCVPercent <- sum(BalancingGrps$Newfactors==UnweightedRes)/length(BalancingGrps$Newfactors)

      ResultsTable$Unweighted.Results[1, K] <- UnweightedCCVPercent
      ResultsTable$Weighted.Results[1, K] <- WeightedCCVPercent
    }


    return(ResultsTable)
  }


  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  a <- 1
  ParResults <- foreach::foreach(a = 1:EqualIter, .combine = ParOutput) %dopar% {
    ParEqualIterStepwise(DistData = DistMat, GrpMem = GroupMembership, ParKmax = Kmax, ParTieBreaker = 'Report', ParSampleSize = NA)
  }


  parallel::stopCluster(clust)


  ResultsTable <- ParResults

  if (PlotResults==TRUE){

    graphics::plot(y = colMeans(ResultsTable$Unweighted.Results), x = 1:Kmax, type = 'n', ylim = c(10,105), ylab = 'CCV %', xlab = 'K')
    graphics::abline(h = seq(from = 20, to = 100, by = 10), v = seq(from = 2, to = Kmax, by =2), lty = '1919')

    WeightRange <- apply(ResultsTable$Weighted.Results, MARGIN = 2, FUN = stats::quantile, probs = c(.05, .95))
    graphics::polygon(x = c(1:Kmax, Kmax:1),
            y = c(WeightRange[1,], WeightRange[2,Kmax:1])*100,
            col = transpar('darkblue', alpha = 25),
            border = NA)

    graphics::lines(y = colMeans(ResultsTable$Weighted.Results*100), x = 1:Kmax, col='darkblue', lwd=3)

    UnweightRange <- apply(ResultsTable$Unweighted.Results, MARGIN = 2, FUN = stats::quantile, probs = c(.05, .95))
    graphics::polygon(x = c(1:Kmax, Kmax:1),
            y = c(UnweightRange[1,], UnweightRange[2,Kmax:1])*100,
            col = transpar('lightblue', alpha = 95),
            border = NA)
    graphics::lines(y = colMeans(ResultsTable$Unweighted.Results*100), x = 1:Kmax, col='lightblue', lwd=3)

    graphics::legend('bottomright', legend = c('Weighted', 'Unweighted'), col = c('darkblue', 'lightblue'), lty=1, lwd=3, bty = 'o')
  }

  return(ResultsTable)
}


