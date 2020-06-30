
#' K-Nearest Neighbour correct cross-validation with distance input
#'
#' This function takes a square matrix of distances among specimens of known group membership
#' and returns the results of a leave-one-out correct cross validation identification for each
#' specimen to provide a correct cross-validation percentage.
#'
#' The function also provides functionality to resample unequal groups to equal sample size a set
#' number of times.
#'
#' This function applies both a weighted approach and an unweighted approach and returns both results.
#'
#'
#'
#' @param DistMat is a square matrix of pairwise distances among all reference specimens.
#' @param GroupMembership a character or factor vector in the same order as the distance data to denote group membership.
#' @param K is the number of nearest neighbours that the method will use for assigning group classification.
#' @param Equal indicates where groups should be sampled to equal sample size
#' @param EqualIter sets the number of iterations resampling to equal sample size will be carried out.
#' @param SampleSize is the sample number that groups will be subsampled to if \code{Equal} is set to TRUE. The default is set to NA and will therefore use the smallest sample size of the groups provided.
#' @param TieBreaker is the method used to break ties if there is no majority resulting from K. Three methods are available('Random', 'Remove' and 'Report'): Random randomly returns one of tied classifications; Remove returns 'UnIDed' for the classification; Report returns a the multiple classifications as a single character string with tied classifications separated by '_'. NOTE: for correct cross-validation proceedures the results of both Report will be considered an incorrect identification even if one of the multiple reported classifications is correct.
#' @param Verbose determines whether the cross-validation results for each reference specimen is returned. Note that if this is set to TRUE and Equal is set to TRUE the funtion will return a list with the results of each iteration which will slow the process dramatically and take a lot of local memory.
#' @param IgnorePrompts if both Verbose and Equal are set to TRUE, then the funciton will ask if you are sure you wish to continue; setting IgnorePrompts to TRUE will ignore this question.
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classificaiton.
#'
#'
#' @keywords shape distances
#' @keywords Geometric morphometric distances
#' @author Ardern Hulme-Beaman
#' @export


KnnDistCV <- function(DistMat, GroupMembership, K, Equal=TRUE, EqualIter=100, SampleSize=NA, TieBreaker=c('Random', 'Remove', 'Report'), Verbose=FALSE, IgnorePrompts=FALSE){
  #K=2
  #DistMat=ProcDTableRes; GroupMembership=Groups; K=10; Equal=TRUE; EqualIter=100
  #Weighted=TRUE; TieBreaker='Report'; Equal=FALSE; EqualIter=100; Verbose=TRUE

  if (dim(DistMat)[1]!=length(GroupMembership)){
    stop('number of specimens in DistMat does not match the number of specimens listed in GroupMembership')
  }

  #Verbose=TRUE; Equal=TRUE
  if (IgnorePrompts==FALSE){
    if (Equal==TRUE && Verbose==TRUE){
      print(' Both Equal and Verbose are set to TRUE.')
      print(' Results of each resampling iteration will be saved and returned as a list.')
      print(' This can dramatically slow the process. ')
      UserInput <- readline(prompt = " Do you wish to continue? Please respond y or n: ")

      while (!UserInput %in% c('y','n')) {
        UserInput <- readline(prompt = " Response was neither n or y. Please respond y or n: ")
      }

      if (UserInput=='n'){
        stop("You have indicated you do not want to continue with Verbose and Equal set to TRUE; Please adjust those arguments and rerun.")
      } else if (UserInput=='y'){
        UserGood <- 'good'
        print("You have indicated you wish to continue, this may take some time.")
      }


    }
  }


  ResultsTable <- NULL

  if (Equal==TRUE){

    IterRes <- list()

    for (iter in 1:EqualIter){
      BalancingGrps <- BalancedGrps(GroupMembership, GroupSize = SampleSize)

      BalancedDistMat <- DistMat[BalancingGrps$IndexedLocations,BalancingGrps$IndexedLocations]

      #this sorts each row of the distance matrix individually into an order of smallest to greatest distance
      #column order is maintained
      AdjSortedDistMat <- SortedDistMat <- BalancedDistMat
      rownames(AdjSortedDistMat) <- rownames(SortedDistMat) <- 1:dim(SortedDistMat)[1]
      for (i in 1:dim(BalancedDistMat)[1]){
        #i <- 2

        AdjSortedDistMat[,i] <- sort(BalancedDistMat[,i])
        SortedDistMat[,i] <- as.character(BalancingGrps$Newfactors[sort(SortedDistMat[,i], index.return=TRUE)$ix])
      }

      #Removing the first row because this will represent the distance 0 i.e. the distance from the column specimen to itself
      KIDmat <- SortedDistMat[-1,]
      Kweightmat <- AdjSortedDistMat[-1,]


      KArray <- array(data = NA, dim = c(dim(KIDmat), 2))
      KArray[,,1] <- KIDmat
      KArray[,,2] <- Kweightmat

      dimnames(KArray) <- list(dimnames(Kweightmat)[[1]], dimnames(Kweightmat)[[2]], c('Call', 'Weight'))

      WeightedRes <- apply(X = KArray, MARGIN = 2, FUN = KVote, K=K, Weighting=TRUE, TieBreaker=TieBreaker)
      UnweightedRes <- apply(X = KArray[,,1], MARGIN = 2, FUN = KVote, K=K, Weighting=FALSE, TieBreaker=TieBreaker)

      ResTable <- cbind(rownames(BalancedDistMat), BalancingGrps$Newfactors, WeightedRes, UnweightedRes)
      colnames(ResTable) <- c('ID', 'True.Classification', 'Weighted.Classification', 'Unweighted.Classification')

      WeightedCCVPercent <- sum(BalancingGrps$Newfactors==WeightedRes)/length(BalancingGrps$Newfactors)
      UnweightedCCVPercent <- sum(BalancingGrps$Newfactors==UnweightedRes)/length(BalancingGrps$Newfactors)

      ResultsTable <- rbind(ResultsTable, c(WeightedCCVPercent, UnweightedCCVPercent))


      if (Verbose==TRUE){
        IterRes[[iter]] <- ResTable
      }

    }

    colnames(ResultsTable) <- c('Weighted', 'Unweighted')

    if (Verbose==TRUE){
      return(list(CCV.Percentages=ResultsTable, Individual.Iteration.Results=IterRes))
    } else {
      return(ResultsTable)
    }

  } else {

    #this sorts each row of the distance matrix individually into an order of smallest to greatest distance
    #column order is maintained
    AdjSortedDistMat <- SortedDistMat <- DistMat
    rownames(AdjSortedDistMat) <- rownames(SortedDistMat) <- 1:dim(SortedDistMat)[1]
    for (i in 1:dim(DistMat)[1]){
      #i <- 2

      AdjSortedDistMat[,i] <- sort(DistMat[,i])
      SortedDistMat[,i] <- as.character(GroupMembership[sort(SortedDistMat[,i], index.return=TRUE)$ix])
    }

    #Removing the first row bbecause this will represent the distance 0 i.e. the distance from the column specimen to itself
    KIDmat <- SortedDistMat[-1,]
    Kweightmat <- AdjSortedDistMat[-1,]


    KArray <- array(data = NA, dim = c(dim(KIDmat), 2))
    KArray[,,1] <- KIDmat
    KArray[,,2] <- Kweightmat

    dimnames(KArray) <- list(dimnames(Kweightmat)[[1]], dimnames(Kweightmat)[[2]], c('Call', 'Weight'))

    WeightedRes <- apply(X = KArray, MARGIN = 2, FUN = KVote, K=K, Weighting=TRUE, TieBreaker=TieBreaker)
    UnweightedRes <- apply(X = KArray[,,1], MARGIN = 2, FUN = KVote, K=K, Weighting=FALSE, TieBreaker=TieBreaker)

    ResTable <- cbind(rownames(DistMat), GroupMembership, WeightedRes, UnweightedRes)
    colnames(ResTable) <- c('ID', 'True.Classification', 'Weighted.Classification', 'Unweighted.Classification')

    WeightedCCVPercent <- sum(GroupMembership==WeightedRes)/length(GroupMembership)
    UnweightedCCVPercent <- sum(GroupMembership==UnweightedRes)/length(GroupMembership)

    ResultsTable <- rbind(ResultsTable, c(WeightedCCVPercent, UnweightedCCVPercent))

    colnames(ResultsTable) <- c('Weighted', 'Unweighted')

    if (Verbose==TRUE){
      return(list(CCV.Percentages=ResultsTable, Individual.Classifications=as.data.frame(ResTable)))
    } else {
      return(list(CCV.Percentages=ResultsTable))
    }
  }
}


#' Stepwise K-Nearest Neighbour correct cross-validation with distance input
#'
#' This function takes a square matrix of distances among specimens of known group membership
#' and returns the results of a leave-one-out correct cross validation identification exercise for
#' each incremental increase in k. The results of the analyses can be plotted to visualise the
#' change in correct identification given changes in k.
#'
#' The function also provides functionality to resample unequal groups to equal sample size a set
#' number of times.
#'
#' This function applies both a weighted approach and an unweighted appraoch and returns both results.
#'
#'
#'
#' @param Kmax This sets the maximum K that K will increase to stepwise.
#' @param PrintProg Only used when resampling to equal sample size is used (i.e. \code{Equal = TRUE}) and shows what iteration the function is on. Set to FALSE to ignore.
#' @param PlotResults logical when set to TRUE the results are plotted. When \code{Equal = TRUE} a polygon is plotted marking the 5th and 9th percentile.
#' @inheritParams KnnDistCV
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classification for both weighted and unweighted approaches.
#' @details When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
#'
#' @keywords shape distances
#' @keywords Geometric morphometric distances
#' @author Ardern Hulme-Beaman
#' @import svMisc
#' @export


KnnDistCVStepwise <- function(DistMat, GroupMembership, Kmax, Equal=TRUE, EqualIter=100, SampleSize=NA, TieBreaker=c('Random', 'Remove', 'Report'), Verbose=FALSE, PrintProg=TRUE, PlotResults=TRUE){
  #DistMat = RatDistMat;  GroupMembership = RatData$Info$Species;  Kmax = 10;  PrintProg = FALSE; Equal = FALSE;  PlotResults = TRUE; TieBreaker = 'Remove'



  MinSamp <- min(table(as.character(GroupMembership)))

  if (Kmax>MinSamp){
    warning('Kmax is set to higher than the smallest sample size.')
  }

  if (Equal==FALSE){
    EqualIter <- 1
  }

  ResultsTable <- list(Unweighted.Results=matrix(NA, nrow = EqualIter, ncol = Kmax), Weighted.Results=matrix(NA, nrow = EqualIter, ncol = Kmax))

  if (Equal==TRUE){

    for (iter in 1:EqualIter){
      BalancingGrps <- BalancedGrps(GroupMembership, SampleSize)

      BalancedDistMat <- DistMat[BalancingGrps$IndexedLocations,BalancingGrps$IndexedLocations]

      #this sorts each row of the distance matrix individually into an order of smallest to greatest distance
      #column order is maintained
      AdjSortedDistMat <- SortedDistMat <- BalancedDistMat
      rownames(AdjSortedDistMat) <- rownames(SortedDistMat) <- 1:dim(SortedDistMat)[1]


      for (i in 1:dim(BalancedDistMat)[1]){
        AdjSortedDistMat[,i] <- sort(BalancedDistMat[,i])
        SortedDistMat[,i] <- as.character(BalancingGrps$Newfactors[sort(SortedDistMat[,i], index.return=TRUE)$ix])
      }

      #Removing the first row bbecause this will represent the distance 0 i.e. the distance from the column specimen to itself

      KArray <- array(data = NA, dim = c(dim(SortedDistMat[-1,]), 2))
      KArray[,,1] <- as.matrix(SortedDistMat[-1,])
      KArray[,,2] <- as.matrix(AdjSortedDistMat[-1,])

      dimnames(KArray) <- list(dimnames(AdjSortedDistMat[-1,])[[1]], dimnames(AdjSortedDistMat[-1,])[[2]], c('Call', 'Weight'))


      for (K in 1:Kmax){
        #K=1
        WeightedRes <- apply(X = KArray, MARGIN = 2, FUN = KVote, K=K, Weighting=TRUE, TieBreaker=TieBreaker)
        UnweightedRes <- apply(X = KArray[,,1], MARGIN = 2, FUN = KVote, K=K, Weighting=FALSE, TieBreaker=TieBreaker)

        WeightedCCVPercent <- sum(BalancingGrps$Newfactors==WeightedRes)/length(BalancingGrps$Newfactors)
        UnweightedCCVPercent <- sum(BalancingGrps$Newfactors==UnweightedRes)/length(BalancingGrps$Newfactors)

        ResultsTable$Unweighted.Results[iter, K] <- UnweightedCCVPercent
        ResultsTable$Weighted.Results[iter, K] <- WeightedCCVPercent
      }

      if (PrintProg==TRUE){
        svMisc::progress(value = iter, max.value = EqualIter, progress.bar = FALSE)
        #Sys.sleep(0.001)

      }

    }

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

    if (Verbose==TRUE){
      return(ResultsTable)
    } else {
      return(list(Unweighted=colMeans(ResultsTable$Unweighted.Results), Weighted=colMeans(ResultsTable$Weighted.Results)))
    }


  } else {

    #this sorts each row of the distance matrix individually into an order of smallest to greatest distance
    #column order is maintained
    AdjSortedDistMat <- SortedDistMat <- DistMat
    rownames(AdjSortedDistMat) <- rownames(SortedDistMat) <- 1:dim(SortedDistMat)[1]
    for (i in 1:dim(DistMat)[1]){
      #i <- 2

      AdjSortedDistMat[,i] <- sort(DistMat[,i])
      SortedDistMat[,i] <- as.character(GroupMembership[sort(SortedDistMat[,i], index.return=TRUE)$ix])
    }

    #Removing the first row bbecause this will represent the distance 0 i.e. the distance from the column specimen to itself
    KArray <- array(data = NA, dim = c(dim(SortedDistMat[-1,]), 2))
    KArray[,,1] <- as.matrix(SortedDistMat[-1,])
    KArray[,,2] <- as.matrix(AdjSortedDistMat[-1,])

    dimnames(KArray) <- list(dimnames(AdjSortedDistMat[-1,])[[1]], dimnames(AdjSortedDistMat[-1,])[[2]], c('Call', 'Weight'))


    for (K in 1:Kmax){
      #K=1

      WeightedRes <- apply(X = KArray, MARGIN = 2, FUN = KVote, K=K, Weighting=TRUE, TieBreaker=TieBreaker)
      UnweightedRes <- apply(X = KArray[,,1], MARGIN = 2, FUN = KVote, K=K, Weighting=FALSE, TieBreaker=TieBreaker)

      WeightedCCVPercent <- sum(GroupMembership==WeightedRes)/length(GroupMembership)
      UnweightedCCVPercent <- sum(GroupMembership==UnweightedRes)/length(GroupMembership)

      ResultsTable$Unweighted.Results[1, K] <- UnweightedCCVPercent
      ResultsTable$Weighted.Results[1, K] <- WeightedCCVPercent
    }

    if (PlotResults==TRUE){

      graphics::plot(y = colMeans(ResultsTable$Unweighted.Results), x = 1:Kmax, type = 'n', ylim = c(10,105), ylab = 'CCV %', xlab = 'K')
      graphics::abline(h = seq(from = 20, to = 100, by = 10), v = seq(from = 2, to = Kmax, by =2), lty = '1919')


      WeightRange <- apply(ResultsTable$Weighted.Results, MARGIN = 2, FUN = stats::quantile, probs = c(.05, .95))
      graphics::polygon(x = c(1:Kmax, Kmax:1),
                        y = c(WeightRange[1,], WeightRange[2,Kmax:1])*100,
                        col = transpar('darkblue', alpha = 75),
                        border = NA)

      graphics::lines(y = colMeans(ResultsTable$Weighted.Results*100), x = 1:Kmax, col='darkblue', lwd=3)

      UnweightRange <- apply(ResultsTable$Unweighted.Results, MARGIN = 2, FUN = stats::quantile, probs = c(.05, .95))
      graphics::polygon(x = c(1:Kmax, Kmax:1),
                        y = c(UnweightRange[1,], UnweightRange[2,Kmax:1])*100,
                        col = transpar('lightblue', alpha = 75),
                        border = NA)
      graphics::lines(y = colMeans(ResultsTable$Unweighted.Results*100), x = 1:Kmax, col='lightblue', lwd=3)

      graphics::legend('bottomright', legend = c('Weighted', 'Unweighted'), col = c('darkblue', 'lightblue'), lty=1, lwd=3, bty = 'o')

      }

      return(list(CCV.Percentages=ResultsTable))
  }

}



