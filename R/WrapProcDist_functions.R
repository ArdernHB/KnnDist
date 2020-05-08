

#' Returns vector of Procrustes distances
#'
#' This function is a wrapper for the \code{\link[shapes]{procdist}} function from the \code{shapes} package.
#' It constructs distances between a target shape and reference shapes using the method specified.
#' @param TargetShape either a matrix where rows are landmarks and columns are dimensions (i.e. x, y, and z) or a vector of coordinated in the format X1, Y1, X2, Y2... etc.
#' @param RefShapes either an array where rows are landmarks and columns are dimensions (i.e. x, y, and z) and slices are specimens or a matrix of coordinated where rows are specimens and columns are landmark coordinates in the format X1, Y1, X2, Y2... etc.
#' @param ShapeArray either TRUE to denote that \code{TargetShape} is in matrix format and \code{Refshape} is an array, or FALSE to denote that \code{TargetShape} is a vector and \code{RefShape} is a matrix. Please do not use a combination of formats as the function does not handle this.
#' @param LMDim set to either 2 or 3 to denote whether the data is 2D or 3D. Default is set to 2.
#' @param Method this is passed to the \code{\link[shapes]{procdist}} function of the \code{shapes} package and must be set to one of the following options: "full" for full Procrustes distance, "partial" for partial Procrustes distance, "Riemannian" for Riemannian shape distance, or "sizeandshape" for size-and-shape Riemannian/Procrustes distance.
#' @return The distance between the point of interest \code{TargetShape} to all other points in the reference material \code{RefShapes}.
#' @section Citations:
#'
#' Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
#' https://CRAN.R-project.org/package=shapes
#'
#'
#' @keywords internal
#' @keywords shape distances
#' @keywords Geometric morphometric distances
#' @author Ardern Hulme-Beaman
#' @import shapes
#' @export

ShapeDist2Specimen <- function(RefShapes, TargetShape, LMDim=2, ShapeArray=TRUE, Method=c('full', 'partial', 'Riemannian', 'sizeandshape')){
  #RefShapes=Rpraetor$LMs[,,-1]; TargetShape=Rpraetor$LMs[,,1]; LMDim=2; ShapeArray=TRUE; Method='full'

  if (ShapeArray==FALSE){
    RefArray <- Mat2Array(RefShapes, LMdim = LMDim)
    TargetArray <- matrix(TargetShape, nrow = length(TargetShape)/LMDim, ncol = LMDim, byrow = TRUE)
  } else {
    RefArray <- RefShapes
    TargetArray <- TargetShape
  }

  ShapeDist <- apply(X = RefArray, MARGIN = 3, FUN = shapes::procdist, y = TargetArray, type=Method)

  return(ShapeDist)
}




#' Wrapper for procdist function to output a distance table
#'
#' This function builds a square matrix of pairwise procrustes distances among specimens using
#' the \code{\link[shapes]{procdist}} function from the \code{shapes} package.
#' @param RefIDs is a vector of the unique identifiers for each of the reference specimens in the reference dataset. These values will be used for naming the columns and rows. Default is set to NA and if a vector is not supplied the columns and rows will remain numbered consecutively.
#' @inheritParams ShapeDist2Specimen
#' @return This function returns a square matrix of Procrustes distances, which is required for both the \code{IDbyDistanceRawDataCCV} and the \code{BoundaryFinder} functions.
#' @section Citations:
#'
#' Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
#' https://CRAN.R-project.org/package=shapes
#'
#' @seealso
#' \code{\link{ProcDistanceTablePar}}
#'
#' @examples
#' #RatDistMat <- ProcDistanceTable(Rpraetor$LMs)
#' @import shapes
#' @export



ProcDistanceTable <- function(RefShapes, LMDim=2, ShapeArray=TRUE, Method=c('full', 'partial', 'Riemannian', 'sizeandshape'), RefIDs=NA){
  #RefIDs <- rownames(Rpraetor$Lat.Long); RefShapes=Rpraetor$LMs
  #ShapeArray=TRUE; Method = 'full'

  if (ShapeArray==FALSE){
    RefArray <- Mat2Array(RefShapes, LMdim = LMDim)
  } else {
    RefArray <- RefShapes
  }

  if (is.na(RefIDs) && length(RefIDs)==1){
    RefIDs <- 1:dim(RefArray)[3]
  } else if (length(RefIDs)<dim(RefArray)[3]){
    stop('Error: RefIDs is not the same length as the number of specimens in RefShapes. Please check your data and make sure the RefIDs correspond with and are in the same order as the specimens in RefShapes')
  }

  if (length(Method)>1 || !Method%in%c('full', 'partial', 'Riemannian', 'sizeandshape')){
    stop('Error: Method argument has not been populated correctly. Please populate it with just one of the following (also make sure not to use capitalisation): full, partial, Riemannian, or sizeandshape')
  }

  ProcDtable <- matrix(0, dim(RefArray)[3], dim(RefArray)[3])
  for (i in 1:dim(RefArray)[3]){
    for (j in 1:dim(RefArray)[3]){
      if (j<i){
        #i <- 1
        #j <- 2
        ProcDtable[i,j] <- shapes::procdist(RefArray[,,i], RefArray[,,j], type=Method)
      }
    }
  }


  ProcDTableRes <- as.matrix(stats::as.dist(ProcDtable))

  colnames(ProcDTableRes) <- rownames(ProcDTableRes) <- RefIDs
  return(ProcDTableRes)
}






#' Wrapper for procdist function to output a distance table using parallel processing
#'
#' This function builds a square matrix of pairwise procrustes distances among specimens using
#' the \code{\link[shapes]{procdist}} function from the \code{shapes} package. This has been set
#' up to run parallel cores using the \code{foreach}, \code{doParallel} and \code{parallel} packages.
#' @inheritParams ShapeDist2Specimen
#' @inheritParams ProcDistanceTable
#' @return This function returns a square matrix of Procrustes distances.
#' @section Citations:
#'
#' Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
#' https://CRAN.R-project.org/package=shapes
#'
#' @examples
#' #See ProcDistanceTable() function example
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import shapes
#' @export



ProcDistanceTablePar <- function(RefShapes, LMDim=2, ShapeArray=TRUE, Method = c('full', 'partial', 'Riemannian', 'sizeandshape'), RefIDs=NA){
  #A <- Rpraetor$LMs

  if (ShapeArray==FALSE){
    RefArray <- Mat2Array(RefShapes, LMdim = LMDim)
  } else {
    RefArray <- RefShapes
  }

  if (length(RefIDs)==1 && is.na(RefIDs)){
    RefIDs <- 1:dim(RefArray)[3]
  } else if (length(RefIDs)<dim(RefArray)[3]){
    stop('Error: RefIDs is not the same length as the number of specimens in RefShapes. Please check your data and make sure the RefIDs correspond with and are in the same order as the specimens in RefShapes')
  }




  Res <- matrix(0, nrow = dim(RefArray)[3], ncol = dim(RefArray)[3])
  Res[lower.tri(Res)] <- 1
  Pairedindex <- which(Res==1, arr.ind = TRUE)


  DistMatloop <- function(X, index){
    ProcDRes <- shapes::procdist(X[,,index[1]], X[,,index[2]], type=Method)
    return(c(index, ProcDRes))
  }

  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)


  a <- 1
  ProcDResComp <- foreach::foreach(a = 1:dim(Pairedindex)[1], .combine = rbind) %dopar%{
    DistMatloop(X=RefArray, index = Pairedindex[a,])
  }

  parallel::stopCluster(clust)

  Res[lower.tri(Res)] <- as.numeric(ProcDResComp[,3])
  ProcDrableRes <- as.matrix(stats::as.dist(Res))

  colnames(ProcDrableRes) <- rownames(ProcDrableRes) <- RefIDs
  return(ProcDrableRes)
}

