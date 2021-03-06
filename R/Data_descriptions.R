
#' Rattini shape dataset
#'
#' A geometric morphometric derived dataset from the data of Hulme-Beaman et al.
#' 2019. The dataset contains an array of 395 specimens each with 105 landmarks.
#' Scale data is also included as is individual specimen info. A matrix of full
#' Procrustes distances calculated among the specimens is also provided.
#'
#' @format A list of 4 objects:
#' \describe{
#'   \item{LMs}{Specimen landmark data in array form, where each row is a landmark, each column a landmark dimension and each slice is an individual specimen.}
#'   \item{Scale}{A numeric vector of scale values for the landmarks.}
#'   \item{Info}{A dataframe corresponding to the LM data which includes specimen genus, species and centroid size (CS) data.}
#'   \item{FullProcrustesDistMat}{A square matrix of pairwise full procrustes distances calculated among also specimens.}
#' }
#' @source \url{https://doi.org/10.1007/s10914-017-9423-8}
"RatData"


#' Dummy Data dataset
#'
#' A dummy dataset for visualisation of difference between weighted and unweighted
#' k-NN analyses
#'
#' @format A list of 3 objects:
#' \describe{
#'   \item{Coords}{A matrix of dummy coordinates}
#'   \item{Groups}{A vector of dummy group classifications}
#'   \item{DistMat}{A distance matrix of Euclidean distances among dummy specimens}
#' }
"DummyData"

#' Blue chaffinch birdsong dataset
#'
#' A birdsong dissimilarity dataset calculated using dynamic time warping. The
#' dataset includes 116 specimens from across Tenerife. These represent at least
#' three regional birdsong groups, with a fourth dispersed group from across the
#' centre of the island. The fourth group may not constitute a real group.
#'
#' @format A list of 2 objects:
#' \describe{
#'   \item{SongDistMat}{A square matrix of dynamic time warping dissimilarities calculated among the songs of 116 individuals.}
#'   \item{Groups}{A dataframe of 4 columns including specimen IDs, Lat Long coordinates and Regional grouping.}
#' }
"BirdsongData"
