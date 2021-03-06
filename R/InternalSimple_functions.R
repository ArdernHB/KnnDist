
#' Internal function for quick class change of a vector to numeric
#'
#' This function takes a factor vector and changes it to a numeric.
#' This is primarily to avoid instances where a csv file has been opened
#' without column classes being specified, which can result in numeric
#' values being accidently read as factors, which in turn when treated
#' with as.numeric() will accidently be given integer values in the order
#' of the factor levels. By treating the values as characters first it avoids
#' this potential error.
#' @return
#'
#' This function returns a numeric vector
#'
#' @param x is a vector to be converted to a numeric vector
#'
#' @keywords internal
#' @author Ardern Hulme-Beaman

chr2nu <- function(X){as.numeric(as.character(X))}





#' Internal function for conversion of a matrix of landmarks to an array
#'
#' This function takes a matrix where rows are specimens and columns are
#' alternating landmark coordinate values and converts it to a 3 dimensional
#' array where rows are landmarks and columns correspond with the dimensions.
#' @return
#'
#' This function returns an array
#'
#' @param Mat is a matrix of landmark data
#' @param LMdim is the number of dimensions of the data, either 2 or 3
#'
#' @keywords internal
#' @author Ardern Hulme-Beaman
#' @export



Mat2Array <- function(Mat, LMdim){
  NewArray <- array(data = NA, dim = c(dim(Mat)[2]/LMdim, LMdim, dim(Mat)[1]))


  for (i in 1:dim(Mat)[1]){
    #i <- 1
    Mat4Array <- matrix(as.numeric(Mat[i,]), nrow = dim(Mat)[2]/LMdim, ncol = LMdim, byrow = TRUE)
    NewArray[,,i] <- Mat4Array
  }


  return(NewArray)
}



#' Internal function for conversion of an array to a matrix
#'
#' This function takes an array and converts it to a matrix where rows
#' specimens and columns are alternating landmark variables
#' (e.g. X1, X2, Y1, Y2...)
#' @return
#'
#' This function returns an array
#'
#' @param Array is an array of landmark data
#'
#' @keywords internal
#' @author Ardern Hulme-Beaman
#' @export


Array2Mat <- function(Array){
  Matrix <- matrix(NA, nrow = dim(Array)[3], ncol = length(c(t(Array[,,1]))))
  for (i in 1:dim(Array)[3]){
    #i <- 1
    Matrix[i,] <- c(t(Array[,,i]))
  }
  return(Matrix)
}


#' Internal function: Transparent named colour
#'
#' This function takes a named colour and returns the transparent equivalent
#' @param Colour A colour name from colours() function which is desired in transparent form.
#' @param alpha The level of transparency from 1 (completely transparent) to 100 (completely opaque) that the returned colour should be.
#' @return The transparent equivalent of a named colour
#' @keywords internal
#' @keywords colour
#' @keywords transparency
#' @author Ardern Hulme-Beaman


transpar<-function(Colour, alpha=100){
  newColour<-grDevices::col2rgb(Colour)
  apply(newColour, 2, function(curcoldata){grDevices::rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
