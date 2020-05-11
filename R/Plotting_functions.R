


#' Plotting function for stepwise results
#'
#' This function plots the results of a stepwise analyses with resampling to equal sample size.
#' The function will plot the line connecting the means of each stepwise incriment increase of K
#' and will also plot the range around the mean (default upper and lower limit is set to the 5th
#' and 95th percentiles).
#'
#' @param StepwiseResultsMat a
#' @param Probs a vector of two values denoting the upper and lower limit of the resampling to be plotted. Values must be between 0 and 1. Default is .05 and .95.
#' @param StepwiseResultsMat
#' @return
#' @details
#' @author Ardern Hulme-Beaman
#'
#' @export

PlotStepwise <- function(StepwiseResultsMat, Probs=c(.05, .95), PlotCol='darkblue', Add=FALSE, Xlabel='K'){

  #StepwiseResultsMat = VoleStepLDA; PlotCol = 'darkblue'; Add = FALSE

  DataDim <- dim(StepwiseResultsMat)[2]

  if (Add==FALSE){
    graphics::plot(y = colMeans(StepwiseResultsMat), x = 1:DataDim, type = 'n', ylim = c(10,105), ylab = 'CCV %', xlab = Xlabel)
    graphics::abline(h = seq(from = 20, to = 100, by = 10), v = seq(from = 2, to = DataDim, by =2), lty = '1919')

  }


  ResRange <- apply(StepwiseResultsMat, MARGIN = 2, FUN = stats::quantile, probs = Probs, na.rm=TRUE)
  graphics::polygon(x = c(1:DataDim, DataDim:1),
                    y = c(ResRange[1,], ResRange[2,DataDim:1])*100,
                    col = transpar(PlotCol, alpha = 75),
                    border = NA)

  graphics::lines(y = colMeans(StepwiseResultsMat*100), x = 1:DataDim, col=PlotCol, lwd=3)

}


