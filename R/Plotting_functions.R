


#' Plotting function for stepwise results
#'
#' This function plots the results of a stepwise analyses with resampling to equal sample size.
#' The function will plot the line connecting the means of each stepwise increment increase of K
#' and will also plot the range around the mean. The default upper and lower limit is set to the 5th
#' and 95th percentiles.
#'
#' @param Percentiles a vector of two values denoting the upper and lower limit of the resampling to be plotted. Values must be between 0 and 1. Default is .05 and .95.
#' @param StepwiseResultsMat a matrix where the rows are the cross-validation percentage result for each iteration of resampling to equal sample size and the columns are the stepwise increase in k.
#' @param PlotCol is the colour of the line and the polygon representing the range of cross-validation percentage at the user defined percentiles.
#' @param Add is a logical value to determine if the the plot should be added to a previous plot. For example if the \code{StepwiseResultsMat} is the weighted result and if the user would like to plot the results of an unweighted analysis on the same plot for comparison.
#' @param Xlabel is the label to be used on the x axis. Default is 'K'; however, this function can be used for stepwise analysis of PCA and so the user may wish to change the x label.
#' @return Plots a graph of mean CCV percentages with the percentiles plotted as a polygon range around the mean.
#'
#'
#' @keywords plotting
#' @author Ardern Hulme-Beaman
#'
#'
#' @export

PlotStepwise <- function(StepwiseResultsMat, Percentiles=c(.05, .95), PlotCol='darkblue', Add=FALSE, Xlabel='K'){

  #StepwiseResultsMat = VoleStepLDA; PlotCol = 'darkblue'; Add = FALSE

  DataDim <- dim(StepwiseResultsMat)[2]

  if (Add==FALSE){
    graphics::plot(y = colMeans(StepwiseResultsMat), x = 1:DataDim, type = 'n', ylim = c(10,105), ylab = 'CCV %', xlab = Xlabel)
    graphics::abline(h = seq(from = 20, to = 100, by = 10), v = seq(from = 2, to = DataDim, by =2), lty = '1919')

  }


  ResRange <- apply(StepwiseResultsMat, MARGIN = 2, FUN = stats::quantile, probs = Percentiles, na.rm=TRUE)
  graphics::polygon(x = c(1:DataDim, DataDim:1),
                    y = c(ResRange[1,], ResRange[2,DataDim:1])*100,
                    col = transpar(PlotCol, alpha = 75),
                    border = NA)

  graphics::lines(y = colMeans(StepwiseResultsMat*100), x = 1:DataDim, col=PlotCol, lwd=3)

}


