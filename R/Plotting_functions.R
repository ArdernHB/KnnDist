


#' Plotting function for stepwise results
#'
#' This function NEEDS COMPLETION
#'
#'
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


