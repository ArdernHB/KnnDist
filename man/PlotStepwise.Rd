% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting_functions.R
\name{PlotStepwise}
\alias{PlotStepwise}
\title{Plotting function for stepwise results}
\usage{
PlotStepwise(
  StepwiseResultsMat,
  Percentiles = c(0.05, 0.95),
  PlotCol = "darkblue",
  Add = FALSE,
  Xlabel = "K"
)
}
\arguments{
\item{StepwiseResultsMat}{a matrix where the rows are the cross-validation percentage result for each iteration of resampling to equal sample size and the columns are the stepwise increase in k.}

\item{Percentiles}{a vector of two values denoting the upper and lower limit of the resampling to be plotted. Values must be between 0 and 1. Default is .05 and .95.}

\item{PlotCol}{is the colour of the line and the polygon representing the range of cross-validation percentage at the user defined percentiles.}

\item{Add}{is a logical value to determine if the the plot should be added to a previous plot. For example if the \code{StepwiseResultsMat} is the weighted result and if the user would like to plot the results of an unweighted analysis on the same plot for comparison.}

\item{Xlabel}{is the label to be used on the x axis. Default is 'K'; however, this function can be used for stepwise analysis of PCA and so the user may wish to change the x label.}
}
\value{
Plots a graph of mean CCV percentages with the percentiles plotted as a polygon range around the mean.
}
\description{
This function plots the results of a stepwise analyses with resampling to equal sample size.
The function will plot the line connecting the means of each stepwise increment increase of K
and will also plot the range around the mean. The default upper and lower limit is set to the 5th
and 95th percentiles.
}
\author{
Ardern Hulme-Beaman
}
\keyword{plotting}
