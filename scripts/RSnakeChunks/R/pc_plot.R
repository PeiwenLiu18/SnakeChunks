#' @title Compute principal components and draw a plot with the first components
#' @author Jacques van Helden & Fanny Casse
#' @description taking as input a data frame, compute principal components and draw XY plot with user-selected component numbers (default: PC1 and PC2). 
#' @param x a data frame with genes as rows and samples as columns
#' @param PCs=c(1,2) user-selected components for the XY plot
#' @param main="PC plot" # Main title for the plot
#' @return the object of prcomp()
PCplot <- function(x,
                   PCs = c(1,2),
                   sample.desc,
                   main = "PC plot",
                   sample.labels=row.names(sample.desc),
                   sample.col = sample.desc$color
                   ) {
  x.pc <- prcomp(t(x))
  plot(x.pc$x[,PCs], panel.first = grid(), type = "n", main = main)
  text(x.pc$x[,PCs], labels = sample.labels, col=sample.col)
}


