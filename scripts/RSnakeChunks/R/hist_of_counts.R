#' @title histogram of RNA-seq counts
#' @author Jacques van Helden
#' @description draw an histogram with suitable parameters to RNA-seq count table.
#' @param counts RNA-seq count table or related dataset.
#' @param breaks=100 approximate number of breaks
#' @param maxPercentile=95 maximum percentile to display on the X axis, in order to avoid the effect of outliers. 
#' @param discardZeroRows=FALSE if TRUE, discard rows containing only zero values before drawing the histogram.
#' @param discardZeros=FALSE if TRUE, discard all zero values before drawing the histogram.
#' @param main="Count histogram" main title.
#' @param ylab="Features" label for Y axis.
#' @param xlab=paste(sep = "", "Counts (axis truncated to percentile ", maxPercentile,")") label for X axis.
#' @param legend.cex=0.85 font size for the legend
#' @param drawMedian = TRUE draw an arrow to mark the median value
#' @param drawMean = TRUE draw an arrow to mark the mean value
#' @param verbose = 0 level of verbosity
#' @param ... additional parameters are passed to hist()
#' @return a list with the computed parameters + the histogram data
#' @export
HistOfCounts <- function(counts,
                         breaks = 100,
                         maxPercentile = 95,
                         discardZeroRows = FALSE,
                         discardZeros = FALSE,
                         main = "Count histogram", 
                         ylab = "Features",
                         xlab = paste(sep = "", "Counts (axis truncated to percentile ", maxPercentile,")"),
                         legend.cex = 0.85,
                         verbose = 0,
                         ...) {
  
  nbFeatures <- nrow(counts)
  nbSamples <- ncol(counts)
  if (verbose >= 1) { message("\thistOfCounts()\t", nbFeatures, " features x ", nbSamples, " samples.") }
  
  ## Count the number of zero values
  nbZeroValues <- sum(!is.na(counts) & (counts == 0))
  
  ## Count number of NA values
  nbNA <- sum(is.na(counts))

  ## Count number of rows with only zeros, and discard them if requested  
  zerosPerRow <- apply(counts == 0, 1, sum, na.rm = TRUE) 
  # hist(zerosPerRow, breaks <- 0:nbSamples)
  zeroOnlyRows <- zerosPerRow == nbSamples
  nbZeroOnlyRows <- sum(zeroOnlyRows)
  if (discardZeroRows) {
    if (verbose >= 2) { message("Discarding rows where all values equal 0") }
    counts <- counts[!zeroOnlyRows,]
  }
    
  x <- as.vector(unlist(counts)) ## Flatten the count table
  x <- na.omit(x) ##  Discard NA values
  if (discardZeros) {
    if (verbose >= 2) { message("Discarding all 0 values") }
    # sum(x==0); length(x); mean(x)
    x <- x[x > 0]
    # sum(x==0); length(x); mean(x)
  }
  ## Restrict the histogram range to the percentile
  histMax <- quantile(x = x, p = maxPercentile/100)
  histBreaks <- pretty(x = x[x < histMax], n = breaks)
  histBreaks <- append(histBreaks, max(x) + 1)
  histData <- hist(x, freq = TRUE,
                   breaks = histBreaks, 
                   xlim = c(0, histMax),
                   main = main, 
                   xlab = xlab,
                   ylab = ylab,
                   ...)
  
  # yMax <- max(histData$counts)
  yMax <- par("usr")[4]
  xMedian <- median(x, na.rm = TRUE)
  xMean <- mean(x, na.rm = TRUE)
  
  if (xMedian < xMean) {
    medPos <- 3
    meanPos <- 4
  }  else {
    medPos <- 4
    meanPos <- 3
  }
  
  ## Draw an arrow at the median value
  arrows(x0 = xMedian, y0 = 0.9 * yMax, x1 = xMedian, y1 = 0.8 * yMax, angle = 30, length = 0.1, col = "#008800", lwd = 2)
  text(x = xMedian, y = 0.9 * yMax, pos = medPos, labels = paste(sep = "", "med=", round(digits=1, xMedian)), cex = 0.8, col =  "#008800")
  
  ## Draw an arrow at the mean value
  arrows(x0 = xMean, y0 = 0.9 * yMax, x1 = xMean, y1 = 0.8 * yMax, angle = 30, length = 0.1, col = "brown", lwd = 2)
  text(x = xMean, y = 0.9 * yMax, pos = meanPos, labels = paste(sep = "", "mean=", round(digits=1, xMean)), cex = 0.8, col = "brown")
  
  legend("topright", cex = legend.cex,
         legend = c(
           paste(sep = "", "Features: ", nbFeatures),
           paste(sep = "", "Samples: ", nbSamples),
           paste(sep = "", "Percentile ", maxPercentile, ": ", round(histMax)),
           paste(sep = "", "NA values: ", nbNA),
           paste(sep = "", "Zero values: ", format(nbZeroValues, big.mark = ",")),
           paste(sep = "", "Rows with all zeros: ", format(nbZeroOnlyRows, big.mark = ",")))
           )
  
  
  ## Prepare the result
  result <- list()
  result$maxPercentile <- maxPercentile
  result$discardZeroRows <- discardZeroRows
  result$discardZeros <- discardZeros
  result$nbFeatures <- nbFeatures
  result$nbZeroValues <- nbZeroValues
  result$nbZeroOnlyRows <- nbZeroOnlyRows
  result$nbNA <- nbNA
  result$histMax <- histMax
  result$histData <- histData
  return(result)
}
