#' @RowStats compute row-wise statistics on a data.frame or matrix
#' @author Jacques van Helden
#' @param x data.frame or matrix
#' @return a data.frame with one row per column of the input matrix, and one column per statistics
RowStats <- function(x, verbose = 0) {
  if (verbose >= 1) {
    message("\t\t", "Computing row-wise statistics")
  }
  result <- data.frame(
    zeros = apply(x == 0, 2, sum, na.rm = TRUE),
    na.values = apply(is.na(x), 2, sum),
    infinite.values = apply(is.infinite(as.matrix(x)), 2, sum),
    non.null = apply(x > 0, 2, sum, na.rm = TRUE), 
    sum = apply(x, 2, sum, na.rm = TRUE),
    mean = apply(x, 2, mean, na.rm = TRUE),
    var = apply(x, 2, var, na.rm = TRUE),
    sd = apply(x, 2, sd, na.rm = TRUE),
    min = apply(x, 2, min, na.rm = TRUE),
    perc01 = apply(x, 2, quantile, na.rm = TRUE, probs = 0.01),
    perc05 = apply(x, 2, quantile, na.rm = TRUE, probs = 0.05),
    perc10 = apply(x, 2, quantile, na.rm = TRUE, probs = 0.10),
    Q1 = apply(x, 2, quantile, na.rm = TRUE, probs = 0.25),
    median = apply(x, 2, quantile, na.rm = TRUE, probs = 0.5),
    Q3 = apply(x, 2, quantile, na.rm = TRUE, probs = 0.75),
    perc90 = apply(x, 2, quantile, na.rm = TRUE, probs = 0.90),
    perc95 = apply(x, 2, quantile, na.rm = TRUE, probs = 0.95),
    perc99 = apply(x, 2, quantile, na.rm = TRUE, probs = 0.90),
    max = apply(x, 2, max, na.rm = TRUE)
  )
  result$max.mean.ratio <- result$max/result$mean
  result$max.sum.ratio <- result$max/result$sum
  result$median.mean.ratio <- result$median/result$mean
  result$below.mean <- apply(t(x) < result$mean, 1, sum, na.rm = TRUE)
  result$fract.below.mean <- result$below.mean/nrow(x)
  return(result)  
}
