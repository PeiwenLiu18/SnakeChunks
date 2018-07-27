#' @title compute margin statistics (per row or column) on a data.frame or matrix.
#' @author Jacques van Helden
#' @param x data.frame or matrix
#' @param margin supported values: 1 (row-wise statistics) or 2 (column-wise)
#' @param verbose=0 level of verbosity
#' @return a data.frame with one row per entry (row or column) of the input matrix, and one column per statistics.
#' #' \describe{
#'  \item{mean}{arithmetic mean}
#'  \item{median}{}
#'  \item{sum}{}
#'  \item{iqr}{inter-quartile range}
#'  \item{sd}{estimation of population standard deviation returned by R sd() function. }
#'  \item{var}{estimation of population variance returned by R var() function. }
#'  \item{min}{}
#'  \item{perc01}{first percentile}
#'  \item{perc05}{5th percentile}
#'  \item{perc10}{10th percentile}
#'  \item{Q1}{first quartile, also called lower quartile (LQ) or 25th percentile}
#'  \item{Q3}{third quartile, also called upper quartile (UQ) or 75th percentile}
#'  \item{perc90}{90th percentile}
#'  \item{perc95}{95th percentile}
#'  \item{perc99}{99th percentile}
#'  \item{zeros}{number of 0 values}
#'  \item{na.values}{number of NA values}
#'  \item{inf.values}{number of infinite values}
#'  \item{non.null}{number of values that are neither null nor NA nor Inf. }
#'  \item{non.null}{number of values that are neither null nor NA nor Inf. }
#'  \item{sdiqr}{standardized IQR (data IQR divided by IQR of the Normal function), which can serve as robust estimator of the standard deviation. }
#'  \item{max.mean.ratio}{ratio between max and mean values. }
#'  \item{max.sum.ratio}{ratio between max and sum. }
#'  \item{mean.median.ratio}{ratio between mean and the median, which can be interpreted as skew indicator. }
#'  \item{below.mean}{number of values below the mean. }
#'  \item{fract.below.mean}{fraction of values below the mean. }
#' }
#' @export
MarginStats <- function(x, margin, verbose = 0) {
    if (margin == 1) {
    margin.name <- "row"
  } else if (margin == 2) {
    margin.name <- "column"
  } else {
    stop('Invalid margin specification for MarginStats(). Supported: 1 (row-wise) or 2 (column-wise)')
  }
  if (verbose >= 1) {
    message("\t\t", "Computing ", margin.name, "-wise statistics")
  }
  margin.stats <- data.frame(
    mean = apply(x, margin, mean, na.rm = TRUE),
    median = apply(x, margin, quantile, na.rm = TRUE, probs = 0.5),
    sum = apply(x, margin, sum, na.rm = TRUE),
    sd = apply(x, margin, sd, na.rm = TRUE),
    iqr = apply(x, margin, IQR, na.rm = TRUE),
    var = apply(x, margin, var, na.rm = TRUE),
    min = apply(x, margin, min, na.rm = TRUE),
    max = apply(x, margin, max, na.rm = TRUE),
    perc01 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.01),
    perc05 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.05),
    perc10 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.10),
    Q1 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.25),
    Q3 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.75),
    perc90 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.90),
    perc95 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.95),
    perc99 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.90),
    zeros = apply(x == 0, margin, sum, na.rm = TRUE),
    na.values = apply(is.na(x), margin, sum),
    infinite.values = apply(is.infinite(as.matrix(x)), margin, sum),
    non.null = apply(x > 0, margin, sum, na.rm = TRUE)
  )
  margin.stats$sdiqr = margin.stats$iqr/(2*qnorm(3/4))
  margin.stats$max.mean.ratio <- margin.stats$max/margin.stats$mean
  margin.stats$max.sum.ratio <- margin.stats$max/margin.stats$sum
  margin.stats$median.mean.ratio <- margin.stats$median/margin.stats$mean
  margin.stats$below.mean <- apply(t(x) < margin.stats$mean, 1, sum, na.rm = TRUE)
  margin.stats$fract.below.mean <- margin.stats$below.mean/nrow(x)

  ## Assign same names as in data table
  if (margin == 1) {
    rownames(margin.stats) = rownames(x)
  } else if (margin == 2) {
    colnames(margin.stats) = colnames(x)
  }

  return(margin.stats)
}

#' @title compute row statistics on a data.frame or matrix.
#' @author Jacques van Helden
#' @param x data.frame or matrix
#' @param verbose=0 level of verbosity
#' @description statistics are computed on each row by passing the data frame/matrux to MarginStats() with marrgin=1.
#' @export
RowStats <- function(x, verbose = 0) {
  MarginStats(x, 1, verbose = 0)
}

#' @title compute column statistics on a data.frame or matrix.
#' @author Jacques van Helden
#' @param x data.frame or matrix
#' @param verbose=0 level of verbosity
#' @description statistics are computed on each column by passing the data frame/matrux to MarginStats() with marrgin=2.
#' @export
ColStats <- function(x, verbose = 0) {
  MarginStats(x, 2, verbose = 0)
}
