#' @title compute margin statistics (per row or column) on a data.frame or matrix.
#' @author Jacques van Helden
#' @param x data.frame or matrix
#' @param margin supported values: 1 (row-wise statistics) or 2 (column-wise)
#' @param verbose=0 level of verbosity
#' @param selected.stats=NULL if specified, restrict the computation to a subset of the supported statistics.
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
#'  \item{Mcounts}{Million counts = round(digits=2, sum/1e6). Classically used for Next-Generation Sequencing libraries.  }
#' }
#' @export
MarginStats <- function(x,
                        margin,
                        verbose = 0,
                        selected.stats=NULL) {

  ## Check margin parameter
  if (margin == 1) {
    margin.name <- "row"
  } else if (margin == 2) {
    margin.name <- "column"
  } else {
    stop('invalid margin specification for MarginStats(). Supported: 1 (row-wise) or 2 (column-wise)')
  }

  ## Check selected statistics
  supported.stats <- c("mean",
                       "median",
                       "sum",
                       "sd",
                       "iqr",
                       "var",
                       "min",
                       "max",
                       "perc01",
                       "perc05",
                       "perc10",
                       "Q1",
                       "Q3",
                       "perc90",
                       "perc95",
                       "perc99",
                       "zeros",
                       "na.values",
                       "infinite.values",
                       "non.null",
                       "sdiqr",
                       "max.mean.ratio",
                       "max.sum.ratio",
                       "mean.median.ratio",
                       "below.mean",
                       "fract.below.mean",
                       "Mcounts"
  )
  if (is.null(selected.stats)) {
    selected.stats <- supported.stats
  } else {
    non.supported <- setdiff(selected.stats, supported.stats)
    if (length(non.supported) > 0) {
      stop("MarginStats()\tInvalid selection of statistics: ",
           paste(collapse = ", ", non.supported),
           "\n\tSupported statistics: ", paste(collapse = ", ", supported.stats)
      )
    }
  }

  ## Verbosity
  if (verbose >= 1) {
    message("\t\t", "Computing ", margin.name, "-wise statistics")
    message("\t\tStats to compute: ", paste(collapse = ", ", selected.stats))
  }


  ## Compute margin statistics
  if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise mean") }
  margin.stats <- data.frame(
    mean = apply(x, margin, mean, na.rm = TRUE)
  )
  if ("median" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise median") }
    margin.stats$median = apply(x, margin, quantile, na.rm = TRUE, probs = 0.5)
  }
  if (length(intersect(selected.stats, c("sum", "Mcounts", "max.sum.ratio"))) > 0) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise sum") }
    margin.stats$sum = apply(x, margin, sum, na.rm = TRUE)
  }
  if ("sd" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise sd") }
    margin.stats$sd = apply(x, margin, sd, na.rm = TRUE)
  }
  if ("iqr" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise iqr") }
    margin.stats$iqr = apply(x, margin, IQR, na.rm = TRUE)
  }
  if ("var" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise var") }
    margin.stats$var = apply(x, margin, var, na.rm = TRUE)
  }
  if ("min" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise min") }
    margin.stats$min = apply(x, margin, min, na.rm = TRUE)
  }
  if ("max" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise max") }
    margin.stats$max = apply(x, margin, max, na.rm = TRUE)
  }
  if ("perc01" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise perc01") }
    margin.stats$perc01 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.01)
  }
  if ("perc05" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise perc05") }
    margin.stats$perc05 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.05)
  }
  if ("perc10" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise perc10") }
    margin.stats$perc10 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.10)
  }
  if ("Q1" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise Q1") }
    margin.stats$Q1 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.25)
  }
  if ("Q3" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise Q3") }
    margin.stats$Q3 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.75)
  }
  if ("perc90" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise perc90") }
    margin.stats$perc90 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.90)
  }
  if ("perc95" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise perc95") }
    margin.stats$perc95 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.95)
  }
  if ("perc99" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise perc99") }
    margin.stats$perc99 = apply(x, margin, quantile, na.rm = TRUE, probs = 0.90)
  }
  if ("zeros" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise zeros") }
    margin.stats$zeros = apply(x == 0, margin, sum, na.rm = TRUE)
  }
  if ("na.values" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise na.values") }
    margin.stats$na.values = apply(is.na(x), margin, sum)
  }
  if ("infinite.values" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise infinite.values") }
    margin.stats$infinite.values = apply(is.infinite(as.matrix(x)), margin, sum)
  }
  if ("non.null" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise non.null") }
    margin.stats$non.null = apply(x > 0, margin, sum, na.rm = TRUE)
  }
  if ("sdiqr" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise sdiqr") }
    margin.stats$sdiqr = margin.stats$iqr/(2*qnorm(3/4))
  }
  if ("max.mean.ratio" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise max.mean.ratio") }
    margin.stats$max.mean.ratio <- margin.stats$max/margin.stats$mean
  }
  if ("max.sum.ratio" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise max.sum.ratio") }
    margin.stats$max.sum.ratio <- margin.stats$max/margin.stats$sum
  }
  if ("mean.median.ratio" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise mean.median.ratio") }
    margin.stats$mean.median.ratio <- margin.stats$mean / margin.stats$median
  }
  if (length(intersect(c("below.mean", "frac.below.mean"), selected.stats)) > 0) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise below.mean") }
    margin.stats$below.mean <- apply(t(x) < margin.stats$mean, 1, sum, na.rm = TRUE)
  }
  if ("fract.below.mean" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise fract.below.mean") }
    margin.stats$fract.below.mean <- margin.stats$below.mean/nrow(x)
  }
  if ("Mcounts" %in% selected.stats) {
    if (verbose >= 2) { message("\t\tComputing ", margin.name, "-wise Mcounts") }
    margin.stats$Mcounts <- round(digits = 3, margin.stats$sum/1e6)
  }

  # ## Assign same names as in data table
  # if (margin == 2) {
  #   rownames(margin.stats) = rownames(x)
  # } else if (margin == 1) {
  #   colnames(margin.stats) = colnames(x)
  # }

  ## Return selected stats in the specified order
  margin.stats <- margin.stats[, selected.stats]

  return(margin.stats)
}

#' @title compute row statistics on a data.frame or matrix.
#' @author Jacques van Helden
#' @param x data.frame or matrix
#' @param verbose=0 level of verbosity
#' @param selected.stats=NULL if specified, restrict the computation to a subset of the supported statistics.
#' @description statistics are computed on each row by passing the data frame/matrux to MarginStats() with marrgin=1.
#' @export
RowStats <- function(x, verbose = 0) {
  MarginStats(x, 1, verbose = 0, selected.stats = NULL)
}

#' @title compute column statistics on a data.frame or matrix.
#' @author Jacques van Helden
#' @param x data.frame or matrix
#' @param verbose=0 level of verbosity
#' @param selected.stats=NULL if specified, restrict the computation to a subset of the supported statistics.
#' @description statistics are computed on each column by passing the data frame/matrux to MarginStats() with marrgin=2.
#' @export
ColStats <- function(x, verbose = 0, selected.stats = NULL) {
  MarginStats(x, 2, verbose, selected.stats)
}


#' @title Descriptive statistics on each row of the input matrix
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Compute descriptive parameters (central tendency, dispersion)
#' for each row of a matrix or data frame.
#'
#' @details
#' First version: 2015-03
#' Replaced by MarginStats: 2018-07-27
#'
#' @param x       A matrix or data frame
#'
#' @return
#' A data.frame with one row per row of the input matrix, and one column
#' per computed statistics.
#' @examples
#' ## Load example data set from Den Boer, 2009
#' library(denboer2009)
#' data(denboer2009.expr)     ## Load expression table
#' data(denboer2009.pheno)    ## Load phenotypic data
#' data(denboer2009.group.labels)    ## Load phenotypic data
#'
#' stats.per.row <- RowStats(denboer2009.expr, verbose = 1)
#' head(stats.per.row)
#' @export
statsPerRow <- function(x,
                        margin,
                        verbose = 0,
                        selected.stats = NULL) {
  message("statsPerRow() has been replaced by RowStats()")
  return(RowStats(x, 1, verbose, selected.stats))
}

