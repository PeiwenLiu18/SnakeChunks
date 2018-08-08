
#' @title Compute sample-wise statistics for a count table
#' @author Jcaques van Helden
#' @param sample.descriptions a data.frame with one row per sample, in the same order as the columns of the count table.
#' @param counts a data.frame with one row per feature (gene, transcript, ...) and one column per sample. The columns of the count table must correspond to the rows of the sample descrtiption file.
#' @param verbose=1 verbosity
calc.stats.per.sample <- function(sample.descriptions,
                                  counts, verbose = 1) {
  if (verbose >= 2) { message("Computing statistics per sample") }

  stats.per.sample <- cbind(
    sample.descriptions[names(counts), ],
    data.frame(
      "zeros" = apply(counts == 0, 2, sum, na.rm = TRUE), ## Number of genes with 0 counts
      "detected" = apply(counts > 0, 2, sum, na.rm = TRUE), ## Number of genes counted at least once
      "sum" = apply(counts, 2, sum, na.rm = TRUE), ## Sum of all counts for the sample
      "mean" = apply(counts, 2, mean, na.rm = TRUE), ## Mean counts per gene
      "min" = apply(counts, 2, min, na.rm = TRUE), ## Min counts per gene
      "perc05" = apply(counts, 2, quantile, probs = 0.05, na.rm = TRUE), ## 5th percentile
      "perc25" = apply(counts, 2, quantile, probs = 0.25, na.rm = TRUE), ## 25th percentile
      "median" = apply(counts, 2, median, na.rm = TRUE), ## median (percentile 50)
      "perc75" = apply(counts, 2, quantile, probs = 0.75, na.rm = TRUE), ## percentile 75
      "perc90" = apply(counts, 2, quantile, probs = 0.90, na.rm = TRUE), ## percentile 95
      "perc95" = apply(counts, 2, quantile, probs = 0.95, na.rm = TRUE), ## percentile 95
      "max" = apply(counts, 2, max, na.rm = TRUE) ## Max counts per gene
    )
  )
  rownames(stats.per.sample) <- colnames(counts)
  stats.per.sample$max.sum.ratio <- stats.per.sample$max / stats.per.sample$sum
  stats.per.sample$median.mean.ratio <- stats.per.sample$media / stats.per.sample$mean
  stats.per.sample$Mreads <- round(stats.per.sample$sum/1e6, digits = 1)

  ## Count number and the fraction of samples with counts below the mean.
  ## This shows the impact of very large counts: in some samples,
  ## 85% of the samples have a value below the mean (i.e. the mean is at the percentile 85 !)
  stats.per.sample$below.mean <- apply(t(counts) < stats.per.sample$mean, 1, sum, na.rm = TRUE)
  stats.per.sample$fract.below.mean <- stats.per.sample$below.mean/nrow(counts)
  # View(stats.per.sample)

  return(stats.per.sample)
}
