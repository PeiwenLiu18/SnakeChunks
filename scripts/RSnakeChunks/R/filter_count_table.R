#' @title Filter features of an RNA-seq count table based on a combination of user-selected criteria. 
#' @author Jacques van Helden
#' @description Filter out the features of an RNA-seq count table based on several filters. 
#'
#' @param counts  A data frame with counts (one row per feature and one column per sample).
#' @param condition=NULL An optional vector indicating the  groups (condition, genotype, treatment) of each sample. Must have the same length as the number of columns of the count table. 
#' @param na.rm=FALSE Remove rows with at least one NA value.
#' @param mean.count=NULL Filter out the features having a cross-sample mean count lower than the specified value
#' @param min.count=NULL Filter out features for which the  max count (across all samples) is lower than the specified threshold.
#' @param mean.per.condition=NULL  ## Filter out features for which none of the conditions has a mean count above the threshold
#' @param black.list=NULL Black-listed" features: a vector of row IDs or row indices indicating a list of features to be filtered out (e.g. 16S, 12S RNA).
#' @return a data.frame with the same columns (samples) as the input count table, and with rows restricted to the features that pass all the thresholds.
#' @export
FilterCountTable <- function(counts,
                             condition = NULL,
                             na.omit = FALSE,
                             mean.count = NULL,
                             min.count = NULL,
                             mean.per.condition = NULL,
                             black.list = NULL
                            ) {
  message("Filtering count table with ", nrow(counts), " features x ", ncol(counts), " samples. ")

  if (na.omit) {
    na.rows <- apply(is.na(counts), 1, sum) > 0
    counts <- counts[!na.rows, ]
    message("\tomitted NA values; discarded features: ", sum(na.rows), "; kept features: ", nrow(counts))
    
  }
  
  ## Min count per rpw  
  if (!is.null(min.count)) {
    discarded <- apply(counts, 1, min, na.rm = TRUE) < min.count
    counts <- counts[!discarded, ]
    message("\tMin count per row >= ", min.count, "; discarded features: ", sum(discarded), "; kept features: ", nrow(counts))
  }
  
  ## Mean count per rpw  
  if (!is.null(mean.count)) {
    discarded <- apply(counts, 1, mean, na.rm = TRUE) < mean.count
    counts <- counts[!discarded, ]
    message("\tMean count per row >= ", mean.count, "; discarded features: ", sum(discarded), "; kept features: ", nrow(counts))
  }


  ## Filter out black-listed features
  if (!is.null(black.list)) {
    kept <- setdiff(row.names(counts), black.list)
    counts <- counts[kept, ]
    message("\tBlack-listed features: ", length(black.list), "; kept features: ", nrow(counts))
  }
  
  message("Returning filtered count table with ", nrow(counts), " features x ", ncol(counts), " samples. ")
  return(counts)
}
