
################################################################
#' @title Initiate a result table with the input counts and derived statistics.
#' @description Initiate a result table with the input counts and derived statistics for the whole data table + for the two sample sets selected for differential analysis.
#' @param count.table table with the counts per reads. This table is used to get feature IDs
#' (gene IDs, gene names, ...) from row names, and thereby ensure consistency between rows
#' of the count tables and the summary of differential expression.
#' @param samples1 vector with  sample IDs of the test condition. Must be a subset of the data table column names.
#' @param samples2 vector with sample IDs of the control condition. Must be a subset of the data table column names.
#' @param gene.info optional table with detailed information about each feature (gene).
#' @param stats=FALSE add columns with statistics (mean, quartiles, var, ...) for the whole count table and for the respective sample groups
#' @return a data frame with one row per feature and some descriptive columns + optional table-wise and group-wise statistics.
init.deg.table <- function(count.table,
                           samples1,
                           samples2,
                           gene.info = NULL,
                           stats = FALSE) {
  message("\tInitializing result table for one differential analysis (two-sample comparison).")
  all.gene.ids <- row.names(count.table)
  #message("\t\tFeatures (genes): " \t, )

  ## Build a minimal gene info table if not provided
  if (is.null(gene.info)) {
    gene.info <- data.frame(
      "id" = all.gene.ids,
      "name" = all.gene.ids,
      "entrez.id" = "",
      "description" = paste("gene_id:", all.gene.ids)
    )
  } else {
    gene.info <- gene.info[all.gene.ids, ] ## Make sure gene.info is in the same order as the count table
  }

  ## Check consistency betwen test/control sample IDs and column names of the data table
  if (length(setdiff(samples1, colnames(count.table))) > 0) {
    stop("Test sample IDs absent from count table column names: ",
         paste(collapse = ", ", setdiff(samples1, colnames(count.table))))
  }

  if (length(setdiff(samples2, colnames(count.table))) > 0) {
    stop("Control sample IDs absent from count table column names: ",
         paste(collapse = ", ", setdiff(samples2, colnames(count.table))))
  }

  ## Build result table
  result.table <- data.frame("gene_id" = all.gene.ids,
                             "name" = gene.info[,"name"])
  row.names(result.table) <- all.gene.ids
  result.table$entrez.id <- gene.info[,"entrez.id"]
  result.table$description <- gene.info[,"description"]

  result.table <- cbind(result.table, count.table) ## Include the original counts in the big result table
  # View(result.table)
  #  result.table <- cbind(result.table, cpm = current.cpms) ## Include CPMs in the big result table

  ## Tag genes detected in less than min.rep samples, which is defined as
  ## the minimal number of replicates per condition.
  min.rep <- min(length(samples1), length(samples2))
  result.table$undetected <- rowSums(count.table > 1) < min.rep
  # table(result.table$undetected)
  # dim(count.table)

  if (stats) {
    message("\tComputing group-wise count statistics")
    result.table$mean <- apply(count.table,1, mean)
    result.table$mean1 <- apply(as.data.frame(count.table[,samples1]),1, mean)
    result.table$mean2 <- apply(as.data.frame(count.table[,samples2]),1, mean)
    result.table$M = log2(result.table$mean1/result.table$mean2)
    result.table$min <-  apply(count.table,1, min)
    result.table$min1 <- apply(as.data.frame(count.table[,samples1]),1, min)
    result.table$min2 <- apply(as.data.frame(count.table[,samples2]),1, min)
    result.table$perc25 <- apply(count.table,1, quantile, probs = 0.75)
    result.table$perc25.1 <- apply(as.data.frame(count.table[,samples1]),1, quantile, probs = 0.75)
    result.table$perc25.2 <- apply(as.data.frame(count.table[,samples2]),1, quantile, probs = 0.75)
    result.table$median <- apply(count.table,1, median)
    result.table$median1 <- apply(as.data.frame(count.table[,samples1]),1, median)
    result.table$median2 <- apply(as.data.frame(count.table[,samples2]),1, median)
    result.table$perc75 <- apply(count.table,1, quantile, probs = 0.75)
    result.table$perc75.1 <- apply(as.data.frame(count.table[,samples1]),1, quantile, probs = 0.75)
    result.table$perc75.2 <- apply(as.data.frame(count.table[,samples2]),1, quantile, probs = 0.75)
    result.table$perc95 <- apply(count.table,1, quantile, probs = 0.75)
    result.table$perc95.1 <- apply(as.data.frame(count.table[,samples1]),1, quantile, probs = 0.95)
    result.table$perc95.2 <- apply(as.data.frame(count.table[,samples2]),1, quantile, probs = 0.95)
    result.table$max <-  apply(count.table,1, max)
    result.table$max1 <- apply(as.data.frame(count.table[,samples1]),1, max)
    result.table$max2 <- apply(as.data.frame(count.table[,samples2]),1, max)
    result.table$sd <-  apply(count.table,1, sd)
    result.table$sd1 <- apply(as.data.frame(count.table[,samples1]),1, sd)
    result.table$sd2 <- apply(as.data.frame(count.table[,samples2]),1, sd)
    result.table$var <-  apply(count.table,1, var)
    result.table$var1 <- apply(as.data.frame(count.table[,samples1]),1, sd)
    result.table$var2 <- apply(as.data.frame(count.table[,samples2]),1, sd)
  }
  return(result.table)
}
