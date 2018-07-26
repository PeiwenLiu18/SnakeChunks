
################################################################
#' @title DESeq2 analysis
#' @author Jacques van Helden
#' @description  Detect differentially expressed genes (DEG) using the package DESeq2,
#' and add a few custom columns (e-value, ...) to the result table.
#' @param counts a count table sent to DESeq2. Must contain raw counts (not normalized).
#' @param condition a vector with the condition associated to each sample. The length of this vector must equal the number of columns of the count table.
#' @param ref.condition=NULL reference condition for the differential analysis
#' @param comparison.prefix a string with the prefix for output files
#' @param title=comparison.prefix main title for the plots
#' @param dir.figures=NULL optional directory to save figures
#' @param ... additional parameters are passed to DESeq2::DESeq() function
#' @export
deseq2.analysis <- function(
  counts,
  condition,
  ref.condition=NULL,
  comparison.prefix,
  title = comparison.prefix,
  dir.figures=NULL,
  ...) {

  require(DESeq2)

  message("\tDESeq2 analysis\t", comparison.prefix)

  ## Check that the length of conditions equals the number of columns of the count table
  if (length(condition) != ncol(counts)) {
    stop("deseq2.analysis\tNumber of columns of count table (", ncol(counts), ") differs from length of condition (", length(condition), ").")
  }


  ## Create a DESeqDataSet object from the count table + condition
  deseq2.dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(condition),
    design = ~ condition)

  ## Indicate that second condition is the reference condition.
  ## If not done, the conditions are considered by alphabetical order,
  ## which may be misleading to interpret the log2 fold changes.
  if (!is.null(ref.condition)) {
    if (!(ref.condition %in% condition)) {
      stop("deseq2.analysis\treference condition (", ref.condition, ") does not exist in sample conditions (",
           paste(collapse = ", ", unique(condition)), ")")
    }
    deseq2.dds$condition <- relevel(deseq2.dds$condition, ref = ref.condition)
  }

  ## Run the differential analysis
  deseq2.dds <- DESeq(deseq2.dds)      ## Differential analysis with negbin distrib
  deseq2.res <- results(deseq2.dds, independentFiltering = FALSE, pAdjustMethod = "BH")  ## Collect the result table

  # names(deseq2.res)
  # VolcanoPlot(multitest.table = deseq2.res, effect.size.col = "log2FoldChange", control.type = "padj", alpha = 0.05, effect.threshold = 1.2)
  # plot(x = deseq2.res$log2FoldChange, y = -log10(deseq2.res$pvalue))

  deseq2.result.table <- data.frame(
    "gene.id" = row.names(deseq2.res),
    "mean" = deseq2.res$baseMean,
    "log2FC" = deseq2.res$log2FoldChange,
    "pvalue" = deseq2.res$pvalue,
    "padj" = deseq2.res$padj)

  ## Add complementary statistics on the DEG table
  deseq2.result.table <- DEGtablePostprocessing(
    deg.table = deseq2.result.table,
    table.name = paste(sep = "_", "DESeq2", comparison.prefix),
    sort.column = "padj",
    thresholds = thresholds,
    round.digits = 3,
    dir.figures = dir.figures)

  result <- list(
    dds = deseq2.dds,
    result.table = deseq2.result.table
  )

  return(result)
}
