
################################################################
#' @title edgeR analysis
#' @author Jacques van Helden
#' @description  Detect differentially expressed genes (DEG) using the package edgeR,
#' and add a few custom columns (e-value, ...) to the result table.
#' @param counts a count table sent to edgeR Must contain raw counts (not normalized).
#' @param condition a vector with the condition associated to each sample. The length of this vector must equal the number of columns of the count table.
#' @param ref.condition=NULL reference condition for the differential analysis
#' @param comparison.prefix a string with the prefix for output files
#' @param title=comparison.prefix main title for the plots
#' @param dir.figures=NULL optional directory to save figures
#' @param norm.method="RLE" normalisation method. This parameter strongly affects the results! See edgeR documentation for a list of supported methods
#' @param ... additional parameters are passed to edgeR::exactTest() function
#' @export
edger.analysis <- function(counts,
                           condition,
                           ref.condition = NULL,
                           comparison.prefix,
                           norm.method = "TMM",
                           title = paste(sep = "_", norm.method, comparison.prefix),
                           dir.figures = NULL,
                           ...) {

  require(edgeR)

  message("\tedgeR analysis\t", comparison.prefix, "\tnormalisation method: ", norm.method)

  ## Check that the length of conditions equals the number of columns of the count table
  if (length(condition) != ncol(counts)) {
    stop("edgeR.analysis\tNumber of columns of count table (", ncol(counts), ") differs from length of condition (", length(condition), ").")
  }


  ## Convert the count table in a DGEList structure and compute its parameters.
  # d <- DGEList(counts = current.counts, group = sample.conditions[names(current.counts)])
  # d$samples$group <- relevel(d$samples$group, ref = ref.condition) ## Ensure that condition 2 is considered as the reference
  # d <- calcNormFactors(d, method="RLE")                 ## Compute normalizing factors
  # d <- estimateCommonDisp(d, verbose=FALSE)             ## Estimate common dispersion
  # d <- estimateTagwiseDisp(d, verbose=FALSE)            ## Estimate tagwise dispersion
  d <- DGEList(counts = counts, group = condition)
  d$samples$group <- relevel(d$samples$group, ref = ref.condition) ## Ensure that condition 2 is considered as the reference
  d <- calcNormFactors(d, method = norm.method)                 ## Compute normalizing factors
  d <- estimateCommonDisp(d, verbose = FALSE)             ## Estimate common dispersion
  d <- estimateTagwiseDisp(d, verbose = FALSE)            ## Estimate tagwise dispersion

  ################################################################
  ## Detect differentially expressed genes by applying the exact
  ## negative binomial test from edgeR package.
  edger.de <- exactTest(d, pair = c(cond2, cond1), ...)      ## Run the exact negative binomial test

  ## Sort genes by increasing p-values, i.e. by decreasing significance
  edger.tt <- topTags(edger.de, n = nrow(d), sort.by = "PValue")

  ## Complete the analysis of edgeR result table
  edger.result.table <- data.frame("gene.id" = row.names(edger.tt$table),
                                   "mean" = edger.tt$table$logCPM,
                                   "log2FC" = edger.tt$table$logFC,
                                   "pvalue" = edger.tt$table$PValue,
                                   "padj" = edger.tt$table$FDR)
  edger.result.table <- DEGtablePostprocessing(
    deg.table = edger.result.table,
    table.name = paste(sep = "_", "edgeR", norm.method, comparison.prefix),
    sort.column = "padj",
    thresholds = thresholds,
    round.digits = 3,
    dir.figures = dir.figures)
  # dim(edger.result.table)
  # dim(edger.result.table)
  # names(edger.result.table)
  # View(edger.result.table)

  result <- list(
    edger.d = d,
    edger.de = edger.de,
    edger.tt = edger.tt,
    result.table = edger.result.table
  )

  par(par.ori)
  return(result)
}
