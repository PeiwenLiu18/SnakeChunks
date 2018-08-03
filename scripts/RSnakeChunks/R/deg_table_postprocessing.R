#' @title Post-process differential expression analysis result.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Starting from a result table produced by DESeq2 or edgeR (or similar), to the following:
#' \itemize{
#'   \item add some relevant columns for further analysis and thresholding
#'   \item sort the table by increasing adjusted p-value
#'   \item optionally, filter the result table by applying thresholds on one or more user-selected columns
#' }
#'
#' @details
#' First version: 2015-08
#' Last modification: 2015-08.
#'
#' @param deg.table   Table with the results of differential expression analysis
#' obtained from RNA-seq data.
#' Such tables can be produced by a variety of packages, e.g. DESeq2, edgeR, etc.
#' However, the column names should be redefined since each of this packages uses
#' different names for similar metrics.
#' The required columns are c("gene.id", mean", "log2FC", "pvalue", "padj").
#' Additional columns can be provided but will be ignored for the analysis.
#'
#' @param table.name  name for the data table, which will be displayed in the plots.
#'
#' @param sort.column="none" Column to sort the result.
#' Supported: c("none", "mean", "log2FC", "pvalue", "padj").
#'
#' @param thresholds=c(padj = 0.05, FC = 2)
#'
#' Thresholds on some specific scores, used for display purpose, and to
#' add some columns to the result table, indicating if the gene passes the
#' threshold or not. Supported threshold fields:
#' \itemize{
#' \itme FC (orientation-insensitive fold change)
#' \item padj (adjusted p-value)
#' \item pvalue (nominal p-value)
#' \item evalue (expected number of false positives)
#' \item log2FC (absolute value of the log2 fold change)
#' }
#'
#' @param round.digits=3 Significant digits to round the values of the
#' output table. Set to NA to avoid rounding.
#'
#' @param dir.figures=NULL if not NULL, figures will be saved in the specified directory.
#' @param verbose=0 level of verbosity
#'
#' @examples
#'  ##-----------------------------------
#'  ## Post-process a DESeq2 result table
#'  ##-----------------------------------
#'  deseq2.result.table <- data.frame(
#'     "gene.id" = row.names(deseq2.res),
#'     "mean" = deseq2.res$baseMean,
#'     "log2FC" = deseq2.res$log2FoldChange,
#'     "pvalue" = deseq2.res$pvalue,
#'     "padj" = deseq2.res$padj)
#' deseq2.result.table <- DEGtablePostprocessing(deseq2.result.table,
#'     table.name="DESeq2", sort.column = "padj")
#'
#' @export
DEGtablePostprocessing <- function(deg.table,
                               table.name="DEG_table",
                               sort.column = "none",
                               thresholds = c(padj = 0.05, FC = 2),
                               round.digits = 3,
                               dir.figures = NULL,
                               verbose = 0) {

  # names(deg.table)
  if (verbose >= 1) {
    message("\tDEGtablePostprocessing()\tTable name: ", table.name, "\t", nrow(deg.table), " features (rows)")
  }

  thresholds <- unlist(thresholds)

  col.descriptions <- vector() ## Initialize vector with column descriptions

  # Check that the input table contains the required columns
  required.columns <- c("gene.id", "mean", "log2FC", "pvalue", "padj")
  for (column in required.columns) {
    if (!column %in% names(deg.table)) {
      stop(paste("DEG table should contain",column,"column"))
    }
  }

  ## Make sure that row names correspond to gene IDs
  row.names(deg.table) <- deg.table$gene.id

  ## Sort DEG table if required
  sort.decreasing <- c("mean" = TRUE,
                       "padj" = FALSE,
                       "pvalue" = FALSE,
                       "log2FC" = TRUE)
  if (sort.column != "none") {
    if (verbose >= 2) { message("\t\tSorting DEG table by ", sort.column) }
    deg.table <- deg.table[order(deg.table[,sort.column],
                                 decreasing = sort.decreasing[sort.column]),]
  }
  # head(deg.table)
  # tail(deg.table)
  # summary(deg.table)

  ## Compute fold-change from log2 fold change
  ## Beware: for the fold-change we ake the absolute value of log2FC,
  ## to have a fold change irrespective of the up- or -down sense
  deg.table$FC <- 2^abs(deg.table$log2FC)
  # hist(deg.table$log2FC, breaks=1000)
  # hist(2^(deg.table$log2FC), breaks = 1000)
  # hist(deg.table$FC, breaks=1000)

  col.descriptions["FC"] <- "Orientation-insensitive fold change (always >= 1). "

  ## Compute E-value
  deg.table$evalue <- deg.table$pvalue * nrow(deg.table)
  col.descriptions["evalue"] <- "Expected nb of FP (=pvalue * number of tests)"

  ## Compute the rank of p-value sorted genes
  deg.table$padj.rank <- rank(deg.table$padj)
  col.descriptions["pval.rank"] <- "Rank of genes sorted by increasing adjusted p-values"

  ## Indicate the sense of the regulation (up-down)
  deg.table$sign <- sign(deg.table$log2FC)
  col.descriptions["sign"] <- "Sign of the regulation (1= up, -1 = down)"

  # as.data.frame(col.descriptions)

  ## Label the genes passing the FDR, E-value and fold-change thresholds
  threshold.type <- c(
    "pvalue" = "le",
    "padj" = "le",
    "evalue" = "le",
    "log2FC" = "abs.ge",
    "FC" = "ge")
  thresholds.to.apply <- intersect(names(thresholds), names(threshold.type))

  if (verbose >= 1) {
    message("\t\tApplying thresholds: ", paste(collapse = ", ", thresholds.to.apply),
            "\tStarting features\t", nrow(deg.table))
  }
  selection.columns <- paste(sep = "", thresholds.to.apply, "_", thresholds[thresholds.to.apply])
  names(selection.columns) <- thresholds.to.apply
  selected.features <- rep(TRUE, length.out = nrow(deg.table))
  s <- thresholds.to.apply[1]
  for (s in thresholds.to.apply) {
    if (threshold.type[s] == "le") {
      threshold.passed <- !is.na(deg.table[, s]) & deg.table[, s] <= thresholds[s]
    } else if (threshold.type[s] == "ge") {
      threshold.passed <- !is.na(deg.table[, s]) & deg.table[, s] >= thresholds[s]
    } else if (threshold.type[s] == "abs.le") {
      threshold.passed <- !is.na(deg.table[, s]) & abs(deg.table[, s]) <= thresholds[s]
    } else if (threshold.type[s] == "abs.ge") {
      threshold.passed <- !is.na(deg.table[, s]) & abs(deg.table[, s]) >= thresholds[s]
    } else if (threshold.type[s] == "abs.log2.ge") {
      threshold.passed <- !is.na(deg.table[, s]) & abs(log2(deg.table[, s])) <= thresholds[s]
    } else if (threshold.type[s] == "abs.log2.le") {
      threshold.passed <- !is.na(deg.table[, s]) & abs(log2(deg.table[, s])) >= thresholds[s]
    } else {
      stop("Invalid threshold criterion: ", s, ". Supported: ", paste(collapse = ", ", names(threshold.type)))
    }
    # summary(threshold.passed)
    selected.features <- selected.features & threshold.passed

    deg.table[, selection.columns[s]] <- threshold.passed * 1
    col.descriptions[selection.columns[s]] <- paste("Passing", threshold.type[s], "threshold on", s)
    if (verbose >= 1) {
      message("\t\t\t", threshold.type[s],
              " threshold on ", s, ": ", thresholds[s],
              "\tPassing features: ", sum(threshold.passed),
              "\tKept features: ", sum(selected.features))
    }
  }

  ## Select genes passing all thresholds
  if (verbose  >= 3) {
    message("\t\tSelection columns: ", paste(collapse = ", ", selection.columns))
  }
  deg.table[,"DEG"] <- 1*(selected.features)
  # table(deg.table[,"DEG"])
  if (verbose >= 3) {
    print(data.frame(col.descriptions))
  }

  ## Round columns to a reasonable number of significant digits
  if (!is.na(round.digits)) {
    verbose(paste("Rounding values to", round.digits, "digits"), 2)
    for (col in setdiff(names(deg.table), selection.columns)) {
      if (is.numeric(deg.table[,col])) {
        verbose(paste("Rounding", col), 2)
        deg.table[,col] <- signif(digits = 3, deg.table[,col])
      }
    }
  }

  ## Export figures
  if (!is.null(dir.figures)) {
    verbose(paste("\tSaving figures in directory", dir.figures), 2)
    dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)

    ## Draw Venn diagram with number of genes declared significant
    ## according to the selection criteria (threshold fields).
    selection.venn.counts <- vennCounts(deg.table[,selection.columns])
    venn.file <- file.path(dir.figures, paste(sep = "", table.name, "selection_Venn.pdf"))
    if (verbose >= 1) { message("\t\tVenn diagram\t", paste(collapse = ", ", selection.columns), "\t", venn.file) }
    pdf(file = venn.file)
    limma::vennDiagram(selection.venn.counts, cex = 1, main = paste(table.name, "selected genes"),
                       circle.col = c("orange", "blue"), mar = c(0,0,5,0))
    silence <- dev.off()

    ## Histogram of the nominal p-values
    pval.hist.file <- file.path(dir.figures, paste(sep = "", table.name, "_pval_hist.pdf"))
    if (verbose >= 1) { message("\t\tP-value histogram\t", pval.hist.file) }
    pdf(file = pval.hist.file, width = 7, height = 5)
    # hist(deg.table$pvalue, breaks = seq(from = 0, to = 1, by = 0.05),
    #      xlab = "Nominal p-value",
    #      ylab = "Number of genes",
    #      main = paste(table.name, "P value distribution"),
    #      col = "#BBBBBB")
    multTestCorr <- multipleTestingCorrections(deg.table$pvalue)
    PlotPvalDistrib.MultiTestTable(multitest.result = multTestCorr,
                                   draw.m0.line = TRUE, draw.mean.line = FALSE, draw.lambda = TRUE,
                                   main = paste(table.name, "\nP value distribution"), legend.cex = 0.9)
    silence <- dev.off()

    ## Histogram of the log2FC
    log2FC.file <- file.path(dir.figures, paste(sep = "", table.name, "_log2FC_hist.pdf"))
    if (verbose >= 1) { message("\t\tLog2-FC histogram\t", log2FC.file) }
    pdf(file = log2FC.file, width = 7, height = 5)
    hist(deg.table$log2FC, breaks = 100,
         xlab = "log2(fold change)",
         ylab = "Number of genes",
         main = paste(table.name, "\nlog2FC distribution"),
         col = "#BBBBBB")
    silence <- dev.off()

    ## Volcano plot
    message("thresholds['FC']\t", thresholds["FC"])
    message("log2(thresholds['FC'])\t", log2(thresholds["FC"]))
    nona.deg <- na.omit(deg.table)
    volcano.file <- file.path(dir.figures, paste(sep = "", table.name, "_volcano_plot_padj.pdf"))
    if (verbose >= 1) { message("\t\tVolcano plot\t", volcano.file) }
    pdf(file = volcano.file, width = 7, height = 7)
    VolcanoPlot.MultiTestTable(
      multitest.table = nona.deg,
      effect.size.col = "log2FC",
      effect.threshold = log2(thresholds["FC"]),
      control.type = "padj",
      alpha = thresholds["padj"],
      main = paste(table.name, "\nVolcano plot"))
    silence <- dev.off()


  }
  return(deg.table)
}
