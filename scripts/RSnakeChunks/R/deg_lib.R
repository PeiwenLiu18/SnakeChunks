#' @title check required libraries for DEG analysis, and install them if required
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @param required.libraries a vector contianing the names of required CRAN libraries, which will be installed with install.packages()
#' @param required.bioconductor a vector containing the required BioConductor libraries, which will be installed with biocLite
CheckRequiredLibraries <- function(required.libraries,
                                   required.bioconductor = NULL) {
  for (lib in required.libraries) {
    if (!require(lib, character.only = TRUE)) {
      install.packages(lib)
      library(lib, character.only = TRUE)
    }
  }

  for (lib in required.bioconductor) {
    if (!require(lib, character.only = TRUE)) {
      ## try http:// if https:// URLs are not supported
      source("https://bioconductor.org/biocLite.R")
      biocLite(lib)
    }
    if (!require(lib, character.only = TRUE)) {
      stop("Missing library: ", lib, " could not be installed")
    }
  }

}

#' @title Load parameters and required libraries
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
LoadDEGparam <- function(yamlFile) {
  library(yaml)
  data <- yaml.load_file(yamlFile)
}


#' @title Display messages at a given verbosity level
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Display messages depending on user-defined verbosity level
#'
#' @details
#' First version: 2015-03.
#' Last modification: 2015-03.
#'
#' @param verbosity   Level of verbosity above which the message should be printed.
#' @param print.date=TRUE   Print date and time
#'
#' @examples
#'
#' verbosity <- 1 ## Define level of verbosity
#'
#' ## This message will be printed because the level is <= verbosity
#' verbose("This is printed", 1)
#'
#' ## This message will not be printed because the verbosity is inferior to the specified level
#' verbose("This is not printed", 2)
#'
#' @export
verbose <- function(message.content,
                    level = 1,
                    print.date = TRUE) {
  if (!exists("verbosity")) {
    verbosity <- 1
  }
  if (verbosity >= level) {
    if (print.date) {
      message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\t", message.content)
    } else {
      message(message.content)
    }
  }
}

#' @title Draw a volcano plot for RNA-seq data.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Draw an equivalent of the t-test volcano plot for RNA-seq
#' data. The main idea is to display simultaneously the effect size (X axis)
#' and the significance (Y axis). The measure of the effect size is in this
#' case the log2 fold change. (logFC).
#'
#' There is a close relationship between this RNA-seq version of the volcano
#' and the volcano plots commonly used to show the result of a t-test with microarrays.
#'
#' @details

########## Draw boxplots with read counts per genes for each sample ################
count.boxplot <- function(count.table,
                          sample.desc,
                          sample.label.col=1,
                          xlab="Raw counts",
                          main="Box plots per sample: raw counts",
                          plot.file=NULL) {

  ## Adapt boxplot size to the number of samples and label sizes
  boxplot.lmargin <- max(nchar(sample.desc$label))/3+5
  boxplot.height <- length(sample.ids)/3+2

  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=8, height=boxplot.height)
  }

  par(mar=c(5,boxplot.lmargin,4,1)) ## adapt axes
  boxplot(count.table, horizontal=TRUE, col=sample.desc$color,
          xlab=xlab, names=sample.desc[, sample.label.col],
          main=main, las=1)
  grid(col="grey", lty="solid",ny = 0)
  if (!is.null(plot.file)) {
    silence <- dev.off()
  }
}

################################################################
## Draw a heatmap with the inter-sample correlation matrix.
count.correl.heatmap <- function(count.table,
                                 main = NULL,
                                 plot.file = NULL,
                                 score = "cor", ## Supported: PCC (Pearson Correlation Coefficient) or SERE (handled by SARTools)
                                 cor.method = "pearson", ## Passed to cor()
                                 log.transform = TRUE, # Perform a log transformation of the values before plotting. Only done for PCC since SERE requires raw counts.
                                 epsilon = 0.1, # Add an epsilon to zero values before log transformation, in order to -Inf values
                                 zlim = NULL,
                                 levels = 100,
                                 gray.palette = TRUE,
                                 gamma = 1, # Gamma parameter passed to gray.colors()
                                 plot.values = TRUE,
                                 ...
) {

  ## Suppress rows with NA values
  count.table <- na.omit(count.table)
  # range(count.table)

  ## Adapt margins to the number of samples and label sizes
  sample.names <- names(count.table)
  if (is.null(sample.names)) {
    margin <- 5
  } else {
    margin <- max(nchar(sample.names))/3+5
  }

  ## Compute comparison score
  if (score == "SERE") {
    # lane.totals <- colSums(count.table)
    # SERE.table <- SERE_fun(
    #   observed = as.matrix(count.table),
    #   laneTotals = lane.totals
    # )
    library(SARTools)
    SERE.table <- tabSERE(as.matrix(count.table))
    SERE.table <- SERE.table/100 ## SARTools returns scores multiplied by 100
    count.cor <- max(SERE.table) - SERE.table ## Invert SERE to get a similarity rather than dissimilarity score
    # count.cor <- SERE.table
  } else if (score == "cor") {
    if (log.transform) {
      # count.table[count.table == 0] <- epsilon ## add epsilon to
      count.table <- log2(count.table + epsilon)
    }
    count.cor <- as.matrix(cor(count.table, method=cor.method))
  } else {
    stop("count.correl.heatmap(): ",
         score,
         " is not a valid score. Supported: cor, SERE. ")
  }
  # count.range <- range(count.table)

  ## Define a color palette for heatmaps.
  if (gray.palette) {
    ## Use a grayscale color
    cols.heatmap <- gray.colors(n = levels, start = 1, end = 0, gamma = gamma)
  } else {
    ## I like this Red-Blue palette because
    ## - it suggests a subjective feeling of warm (high correlation)/cold (low correlation)
    ## - it can be seen by people suffering from red–green color blindness.
    cols.heatmap <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(levels))
  }

  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=8, height=boxplot.height)
  }

  # Define main title
  if (is.null(main)) {
    if (score == "cor") {
      main <- paste(cor.method, " correlation")
    } else if (score == "SERE") {
      main <- paste(score, "scor")
    } else {
      main <- paste(score)
    }
  }

  if (is.null(zlim)) {
    zlim <- range(count.cor)
  }
  breaks <- seq(from=zlim[1], to = zlim[2], length.out = levels + 1)

  if (plot.values) {
    hm <- heatmap.2(count.cor,  scale="none", trace="none",
                    #breaks=c(-1, seq(0,1,length.out = 100)),
                    main=main, margins=c(margin,margin),
                    col=cols.heatmap, breaks = breaks,
                    cellnote = signif(digits=2, count.cor),
    ...)
  } else {
    hm <- heatmap.2(count.cor,  scale="none", trace="none",
                    #breaks=c(-1, seq(0,1,length.out = 100)),
                    main=main, margins=c(margin,margin),
                    col=cols.heatmap, breaks = breaks,
                    ...)
  }


  if (!is.null(plot.file)) {
    silence <- dev.off()
  }

  return(count.cor)
}

# ## Draw a heatmap with the inter-sample correlation matrix. ###########
# count.correl.heatmap <- function(count.table,
#                                  main="Correlation between raw counts",
#                                  plot.file=NULL,
#                                  log.transform=TRUE, # Perform a log transformation of the values before plotting
#                                  epsilon=0.1, # Add an epsilon to zero values before log transformation, in order to -Inf values
#                                  zlim = NULL,
#                                  grey = FALSE,
#                                  plot.values = TRUE, # plot correlation values on the heatmap
#                                  ...
#                                  ) {
#
#
#   ## Adapt boxplot size to the number of samples and label sizes
#   margin <- max(nchar(names(count.table)))/3+5
#
#   if (log.transform) {
#     range(count.table)
#     count.table[count.table==0] <- epsilon
#     count.table <- log10(count.table)
#   }
#   count.cor <- as.matrix(cor(count.table))
#
#   ## Limits for the color scale
#   if (is.null(zlim)) {
#     zlim <- range(count.cor)
#   }
#
#   ## Use a grayscale color
#   if (grey) {
#     cols.heatmap <- gray.colors(256, start = 1, end = 0, gamma = 3, alpha = NULL)
#   } else {
#     ## Define a color palette for heatmaps. I like this Red-Blue palette because
#     ## - it suggests a subjective feeling of warm (high correlation)/cold (low correlation)
#     ## - it can be seen by people suffering from red–green color blindness.
#     cols.heatmap <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))
#
#   }
#
#   ## Sample-wise library sizes
#   if (!is.null(plot.file)) {
#     message("Generating plot", plot.file)
#     pdf(file=plot.file, width=8, height=boxplot.height)
#   }
#
#
#   if (plot.values) {
#     cellnote <- signif(digits=2, count.cor)
#   } else {
#     cellnote <- FALSE
#   }
#   hm <- heatmap.2(count.cor,
#                   scale="none",
#                   trace="none",
#                   zlim = ,
#                   #breaks=c(-1, seq(0,1,length.out = 100)),
#                   main=main,
#                   margins=c(margin,margin),
#                   col=cols.heatmap)
#
# #                  zlim = zlim,
#                   cellnote = cellnote,
#                   ...
#                   )
#
#
#   if (!is.null(plot.file)) {
#     silence <- dev.off()
#   }
#
#   return(count.cor)
# }

## Generate a set of plots displaying some sample-wise statistics
sample.description.plots <- function (sample.desc,
                                      stats.per.sample,
                                      dir.figures,
                                      exploratory.plots=FALSE) {

  plot.files <- c() ## To retun: list of files with plots

  par.ori <- par() ## Save original plot parameters


  ## Library size barplots
  plot.files["Mreads_barplot"] <- file.path(dir.figures, "sample_libsum_barplot.pdf")
  libsize.barplot(stats.per.sample, plot.files["Mreads_barplot"])

  ################################################################
  ## Boxplots of raw counts and CPMs, in linear + log axes.
  ## These plots give a pretty good intuition of the raw data per sample:
  ## library sizes, outliers, dispersion of gene counts.

  ## Boxplot of raw counts
  plot.files["sample_boxplot_counts"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_counts.pdf"))
  count.boxplot(count.table, stats.per.sample,
                xlab="Raw counts", main="Box plots per sample: raw counts",
                plot.file=plot.files["sample_boxplot_counts"])

  ## Boxplot of log10-transformed counts
  plot.files["sample_boxplot_counts_log10"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_counts_log10.pdf"))
  count.boxplot(count.table.log10, stats.per.sample,
                xlab="log10(counts)", main="Box plots per sample: log10(counts)",
                plot.file=plot.files["sample_boxplot_counts_log10"])

  ## Boxplot of CPMs
  plot.files["sample_boxplot_CPM"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_CPM.pdf"))
  count.boxplot(cpms, stats.per.sample,
                xlab="CPM", main="Box plots per sample: counts per million reads (CPM)",
                plot.file=plot.files["sample_boxplot_CPM"])

  ## Boxplot of log10-transformed CPMs
  plot.files["sample_boxplot_CPM_log10"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_CPM_log10.pdf"))
  count.boxplot(cpms.log10, stats.per.sample,
                xlab="log10(CPM)", main="Box plots per sample: counts per million reads (CPM)",
                plot.file=plot.files["sample_boxplot_CPM"])

  par <- par.ori ## Restore original plot parameters
  par(mar=c(4.1,5.1,4.1,1.1))

  ## Draw sample correlation heatmaps for the raw read counts
  plot.files["sample_correl_heatmap_counts"] <- paste(sep = "", prefix["general.file"],"_sample_correl_heatmap_counts.pdf")
  count.correl.heatmap(count.table, plot.file=plot.files["sample_correl_heatmap_counts"])
#   hm <- heatmap.2(,  scale="none", trace="none",
#                   main="Correlation between raw counts", margins=c(8,8),
#                   col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))

  ## Draw sample correlation heatmaps for CPM. Actually it gives exactly the
  ## same result as correlation between raw counts, since the correlation has a
  ## standardizing effect.
  # pdf(file=paste(sep = "", prefix["general.file"],"_sample_correl_heatmap_cpms.pdf"))
  # hm <- heatmap.2(as.matrix(cor(cpms)),  scale="none", trace="none",
  #                 main="Correlation between CPM",
  #                 col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))
  # quiet <- dev.off()

  ## Plot the first versus second components of samples
  cpms.pc <- prcomp(t(cpms))
  plot.file <- paste(sep = "", prefix["general.file"],"_CPM_PC1-PC2.pdf")
  plot.files["CPM_PC1-PC2"] <- plot.file
  message("Generating plot", plot.file)
  pdf(file=plot.file)
  plot(cpms.pc$x[,1:2], panel.first=grid(), type="n", main="First components from PCA-transformed CPMs")
  text(cpms.pc$x[,1:2], labels = sample.conditions, col=sample.desc$color)
  quiet <- dev.off()



  ## Exploratory plots, should not be done for all projects.
  if (exploratory.plots) {
    verbose("Drawing generic plots from the whole count table", 1)

    ## Plot the impact of the normalization factor (library sum , median or percentile 75)
    plot.file <- file.path(dir.diffexpr, paste(sep = "", "CPM_libsum_vs_median_vs_perc75.png"))
    plot.files["CPM_libsum_vs_median_vs_perc75"] <- plot.file
    message("Generating plot", plot.file)
    png(file= plot.file, width=1000, height=1000)
    cols.counts <- as.data.frame(matrix(sample.desc$color, nrow=nrow(count.table), ncol=ncol(count.table), byrow = TRUE))
    colnames(cols.counts) <- names(count.table)
    rownames(cols.counts) <- rownames(count.table)
    plot(data.frame("libsum" = as.vector(as.matrix(cpms.libsum)),
                    "median" = as.vector(as.matrix(cpms.median)),
                    "perc75" = as.vector(as.matrix(cpms.perc75))),
         col=as.vector(as.matrix(cols.counts)))
    quiet <- dev.off()

    ## Plot some sample-wise statistics
    plot.file <- file.path(dir.diffexpr, paste(sep = "", "sample_statistics_plots.pdf"))
    plot.files["sample_statistics_plots"] <- plot.file
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=10, height=10)
    par(mar=c(5,5,1,1)) ## adpt axes
    par(mfrow=c(2,2))
    ## Median versus mean
    plot(stats.per.sample[,c("mean", "median")],
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)

    ## First versus third quartile
    plot(stats.per.sample[,c("perc25", "perc75")],
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)

    ## Sum versus third quartile.
    plot(stats.per.sample[,c("sum", "perc75")],
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)

    ## Mean versus third quartile.
    plot(stats.per.sample[,c("mean", "perc75")],
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)
    par(mfrow=c(1,1))
    quiet <- dev.off()
  }

  return(plot.files)
}

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
  deseq2.result.table <- complete.deg.table(
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
edger.analysis <- function(counts,
                           condition,
                           ref.condition=NULL,
                           comparison.prefix,
                           title = comparison.prefix,
                           dir.figures=NULL,
                           norm.method = "TMM",
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
  edger.result.table <- complete.deg.table(
    deg.table = edger.result.table,
    table.name = paste(sep = "_", "edgeR", prefix["comparison"]),
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



