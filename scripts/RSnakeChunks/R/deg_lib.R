#' @title check required libraries for DEG analysis, and install them if required
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @param required.libraries a vector contianing the names of required CRAN libraries, which will be installed with install.packages()
#' @param required.bioconductor a vector containing the required BioConductor libraries, which will be installed with biocLite
#' @param verbose=1 verbosity
#' @export
CheckRequiredLibraries <- function(required.libraries,
                                   required.bioconductor = NULL,
                                   verbose = 1) {
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
#' @export
LoadDEGparam <- function(yamlFile) {
  library(yaml)
  data <- yaml.load_file(yamlFile)
  return(data)
}


#' @title Display messages at a given verbosity level
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Display messages depending on user-defined verbosity level
#'
#' @param verbosity   Level of verbosity above which the message should be printed.
#' @param print.date=TRUE   Print date and time
#'
#' @examples
#'
#' verbosity <- 1 ## Define level of verbosity
#'
#' ## This message will be printed because the level is <= verbosity
#' if (verbose >= 2) { message("This is printed")
#'
#' ## This message will not be printed because the verbosity is inferior to the specified level
#' if (verbose >= 2) { message("This is not printed")
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

#' @title Draw a heatmap with the inter-sample correlation matrix.
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
#' @export
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
    if (verbose >= 2) { message("Drawing generic plots from the whole count table") }

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



