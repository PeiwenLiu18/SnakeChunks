#' @title Normalisation and differential analysis of an RNA-seq count table
#' @author Jacques van Helden
#' @description Workflow for the analysis of differential expression with RNA-seq transcripttome data.
#' Includes
#' \itemize{
#' \item sample-wiste descriptive statistics (central tendency, dispersion, percentiles, ...)
#' \item between-sample normalisation with different options
#' \item detection of differentially expressed genes with DESeq2 and edgeR
#' \item comparisons between results obtained with alternative methods for between-sample normalisation
#' \item export of all results in tab-separated value files
#' \item export of figures in pdf and png formats
#' \item automatic generation of an HTML report (via Rmd) with links to result files and figures
#' }
#' @param countFile file containing a table with counts of reads per feature (row) in each sample (column).
#' @param configFile a yaml-formatted file defining the mandatory + some optional parameteres
#' @param main.dir=getwd()  directory from which the script runs. Paths are defined relative to this directory.
#' @param result.dir="results" directory where the results will be stored. Should be defined relative to main directory.
#' @param rmd.report="rnaseq_deg.Rmd" path to the R markdown report, which is also converted to html and pdf reports
#' @param verbose=1 level of verbosity
#'
#' @export
RNAseqAnalysis <- function(countFile,
                           configFile,
                           checkSampleIDs=TRUE,
                           main.dir = getwd(),
                           result.dir = "results",
                           rmd.report = "rnaseq_deg.Rmd",
                           verbose = 1) {

  

  ## ---- TO DO ----
  ##
  ## Compute a table with the mean coutns per condition after normalisation
  ## Generate scatter plots comparing mean counts per condition
  ##
  ## Check non-used variables: suppress or use them
  ## - cols.heatmap
  ## - DEG.selection.criterion
  ## - filteredCounts.log10
  ## - filteredCounts.log2
  ## - scaledCounts.libsum
  ## - scaledCounts.perc95
  ## - scaledCounts.median
  ## - current.labels
  ##
  ## - Add PCA
  ## - Add multi-group differential analysis
  ## - Comparison between different anlayses of the design file
  ## - Clustering


  ## ---- Load required libraries ----
  required.cran.libraries <- c("knitr",
                               "yaml",
                               "pander",
                               # "xlsx",
                               "ascii",
                               "xtable",
                               "gplots",
                               "codetools",
                               "RColorBrewer",
                               'rmarkdown',
                               "devtools"#,
                               #                        "stats4bioinfo" ## Generic library from Jacques van Helden
  )
  LoadRequiredCRANPackages(required.cran.libraries, verbose = 1)

  required.bioconductor <- c(
    "edgeR",
    "DESeq2",
    "limma",
    #  "SARTools", ## for SERE coefficient
    "GenomicFeatures")
  LoadRequiredBioconductorPackages(required.bioconductor, verbose = 1)
  ## ---- knitr_setup, include=FALSE,  eval=TRUE, echo=FALSE, warning=FALSE----

  ## To do: check if these options are taken into account at rendering
  knitr::opts_chunk$set(
    fig.path = "figures/",
    echo = FALSE,
    eval = TRUE,
    cache = FALSE,
    message = FALSE,
    warning = FALSE)


  ## ---- Initialize lists of input and output files ----
  ## This will later serve to generate the report.
  message("\tInitializing indexes and main parameters")
  dirs <- vector()
  dirs["main"] <- main.dir
  if (!dir.exists(main.dir)) {
    stop("Cannot set working directory to\t", main.dir)
  }
  setwd(main.dir)
  message("\tMain directory: ", main.dir)


  ## Directory to store differential expression results
  dirs["output"] <- result.dir ## Index dir for the report
  message("\tOutput directory: ", result.dir)
  dir.create(result.dir, showWarnings = FALSE, recursive = TRUE)

  ## Input files
  infiles <- vector()   ## Input files

  ## Configuration file
  infiles["config"] <- configFile
  message("\tConfiguration file: ", configFile)

  ## Count table
  infiles["count_table"] <- countFile
  message("\tCount table file: ", countFile)



  ## Output files
  outfiles <- vector()  ## For tab-separated value files


  ## ---- Load RNA-seq data, metadata and configuration ----
  dataset <- LoadRNAseqDataset(countFile, configFile, verbose = verbose)
  parameters <- dataset$parameters
  
  ## Make the dataset fields available as variables (to avoid prepending "dataset" everywhere)
  attach(dataset)
  # detach(dataset)
  
  
  # dim(rawCounts)
  # dim(filteredCounts)
  
  ## ---- log2 transformation for raw and filtered counts ----
  rawCounts.log2 <- log2(rawCounts + parameters$DEG$epsilon)
  filteredCounts.log2 <- log2(filteredCounts + parameters$DEG$epsilon)
  
  
  ## Instantiate some local variables for some frequently used parameters,
  ## for the readability.
  norm.methods <- parameters$DEG$norm_method
  norm.percentile <- parameters$DEG$norm_percentile
  thresholds <- as.vector(parameters$DEG$thresholds)
  epsilon <- parameters$DEG$epsilon

  ## Index input files
  infiles["sample descriptions"] <- parameters$metadata$samples
  infiles["design"] <- parameters$metadata$design
  if ((!is.null(parameters$DEG$blacklist)) & (parameters$DEG$blacklist != "")) {
    infiles["black_list"] <- parameters$DEG$blacklist
  }
  
  ## ---- Define some output parameters -----------------------------------------------------

  ## Compute count prefix, i.e. the basename to export various transformations of the count tables
  if (is.null(parameters$dir$count_prefix)) {
    ## use the basename of the count table as prefix for output file
    parameters$dir$count.prefix <- basename(path = countFile)

    ## Suppress usual extensions from the basename
    ext <- ".tsv"
    for (ext in c(".tsv", ".tab", ".txt", ".csv")) {
      parameters$dir$count_prefix <- sub(
        pattern = ext, replacement = "", x = parameters$dir$count.prefix)
    }
  }
  count.prefix <- parameters$dir$count.prefix
  message("\tFile prefix for normalized count tables ", count.prefix)


  ## Directory for Figures
  dir.figures.samples <- file.path(result.dir, "figures")
  dirs["sample_figures"] <- dir.figures.samples
  message("\tDirectory for the sample-related figures: ", dir.figures.samples)
  dir.create(dir.figures.samples, showWarnings = FALSE, recursive = TRUE)

  ## Create output file lists for the different figure formats
  figure.formats <- parameters$figure_formats
  figure.files <- list()
  for (fig.format in figure.formats) {
    figure.files[[fig.format]] <- vector()
  }


  ## Directory to export result files in tba-separated value (tsv) format
  dir.tables.samples <- file.path(result.dir, "samples/tables")
  dirs["tsv"] <- dir.tables.samples ## Index directory for the report
  message("\tDirectory to export sample-related tables:\t", dir.tables.samples)
  dir.create(dir.tables.samples, showWarnings = FALSE, recursive = TRUE)

  ## ---- Set some default parameters ----

  ## Color palette for heatmaps. I like this Red-Blue palette because
  ## - it suggests a subjective feeling of warm (high correlation)/cold (low correlation)
  ## - it can be seen by people suffering from red–green color blindness.
  if (!exists("cols.heatmap")) {
    cols.heatmap <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))
  }



  ## ---- Initialize the Rmd report (index of input/output file) ----
  #  rmd.report <- "rmd.report"
  dir.report <- dirname(rmd.report)
  dirs["report"] <- dir.report ## Index directory for the report
  dir.create(dir.report, showWarnings = FALSE, recursive = TRUE)
  
  report.socket <- file(rmd.report)


  Rmd.header <- paste(
    sep = '',
    '---
title:  "', parameters$title,'"
author: "', parameters$author, ' <', parameters$author_email, '>"
date: Last update:`r format(Sys.time())`
output:
  html_document:
    fig_caption: yes
    highlight: kate
    self_contained: yes
    theme: cerulean
    toc: yes
    toc_depth: 4
    toc_float: yes
  pdf_document:
    fig_caption: yes
    highlight: kate
    toc: yes
    toc_depth: 3
  word_document:
    toc: yes
    toc_depth: 2
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
    echo = FALSE,
    eval = TRUE,
    cache = TRUE,
    message = FALSE,
    warning = FALSE,
    comment = "",
    fig.align = "center",
    fig.path = "figures/")
```
')

  report.text <- Rmd.header

  ## Print the description
  if (!is.null(parameters$description)) {
    report.text <- append(report.text, "\n\n## Description\n")
    report.text <- append(report.text, parameters$description)
  }
  
  ## Print the data source
  if (!is.null(parameters$dataset)) {
    report.text <- append(report.text, paste(sep = "","\n\n- Dataset: ", parameters$dataset))
  }
  
  ## Print the original publication
  if (!is.null(parameters$citation)) {
    report.text <- append(report.text, paste(sep = "","\n- Publication: ", parameters$citation))
  }
  
  ## Print the threshold tables
  ## NOTE: I should evaluate what I do with the kable calls
  report.text <- append(report.text, "\n\n## Parameters\n")

  report.text <- append(report.text, 
                       paste(sep = "", "- Congiguration file: ", 
                             ReportLink(source = rmd.report, target = configFile)))
                       
  report.text <- append(report.text, "\n\n### Thresholds\n")
  report.text <- append(
    report.text,
    kable(t(as.data.frame(thresholds)), col.names = "Threshold",
          caption = "Thresholds for the selection of differentially expressed genes. "))

  outfiles["thresholds"] <- file.path(dir.tables.samples, "thresholds.tsv")
  message("\t\tExporting threshold table\t", outfiles["thresholds"])
  write.table(x = t(as.data.frame(thresholds)),
              file = outfiles["thresholds"],
              sep = "\t", row.names = TRUE, col.names = FALSE)
  # list.files(dir.tables.samples)
  # system(paste("open", dir.tables.samples))



  ## Print the sample descriptons
  report.text <- append(report.text, "\n\n## Sample descriptions\n")
  report.text <- append(report.text, kable(sample.desc, caption = "Sample description table"))


  samples.per.condition <- as.data.frame.table(table(sample.conditions))
  report.text <- append(report.text, "\n\n### Samples per condition\n")
  report.text <- append(report.text, kable(samples.per.condition, caption = "Samples per condition",
                                         col.names = c("Condition", "Nb samples")))
  if (verbose >= 2) { print(as.data.frame(samples.per.condition)) }


  ## Print out the design table (pairs of conditions to be compared)
  report.text <- append(report.text, "\n\n## Design\n")
  report.text <- append(report.text, kable(
    comparison.summary,
    row.names = TRUE,
    caption = "**Design**. Each row describes one comparison between two conditions."))


  # ## Check that the header of rawCounts match the sample IDs
  # ids.not.found <- setdiff(sample.ids, names(rawCounts)) ## Identify sample IDs with no column in the count table
  # if (length(ids.not.found) == length(sample.ids)) {
  #   colnames(rawCounts) <- sample.ids
  #   ids.not.found <- setdiff(sample.ids, names(rawCounts)) ## Identify
  # } else if (length(ids.not.found) > 0) {
  #   stop(length(ids.not.found), " missing columns in count table\t", countFile,
  #        "\n\tMissing columns: ", paste(collapse = "; ", ids.not.found))
  # }
  # 
  # ## Restrict the count table to the sample IDs found in the sample description file
  # rawCounts <- rawCounts[, sample.ids]
  # # names(rawCounts)
  # # dim(rawCounts)
  # # (zeros <- sum(rawCounts == 0))
  # # (zero.rows <- sum(apply(rawCounts, 1, sum, na.rm=TRUE) == 0))
  # # zeros/ zero.rows


  ## ---- Compute sample-wise statistics on counts -----

  ## Selected sample statistics
  selected.stats <- c("Mcounts",
                      "sum",
                      "min",
                      "zeros",
                      "non.null",
                      "perc05",
                      "Q1",
                      "mean",
                      "median",
                      "Q3",
                      "perc95",
                      "max",
                      "max.sum.ratio",
                      "mean.median.ratio",
                      "fract.below.mean")

  message("\tComputing sample-wise statistics on all counts (non-filtered)")
  #stats.per.sample <- calc.stats.per.sample(sample.desc, rawCounts)
  # View(stats.per.sample.all)
  # dim(rawCounts)
  # dim(sample.desc)
  stats.per.sample.all <- cbind(
    sample.desc,
    ColStats(x = rawCounts, verbose = verbose, selected.stats = selected.stats))
  #stats.per.sample.all$Mcounts <- stats.per.sample.all$sum / 1e6
  # View(stats.per.sample.all.all)
  outfiles["stats_per_sample_all_features"] <- file.path(dir.tables.samples, "stats_per_sample_all_features.tsv")
  message("\t\t", outfiles["stats_per_sample_all_features"])
  write.table(x = stats.per.sample.all, file = outfiles["stats_per_sample_all_features"], quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

  ## Compute statistics ommitting zero values
  message("\tComputing sample-wise statistics for non-zero counts")
  rawCounts.nozero <- rawCounts
  rawCounts.nozero[rawCounts.nozero == 0] <- NA
  stats.per.sample.nozero <- cbind(
    sample.desc,
    ColStats(rawCounts.nozero, verbose = verbose, selected.stats = selected.stats))
  #stats.per.sample.nozero$Mcounts <- stats.per.sample.nozero$sum / 1e6
  outfiles["stats_per_sample_no-zero"] <- file.path(dir.tables.samples, "stats_per_sample_no-zero.tsv")
  message("\t\t", outfiles["stats_per_sample_no-zero"])
  write.table(x = stats.per.sample.nozero, file = outfiles["stats_per_sample_no-zero"], quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
  # names(stats.per.sample)
  # View(stats.per.sample.nozero)

  ## Compute statistics ommitting zero values
  message("\tComputing sample-wise statistics for filtered counts")
  stats.per.sample.filtered <- cbind(
    sample.desc,
    ColStats(filteredCounts, verbose = verbose, selected.stats = selected.stats))
  #stats.per.sample.filtered$Mcounts <- stats.per.sample.filtered$sum / 1e6
  # names(stats.per.sample)
  # View(stats.per.sample.nozero)
  outfiles["stats_per_sample_filtered_features"] <- file.path(dir.tables.samples, "stats_per_sample_filtered_features.tsv")
  message("\t\t", outfiles["stats_per_sample_filtered_features"])
  write.table(x = stats.per.sample.filtered, file = outfiles["stats_per_sample_filtered_features"], quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

  ## ---- Compute the counts per million reads with edgeR, using different lib sizes ----
  message("\tComputing normalized values with edgeR::cpm()")
  ## Note: the default normalization criterion (scaling by libbrary sum)
  ## is questionable because it is stronly sensitive to outliers
  ## (very highly expressed genes).  A more robust normalisation criterion
  ## is to use the 75th percentile, or the median. We use the median, somewhat arbitrarily,
  ## beause it gives a nice alignment on the boxplots.
  # scaledCounts.libsum <- cpm(filteredCounts.epsilon)    ## Counts per million reads, normalised by library sum
  # scaledCounts.perc75 <- cpm(filteredCounts.epsilon, lib.size = stats.per.sample.filtered$Q3)    ## Counts per million reads, normalised by 75th percentile
  # scaledCounts.perc95 <- cpm(filteredCounts.epsilon, lib.size = stats.per.sample.filtered$perc95)    ## Counts per million reads, normalised by 95th percentile
  # scaledCounts.median <- cpm(filteredCounts.epsilon, lib.size = stats.per.sample.filtered$median)    ## Counts per million reads, normalised by sample-wise median count


  ## ---- Compute standardized counts (percentile method) ----
  sample.norm.percentiles <- apply(filteredCounts.epsilon, 2, quantile, p = parameters$edgeR$norm.p)
  scaledCounts <- cpm(filteredCounts.epsilon, lib.size = sample.norm.percentiles)
  scaledCounts.log10 <- log10(scaledCounts) ## Log-10 transformed scaledCounts, xwith the epsilon for 0 counts
  scaledCounts.log2 <- log2(scaledCounts) ## Log-10 transformed scaledCounts, with the epsilon for 0 counts

  ## ???? WHY IS THIS SO DIFFERENT ????
  # sample.norm.factors <- calcNormFactors(filteredCounts.epsilon, method = "upperquartile") #, p = parameters$edgeR$norm.p)
  # plot(sample.norm.percentiles, sample.norm.factors)

  ## Export normalized counts (in log2-transformed counts per million reads)
  outfiles["scaledCounts"] <- file.path(dir.tables.samples, paste(sep = "", count.prefix, "_scaledCounts.tsv"))
  message("\tExporting standardized counts: ", outfiles["scaledCounts"])
  write.table(x = scaledCounts, row.names = TRUE, col.names = NA,
              file = outfiles["scaledCounts"], sep = "\t", quote = FALSE)

  outfiles["log2scaledCounts"] <- file.path(dir.tables.samples, paste(sep = "", count.prefix, "_scaledCounts_log2.tsv"))
  message("\tExporting log2-transformed standardized counts: ", outfiles["log2scaledCounts"])
  write.table(x = scaledCounts.log2, row.names = TRUE, col.names = NA,
              file = outfiles["log2scaledCounts"], sep = "\t", quote = FALSE)


  ## Compute Trimmed Means of M Values (TMM): TO BE DONE
  stats.per.sample.filtered$mean.std.count <- apply(scaledCounts, 2, mean, na.rm = TRUE)
  stats.per.sample.filtered$log2.mean.std.count <- apply(scaledCounts.log2, 2, mean, na.rm = TRUE)
  stats.per.sample.filtered$log10.mean.std.count <- apply(scaledCounts.log10, 2, mean, na.rm = TRUE)

  ## ---- Export stats per sample ----
  outfiles["stats_per_sample_filtered_features"] <- file.path(dir.tables.samples, "stats_per_sample_filtered_features.tsv")
  message("\t", "Exporting stats per sample for filtered features")
  message("\t\t", outfiles["stats_per_sample_filtered_features"])
  write.table(x = stats.per.sample.nozero, file = outfiles["stats_per_sample_filtered_features"], quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)


  ## ----print_sample_stats--------------------------------------------------
  ## Statistics per sample
  # setdiff(selected.stats, names(stats.per.sample))

  # names(stats.per.sample.all)
  report.text <- append(report.text, "\n\n## Sample-wise statistics\n")

  report.text <- append(report.text, "\n\n### All feature counts (zeros included)\n")
  report.text <- append(report.text, kable(stats.per.sample.all[,selected.stats], digits = 2,
                                         format.args = list(big.mark = ",", decimal.mark = "."),
                                         caption = "Sample-wise statistics for all features (zeros included)"))

  ## ----library_sizes_barplot, fig.width=6, fig.height=6, fig.cap="**Barplot of assigned reads per sample. ** Bars indicate the sum of read counts assigned to features (genes) per sample (library)."----
  par.ori <- par(no.readonly = TRUE) # Store original parameters
  figname <- "libsize_barplot_all_features"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)
  for (f in 1:length(figure.formats)) {
    fig.format <- figure.formats[f]
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 8, height = 8)

    LibsizeBarplot(counts = rawCounts, sample.labels = sample.labels, sample.colors = sample.desc$color, main = "All features", cex.axis = 0.8)

    silence <- dev.off(); rm(silence)
    if (f == 1) {
      report.text <- ReportFigure(name = figname, 
                                  figureFile = figure.file, 
                                  reportFile = rmd.report, 
                                  report.text = report.text, 
                                  out.width = parameters$DEG$out_width)
    }
    # system(paste("open", figure.file))
  }

  report.text <- append(report.text, "\n\n### All features, non-zero counts\n")
  report.text <- append(report.text, kable(stats.per.sample.nozero[,selected.stats], digits = 2,
                                         format.args = list(big.mark = ",", decimal.mark = "."),
                                         caption = "Sample-wise statistics for all features (zeros excluded)"))

  report.text <- append(report.text, "\n\n### Filtered features\n")
  report.text <- append(report.text, kable(stats.per.sample.filtered[,selected.stats], digits = 2,
                                         format.args = list(big.mark = ",", decimal.mark = "."),
                                         caption = "Sample-wise statistics for filtered features"))


  figname <- "libsize_barplot_filtered_features"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)
  for (f in 1:length(figure.formats)) {
    fig.format <- figure.formats[f]
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 7, height = 8)

    LibsizeBarplot(counts = filteredCounts, sample.labels = sample.labels, sample.colors = sample.desc$color, main = "After filtering", cex.names = 0.8)

    silence <- dev.off(); rm(silence)
    if (f == 1) {
      report.text <- ReportFigure(name = figname, 
                                  figureFile = figure.file, 
                                  reportFile = rmd.report, 
                                  report.text = report.text, 
                                  out.width = parameters$DEG$out_width)
      #      report.text <- ReportFigure(figname, figure.file, report.text, out.width = parameters$DEG$out_width)
    }
    # system(paste("open", figure.file))
  }

  ## ---- Barplot comparing library sizes before and after filtering ----
  figname <- "libsize_barplot_all_vs_filtered_features"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)
  for (f in 1:length(figure.formats)) {
    fig.format <- figure.formats[f]
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    message("\t\t\t", fig.format, " plot\t", figname)
    plot.height <- max(2 + 0.35 * nrow(sample.desc), 6)

    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 8, height = plot.height)

    par.ori <- par(no.readonly = TRUE)

    ## Prepare data frame for the comparative barplot
    x <- data.frame(
      "filtered" = stats.per.sample.filtered$Mcounts,
      "all" = stats.per.sample.all$Mcounts)
    row.names(x) <- sample.labels

    label.sizes <- nchar(row.names(x))
    left.margin <- 2 + max(label.sizes)*0.5
    # par("mar")
    par(mar = c(5.1, left.margin, 4.1, 1))

    x <- x[nrow(x):1,]
    indices <- sort(rep(1:nrow(x), times = 2))
    barplot.densities <- rep(c(50, -1), length.out = 2 * nrow(x))
    bplt <- barplot(
      as.matrix(t(x)),
      horiz = TRUE,
      beside = TRUE,
      las = 1,
      cex.names = 0.9,
      col = rev(as.vector(unlist(sample.desc$color[indices]))),
      density = barplot.densities, xlab = "Million reads",
      main = "Library sizes\nBefore vs after filtering",
      xlim = c(0, max(x) * 1.3))
    barplot.Mreads <- round(digits = 2, as.vector(unlist(t(x))))
    text(x = pmax(barplot.Mreads, 3), labels =  barplot.Mreads, y = bplt, pos = 2, font = 2, cex = 0.8)
    #sample.type <- rep(c("Filtered", "All"), length.out = 2 * nrow(x))
    #text(x = 0, labels = sample.type, y = bplt, pos = 4, font = 1, cex=0.5)
    legend("topright", legend = c("All genes", "Filtered"), col = "gray", density = c(-1, 50), cex = 0.7)

    par(par.ori)

    silence <- dev.off(); rm(silence)
    if (f == 1) {
      report.text <- ReportFigure(name = figname, 
                                  figureFile = figure.file, 
                                  reportFile = rmd.report, 
                                  report.text = report.text, 
                                  out.width = parameters$DEG$out_width)
      # report.text <- ReportFigure(
      #   figname, figure.file, report.text, out.width = parameters$DEG$out_width)
      #, chunk.opt = paste(sep = "", ", out.height = ", plot.height))
    }
    # system(paste("open", figure.file))
  }


  ## ---- Between-sample normalisation -------------------------------------------------------

  norm.comparison <- NormalizeCountTable(
    counts = filteredCounts, class.labels = sample.conditions, nozero = TRUE,
    method = norm.methods, percentile = norm.percentile, log2 = FALSE, epsilon = epsilon, detailed.sample.stats = TRUE,
    verbose = verbose)
  # names(norm.comparison)

  ## ----size_factors, fig.width=8, fig.height=8, fig.cap="Sample size factors for different normalisation methods. "----
  # m <- "median"
  size.factors <- data.frame(matrix(nrow = length(sample.ids), ncol = length(norm.methods)))
  colnames(size.factors) <- norm.comparison$method.name
  rownames(size.factors) <- sample.ids
  for (m in norm.comparison$method.name) {
    size.factors[,m] <- norm.comparison[[m]]$size.factor
  }

  figname <- "size_factors"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)
  for (f in 1:length(figure.formats)) {
    fig.format <- figure.formats[f]
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    # message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 7, height = 8)

    plot(size.factors, main = "Sample size factors", col = sample.desc$color)

    silence <- dev.off(); rm(silence)
    if (f == 1) {
      report.text <- ReportFigure(name = figname, 
                                  figureFile = figure.file, 
                                  reportFile = rmd.report, 
                                  report.text = report.text, 
                                  out.width = parameters$DEG$out_width)
      #      report.text <- ReportFigure(figname, figure.file, report.text, out.width = parameters$DEG$out_width)
    }
    # system(paste("open", figure.file))
  }

  ## ---- Read count distributions per sample ----
  report.text <- append(report.text, "\n\n## Distribution of read counts\n\n")
  
  ### Histograms: counts per gene (all samples together)
  report.text <- append(report.text, "\n\n### Histograms: counts per gene\n\n")
  
  figname <- "count_histograms"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)
  
  for (f in 1:length(figure.formats)) {
    par.ori <- par(no.readonly = TRUE)
    fig.format <- figure.formats[f]
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    # message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 12, height = 9)
    
    par(mfrow = c(2,2))
    h1 <- HistOfCounts(counts = rawCounts, maxPercentile = 95, 
                       las = 1, cex.axis = 0.8, legend.cex = 0.7,
                       main = "Raw counts", ylab = "Genes", 
                       col = "#FFEEDD", border = "#AA8877")
    h2 <- HistOfCounts(counts = filteredCounts, maxPercentile = 95, discardZeroRows = FALSE,
                       las = 1, cex.axis = 0.8, legend.cex = 0.7,
                       main = "Filtered counts", ylab = "Genes", 
                       col = "#DDEEFF", border = "#7788AA")
    ## log2(raw counts) histograme
    h3 <- HistOfCounts(counts = rawCounts.log2, 
                       maxPercentile = 100, discardZeroRows = FALSE, 
                       las = 1, cex.axis = 0.8, legend.cex = 0.7,
                       col = "#FFEEDD", border = "#AA8877",
                       xlab = "log2(counts)",
                       main = "Yeast Bdf1 vs WT\nRaw counts (log2)", ylab = "Genes")
    
    ## log2(filtered counts) histogram
    h4 <- HistOfCounts(counts = filteredCounts.log2, maxPercentile = 100,  
                       las = 1, cex.axis = 0.8, legend.cex = 0.7,
                       col = "#DDEEFF", border = "#7788AA",
                       xlab = "log2(counts)",
                       main = "Yeast Bdf1 vs WT\nfiltered counts (log2)", ylab = "Genes")
    
    par(mfrow = c(1,1))
    par(par.ori)
    silence <- dev.off(); rm(silence)
    if (f == 1) {
      report.text <- ReportFigure(name = figname, 
                                  figureFile = figure.file, 
                                  reportFile = rmd.report, 
                                  report.text = report.text, 
                                  out.width = parameters$DEG$out_width)
#      report.text <- ReportFigure(figname, figure.file, report.text, out.width = "95%")
    }
    # system(paste("open", figure.file))
  }
  
  
  ### Sample-wise read count box plots
  report.text <- append(report.text, "\n\n### Sample-wise read count box plots\n\n")
  
###  ```{r count_boxplots, fig.width=10, fig.height=10, fig.cap="Read count distributions. Top: raw counts. Bottom: counts per millon reads (scaledCounts). Left panels: linear scale, which emphasizes  outlier features denoted by very high counts. Rigt panels log counts permit to perceive the distribution of its whole range, including small count values. Null counts are replaced by an epsilon < 1, and appearas negative numbers after log transformation."}

  figname <- "boxplots_per_sample"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)

  ## Select a reasonable number of samples if there are too any of them
  boxplot.max.samples <- 30
  if (length(sample.ids) > boxplot.max.samples) {
    message("Selecting a subset of ", boxplot.max.samples, " samples",
            " among ", length(sample.ids), " for the boxplots. ")
    selected.samples <- sample(x = sample.ids, size = boxplot.max.samples)
  } else {
    selected.samples <- sample.ids
  }

  for (f in 1:length(figure.formats)) {
    par.ori <- par(no.readonly = TRUE)
    fig.format <- figure.formats[f]
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    # message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 10, height = 12)

    par(mfrow = c(2,2))

    ## Boxplot of raw counts
    BoxplotsPerSample(
      counts = rawCounts[, selected.samples],
      sample.desc = sample.desc[selected.samples,], 
      sample.label.col = "Label",
      xlab = "Raw counts",
      main = "Raw counts")

    ## Boxplot of log10-transformed counts
    BoxplotsPerSample(
      counts = rawCounts.log2[, selected.samples], 
      sample.desc = sample.desc[selected.samples,], 
      sample.label.col = "Label",
      xlab = "log2(counts)",
      main = "log2 raw counts")

    ## Boxplot of scaledCounts
    BoxplotsPerSample(
      counts = scaledCounts[, selected.samples], 
      sample.desc = sample.desc[selected.samples,], 
      sample.label.col = "Label",
      xlab = "std counts",
      main = "Normalized counts")
    
    ## Boxplot of log10-transformed scaledCounts
    BoxplotsPerSample(
      counts = scaledCounts.log2[, selected.samples], 
      sample.desc = sample.desc[selected.samples,], 
      sample.label.col = "Label",
      xlab = "log2(counts)",
      main = "log2 normalised counts")

    par(par.ori)

    silence <- dev.off(); rm(silence)
    if (f == 1) {
      report.text <- ReportFigure(name = figname, 
                                  figureFile = figure.file, 
                                  reportFile = rmd.report, 
                                  report.text = report.text, 
                                  out.width = parameters$DEG$out_width)
      #      report.text <- ReportFigure(figname, figure.file, report.text, out.width = parameters$DEG$out_width)
    }
    # system(paste("open", figure.file))
  }


  ## ---- Differential analysis  --------
  report.text <- append(report.text, "\n\n## Differential analysis\n\n")

  # i <- 1  # for quick test
  for (i in 1:nrow(design)) {


    prefix <- list() ## list for output file prefixes

    deg.results <- list()

    ## Identify samples for the first condition
    ref.condition <- as.vector(design[i, reference.column])  ## First condition for the current comparison
    ref.samples <- sample.ids[sample.conditions == ref.condition]
    if (length(ref.samples) < 2) {
      stop(paste("Cannot perform differential analysis. The count table contains less than 2 samples for condition", ref.condition))
    }

    ## Identify samples for the second condition
    test.condition <- as.vector(design[i,test.column])  ## Second condition for the current comparison
    test.samples <- sample.ids[sample.conditions == test.condition]
    if (length(test.samples) < 2) {
      stop(paste("Cannot perform differential analysis. The count table contains less than 2 samples for condition", test.condition))
    }


    report.text <- append(report.text, paste(sep= "", "\n\n### ", test.condition, " versus ", ref.condition, "\n\n"))
#    report.text <- append(report.text, paste(sep= "", "Test condition : ",  test.condition, "\n"))
#    report.text <- append(report.text, paste(sep= "", "Reference condition :",  ref.condition, "\n"))
    compa.table <- data.frame(
      type = c("Reference", "Test"),
      condition = c(ref.condition, test.condition),
      samples = c(length(ref.samples), length(test.samples))
    )
    report.text <- append(report.text,
                         kable(compa.table))


    message("\tDifferential analysis\t", i , "/", nrow(design), "\t", ref.condition, " vs ", test.condition)

    ## Create a specific directory for the results of this comparison
    comparison.prefix <- comparison.summary$prefixes[i]
    dir.results.diffexpr <- file.path(result.dir, paste(sep = "", comparison.prefix))
    dirs[comparison.prefix] <- dir.results.diffexpr ## Index current result dir for the list of directories
    comparison.summary[i, "result.dir"] <- dir.results.diffexpr ## Include current result dir to the comparison summary table
    dir.create(path = file.path(dirs["main"], dir.results.diffexpr), showWarnings = FALSE, recursive = TRUE)
    prefix["comparison_file"] <- file.path(dir.results.diffexpr, comparison.prefix)
    message("\t\tDifferential expression results:\t", dir.results.diffexpr)

    ## Create a specific directory for the figures of this comparison
    dir.figures.diffexpr <-  file.path(dirs[comparison.prefix], "figures")
    dirs[paste(sep = "_", comparison.prefix, "figures")] <- dir.figures.diffexpr ## Index current figures dir for the list of directories
    comparison.summary[i, "figures"] <- dir.figures.diffexpr ## Include current figures dir to the comparison summary table
    dir.create(path = file.path(dirs["main"], dir.figures.diffexpr), showWarnings = FALSE, recursive = TRUE)
    prefix["comparison_figure"] <- file.path(dir.figures.diffexpr, comparison.prefix)
    message("\t\tfigures:\t", dir.figures.diffexpr)
    #    paste(sep = "", comparison.prefix, "_",  suffix.deg))


    ## Select counts for the samples belonging to the two conditions
    current.samples <- c(ref.samples, test.samples)
    current.counts <- data.frame(filteredCounts[,current.samples])
    # dim(current.counts)  ## For test
    # names(current.counts)

    if (sum(!names(current.counts) %in% sample.ids) > 0) {
      stop("Count table contains column names without ID in sample description file.")
    }

    ## Define conditions and labels for the samples of the current analysis
    current.conditions <- sample.conditions[current.samples]
    current.labels <- paste(current.conditions, names(current.counts), sep = "_")

    ## A big result table with all features and all statistics
    result.table <- init.deg.table(scaledCounts, ref.samples, test.samples, stats = FALSE)
    # View(result.table)
    # dim(result.table)

    ## A synthetic table with significant featrues only
    result.table.synthetic <- data.frame(row.names = row.names(result.table),"feature_id" = row.names(result.table))

    ## ---- DESeq2 analysis ----
    message("\tDESeq2 analysis\t", comparison.prefix)
    deseq2.result <- deseq2.analysis(
      counts = current.counts,
      # head(counts)
      condition = current.conditions,
      ref.condition = test.condition,
      comparison.prefix = comparison.prefix,
      thresholds = as.vector(thresholds),
      title = comparison.prefix,
#      verbose = verbose,
      dir.figures = dir.figures.diffexpr)
    deg.results[["DESeq2"]] <- deseq2.result
    # names(deg.results[["DESeq2"]])
    #  attributes(deg.results[["DESeq2"]]$dds)
    # View(deg.results[["DESeq2"]]$result.table)

    #  head(rownames(deseq2.result$result.table))
    # head(rownames(result.table))
    # x <- rownames(result.table)
    # y <- rownames(deseq2.result$result.table[row.names(result.table),])
    # sum(x != y)
    # names(result.table)
    result.table <- cbind(
      result.table,
      "DESeq2" = deseq2.result$result.table[row.names(result.table),])

    result.table.synthetic <- cbind(
      result.table.synthetic,
      "DESeq2" = deseq2.result$result.table[row.names(result.table),c("padj", "FC")])
    # names(result.table)
    # names(result.table.synthetic)
    # dim(deseq2.result$result.table)
    # dim(deseq2.result$result.table)
    # names(deseq2.result$result.table)
    # View(deseq2.result$result.table)
    # View(result.table)


    ## Save the completed DESeq2 result table
    deseq2.result.file <- paste(sep = "", prefix["comparison_file"], "_DESeq2.tsv")
    comparison.summary[i,"deseq2"] <- deseq2.result.file
    message("\tExporting DESeq2 result table (tab): ", deseq2.result.file)
    outfiles[paste(sep = "_", comparison.prefix, "DESeq2")]  <- deseq2.result.file
    write.table(
      x = deseq2.result$result.table, row.name    = FALSE,
      file = deseq2.result.file,
      sep = "\t", quote = FALSE)

    ## ---- edgeR analysis ----
    # norm.method <- "TMM" ## For quick test and debugging
    for (norm.method in parameters$edgeR$norm_method) {

      edgeR.prefix <- paste(sep = "_", "edgeR", norm.method)

      edger.result <- edger.analysis(
        counts = current.counts,
        condition = current.conditions,
        test.condition = test.condition,
        ref.condition = ref.condition,
        thresholds = as.vector(thresholds),
        comparison.prefix = comparison.prefix,
        norm.method = norm.method,
        title = paste(sep = "_", norm.method, comparison.prefix),
#        verbose = verbose,
        dir.figures = dir.figures.diffexpr)
      deg.results[[edgeR.prefix]] <- edger.result

      ## A tricky way to add edgeR with normalisation in column names
      edger.to.bind <- edger.result$result.table[row.names(result.table),]
      colnames(edger.to.bind) <- paste(sep = "_", edgeR.prefix, colnames(edger.to.bind))
      # names(edger.to.bind)
      # View(edger.to.bind)
      result.table <- cbind(
        result.table,
        edger.to.bind)
      # names(result.table)

      colnames.synthetic <- paste(sep = "_", edgeR.prefix, c("padj", "FC"))
      result.table.synthetic <- cbind(
        result.table.synthetic,
        edger.to.bind[, colnames.synthetic])
      # names(edger.to.bind)

      ## Export edgeR result table
      edger.result.file <- paste(sep = "", prefix["comparison_file"], "_", edgeR.prefix, ".tsv")
      comparison.summary[i,"edger"] <- edger.result.file
      message("\tExporting edgeR result table (tab): ", edger.result.file)
      outfiles[paste(sep="_", comparison.prefix, edgeR.prefix)] <- edger.result.file
      write.table(x = edger.result$result.table,
                  file = edger.result.file,
                  row.names = FALSE,
                  sep = "\t", quote = FALSE)
    }

    ## ----  Select significant genes by combining DESeq2 and edgeR results ----
    padj.columns <- grep(colnames(result.table.synthetic), pattern = "padj", value = TRUE)
    FC.columns <- grep(colnames(result.table.synthetic), pattern = "FC", value = TRUE)
    if (is.null(parameters$DEG$selection_criterion)) {
      parameters$DEG$selection_criterion <- "union"
    }

    pval.passed <- result.table.synthetic[, padj.columns] <= thresholds$padj
    # table(is.na(pval.passed))
    pval.passed[is.na(pval.passed)] <- FALSE
    #  table(is.na(pval.passed))
    # table(unlist(pval.passed))
    FC.passed <- result.table.synthetic[, FC.columns] >= thresholds$FC
    positive <- pval.passed & FC.passed
    colnames(positive) <- names(deg.results)
    # View(positive)
    # table(positive[,1], positive[,2])

    positive.col.name <- paste(sep = "", "positive_", parameters$DEG$selection_criterion)

    if (parameters$DEG$selection_criterion == "union") {
      combined.positive <- apply(positive, 1, sum) > 0
    } else if (parameters$DEG$selection_criterion == "intersection") {
      combined.positive <- apply(!positive, 1, sum) == 0
    } else if (parameters$DEG$selection_criterion == "DESeq2") {
      combined.positive <- positive[, "DESeq2"]
    } else if (parameters$DEG$selection_criterion == "edgeR") {
      ## In case several normalisation methods would have been selected for edgeR, take the intersection between them (but ignore DESeq2)
      combined.positive <- apply(as.matrix(!positive[, setdiff(colnames(positive), "DESeq2")]), 1, sum) == 0
    }
    # nrow(result.table.synthetic)
    # table(result.table.synthetic[, positive.col.name])
    result.table[, positive.col.name] <- combined.positive
    result.table.synthetic[, positive.col.name] <- combined.positive
    result.table.synthetic <- result.table.synthetic[result.table.synthetic[, positive.col.name],]
    # nrow(result.table.synthetic)

    ## ---- Export DEG result tables ----

    ## Export full result table (DESeq2 + edgeR with different normalisation methods)
    ## in a tab-separated values (tsv) file
    result.file <- paste(sep = "",
                         prefix["comparison_file"],
                         "_diffexpr_DESeq2_and_edgeR.tsv")
    # comparison.summary[i,"result.table"] <- paste(sep=".", result.file, "tsv")
    message("\tExporting result table (tsv): ", result.file)
    outfiles[paste(sep = "", comparison.prefix, "_complete_result_table")] <- result.file
    write.table(x = result.table, row.names = FALSE,
                file = result.file, sep = "\t", quote = FALSE)

    ## Export synthetic table with positive features, and only the most relevant stats
    deg.file <- paste(sep = "",
                         prefix["comparison_file"],
                         "_diffexpr_DESeq2_and_edgeR_DEG.tsv")
    # comparison.summary[i,"result.table"] <- paste(sep=".", result.file, "tsv")
    message("\tExporting table of differentially expressed genes (tsv): ", deg.file)
    outfiles[paste(sep = "", comparison.prefix, "_synthetic_result_table")] <- deg.file
    write.table(x = result.table.synthetic, row.names = FALSE,
                file = deg.file, sep = "\t", quote = FALSE)



    ## Collect results by output statistics
    deg.compa <- list()
    feature.ids <- row.names(current.counts)
    stats.to.collect <- c("padj", "FC", "log2FC")
    for (stat in stats.to.collect) {
      message("Collecting ", stat, " from alternative DEG results. ")
      deg.compa[[stat]] <- data.frame(
        matrix(nrow = nrow(current.counts),
               ncol = length(names(deg.results))))
      colnames(deg.compa[[stat]]) <- names(deg.results)
      rownames(deg.compa[[stat]]) <- feature.ids
      # deg.name <- "DESeq2"
      for (deg.name in names(deg.results)) {
        deg.compa[[stat]][feature.ids, deg.name] <-
          as.vector(deg.results[[deg.name]]$result.table[feature.ids,stat])
      }
      #    View(deg.compa[[stat]])
    }




    ## Define feature colors according to their level of expression (count means)
    # feature.scores <- log2(apply(filteredCounts.epsilon, 1, median))

    # hist(feature.scores, breaks = 100)

    # View(deg.compa$padj)
    ## compare DESeq2 and edgeR normalisatio results

    ## ---- Plot comparing DEGs obtained with DESeq2 and edgeR ----

    ## Comparison between adjusted p-values
    report.text <- append(report.text, "\n\n#### Adjusted p-values\n\n")

    figname <- paste(sep = "", comparison.prefix, "_norm_compa_padj")
    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (f in 1:length(figure.formats)) {
      fig.format <- figure.formats[f]
      figure.file <- paste(sep = "", file.prefix, ".", fig.format)
      figure.files[[fig.format]][figname] <- figure.file
      # message("\t\t\t", fig.format, " plot\t", figname)
      OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 7, height = 8)

      plot(deg.compa$padj, log = "xy",
           #       col = FeatureColors(palette.type = "2col", scores = feature.scores),
           col = FeatureColors(palette.type = "dens",
                               x = deg.compa$padj[,1], y = deg.compa$padj[,2]),
           main = paste(sep = "", comparison.prefix, "\nAdjusted p-values"))

      silence <- dev.off(); rm(silence)
      if (f == 1) {
        report.text <- ReportFigure(name = figname, 
                                    figureFile = figure.file, 
                                    reportFile = rmd.report, 
                                    report.text = report.text, 
                                    out.width = parameters$DEG$out_width)
        #        report.text <- ReportFigure(figname, figure.file, report.text, out.width = parameters$DEG$out_width)
      }
      # system(paste("open", figure.file))
    }

    report.text <- append(report.text, "\n\n#### Log2(fold changes)\n\n")
    figname <- paste(sep = "", comparison.prefix, "_norm_compa_log2FC")
    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (f in 1:length(figure.formats)) {
      fig.format <- figure.formats[f]
      figure.file <- paste(sep = "", file.prefix, ".", fig.format)
      figure.files[[fig.format]][figname] <- figure.file
      # message("\t\t\t", fig.format, " plot\t", figname)
      OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 7, height = 8)

      plot(deg.compa$log2FC,
           #       col = FeatureColors(palette.type = "2col", scores = feature.scores),
           col = FeatureColors(palette.type = "dens",
                               x = deg.compa$log2FC[,1], y = deg.compa$log2FC[,2]),
           main = paste(sep = "", comparison.prefix, "\nlog2(fold change)"))

      silence <- dev.off(); rm(silence)
      if (f == 1) {
        report.text <- ReportFigure(name = figname, 
                                    figureFile = figure.file, 
                                    reportFile = rmd.report, 
                                    report.text = report.text, 
                                    out.width = parameters$DEG$out_width)
        #        report.text <- ReportFigure(figname, figure.file, report.text, out.width = parameters$DEG$out_width)
      }
      # system(paste("open", figure.file))
    }


    ## ---- Draw Volcano plots -----
    report.text <- append(report.text, "\n\n#### Volcano plots\n\n")
    deg.names <- names(deg.results)
    nb.panels <- n2mfrow(length(deg.names))
    if (length(deg.names) <= 5) {
      nb.panels <- rev(nb.panels)
    }
    figname <- paste(sep = "", comparison.prefix, "_norm_compa_volcano_plots")

    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (f in 1:length(figure.formats)) {
      fig.format <- figure.formats[f]
      figure.file <- paste(sep = "", file.prefix, ".", fig.format)
      figure.files[[fig.format]][figname] <- figure.file
      # message("\t\t\t", fig.format, " plot\t", figname)
      height <- 2 + nb.panels[1]*3
      width <- 2 + nb.panels[2]*3.5
      OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, height = height, width = width)

      par.ori <- par(no.readonly = TRUE)
      par(mfrow = nb.panels)
      for (deg.name in deg.names) {
        message("\tDrawing volcano plot\t", deg.name)
        deg.table <- deg.results[[deg.name]]$result.table
        VolcanoPlot.MultiTestTable(
          multitest.table = deg.results[[deg.name]]$result.table,
          main = deg.name,
          effect.size.col = "log2FC",
          control.type = "padj",
          alpha = parameters$DEG$thresholds$padj,
          effect.threshold = log2(parameters$DEG$thresholds$FC),
          legend.cex = 0.9,
          legend.corner = "top")
      }
      par(mfrow = c(1,1))
      par(par.ori)
      silence <- dev.off(); rm(silence)
      if (f == 1) {
        report.text <- ReportFigure(name = figname, 
                                    figureFile = figure.file, 
                                    reportFile = rmd.report, 
                                    report.text = report.text, 
                                    out.width = parameters$DEG$out_width)
        #        report.text <- ReportFigure(figname, figure.file, report.text, out.width = parameters$DEG$out_width)
      }
      # system(paste("open", figure.file))
    }


    ## system(paste("ls -ltr ", pdf.files[prefix]))



    ## ---- Comparison between p-value histograms -----
    ## Draw Volcano plots
    # deg.name <- "DESeq2"
    # deg.name <- "edgeR_TMM"
    report.text <- append(report.text, "\n\n#### P-value histograms\n\n")
    figname <- paste(sep = "", comparison.prefix, "_norm_compa_pvalue_histograms")
    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (f in 1:length(figure.formats)) {
      fig.format <- figure.formats[f]
      figure.file <- paste(sep = "", file.prefix, ".", fig.format)
      figure.files[[fig.format]][figname] <- figure.file
      # message("\t\t\t", fig.format, " plot\t", figname)
      OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, height = 2 + nb.panels[1]*3, width = 2 + nb.panels[2]*3.5)


      par.ori <- par(no.readonly = TRUE)
      par(mfrow = nb.panels)
      # deg.name <- "DESeq2"
      for (deg.name in deg.names) {
        message("\tDrawing P value histogram\t", deg.name)
        deg.table <- deg.results[[deg.name]]$result.table
        degMultiTest <- multipleTestingCorrections(p.values = deg.results[[deg.name]]$result.table$pvalue)
        PlotPvalDistrib.MultiTestTable(
          multitest.result = degMultiTest,
          main = paste(sep = "", deg.name, "\nP-value histogram"),
          ylab = "Number of features",
          legend.cex = 0.8,
          draw.mean.line = TRUE, legend.corner = "topright"
        )
      }
      par(mfrow = c(1,1))
      par(par.ori)

      silence <- dev.off(); rm(silence)
      if (f == 1) {
        report.text <- ReportFigure(name = figname, 
                                    figureFile = figure.file, 
                                    reportFile = rmd.report, 
                                    report.text = report.text, 
                                    out.width = parameters$DEG$out_width)
#        report.text <- ReportFigure(figname, figure.file, report.text, out.width = parameters$DEG$out_width)
      }
      # system(paste("open", figure.file))
    }


    ## ---- Draw Venn diagram with number of genes declared significant ----
    ## according to the selection criteria (threshold fields).
    report.text <- append(report.text, "\n\n#### Venn diagrams: padj versus FC\n\n")
    selection.fields <- c("padj", "FC")
    selection.thresholds <- thresholds[selection.fields]
    selection.columns <- paste(sep = "", selection.fields, "_", selection.thresholds)

    figname <- paste(sep = "", comparison.prefix, "_norm_compa_Venn_", paste(collapse = "_", selection.fields))

    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (f in 1:length(figure.formats)) {
      fig.format <- figure.formats[f]
      figure.file <- paste(sep = "", file.prefix, ".", fig.format)
      figure.files[[fig.format]][figname] <- figure.file
      # message("\t\t\t", fig.format, " plot\t", figname)
      height <- 2 + nb.panels[1]*3
      width <- 2 + nb.panels[2]*3.5
      OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, height = height, width = width)

      par.ori <- par(no.readonly = TRUE)
      par(mfrow = nb.panels)
      # deg.name <- "DESeq2"
      for (deg.name in deg.names) {
        message("\tDrawing Venn diagram\t", deg.name)
        message("\t\tSelected columns\t", selection.columns)
        deg.table <- deg.results[[deg.name]]$result.table
        selection.venn.counts <- vennCounts(deg.table[,selection.columns])
        limma::vennDiagram(selection.venn.counts, cex = 1,
                           main = paste(comparison.prefix, "\n", deg.name, "selected features"),
                           circle.col = c("orange", "blue"), mar = c(0,0,5,0))
      }
      silence <- dev.off(); rm(silence)
      if (f == 1) {
        report.text <- ReportFigure(name = figname, 
                                    figureFile = figure.file, 
                                    reportFile = rmd.report, 
                                    report.text = report.text, 
                                    out.width = parameters$DEG$out_width)
        #        report.text <- ReportFigure(figname, figure.file, report.text, out.width = parameters$DEG$out_width)
      }
      # system(paste("open", figure.file))
    }


  }


  ## ---- Index input / output files and directories -----
  index <- data.frame()
  index <- rbind(index,
                 data.frame(
                   type = "directory",
                   name = names(dirs),
                   path = dirs
                 ))
  index <- rbind(index,
                 data.frame(
                   type = "input file",
                   name = names(infiles),
                   path = infiles
                 ))
  index <- rbind(index,
                 data.frame(
                   type = "output file",
                   name = names(outfiles),
                   path = outfiles
                 ))

  index$html.link <- paste("<a href='", index$path, "'>", index$name, "</a>")
  index.file <- "file.index.tsv"
  write.table(x = index, file = index.file, quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = TRUE)



  ## ---- Build the file index ----

  report.text <- append(report.text, "\n\n## Directories and files\n")

  report.text <- append(report.text, paste(sep = "", "- Config file : ", ReportLink(source = rmd.report, target = configFile)))
  report.text <- append(report.text, paste(sep = "", "- Count table : ", ReportLink(source = rmd.report, target = countFile)))
  report.text <- append(report.text, paste(sep = "", "- Sample descriptions : ", ReportLink(source = rmd.report, target = infiles["sample descriptions"])))
  report.text <- append(report.text, paste(sep = "", "- Design : ", ReportLink(source = rmd.report, target = infiles["design"])))

  report.text <- append(report.text, paste(sep = "", "\n\n### Directories\n"))
  report.text <- append(report.text, paste(sep = "", "| Content | Path |"))
  report.text <- append(report.text, paste(sep = "", "|:----------------------|--------------------------------------------------|"))
  for (dirname in names(dirs)) {
    dir <- dirs[dirname]
    report.text <- append(report.text, paste(sep = "", "| ", dirname, " | ", ReportLink(source = rmd.report, target = dir), " | "))
  }

  report.text <- append(report.text, paste(sep = "", "\n\n### Input files\n"))
  report.text <- append(report.text, paste(sep = "", "| Content | Path |"))
  report.text <- append(report.text, paste(sep = "", "|:----------------------|--------------------------------------------------|"))
  for (filename in names(infiles)) {
    file <- infiles[filename]
    report.text <- append(report.text, paste(sep = "", "| ", filename, " | ", ReportLink(source = rmd.report, target = file), " | "))
  }

  report.text <- append(report.text, paste(sep = "", "\n\n### Output files"))
  report.text <- append(report.text, paste(sep = "", "| Content | Path |"))
  report.text <- append(report.text, paste(sep = "", "|:----------------------|--------------------------------------------------|"))
  for (filename in names(outfiles)) {
    file <- outfiles[filename]
    report.text <- append(report.text, paste(sep = "", "| ", filename, " | ", ReportLink(source = rmd.report, target = file), " | "))
  }

  f <- 1
  for (f in 1:length(figure.formats)) {
    fig.format <- figure.formats[f]
    figfiles <- unlist(figure.files[fig.format])
    report.text <- append(report.text, paste(sep = "", "\n\n### Figures (", fig.format,")"))
    report.text <- append(report.text, paste(sep = "", "| Content | Path |"))
    report.text <- append(report.text, paste(sep = "", "|----------------------|--------------------------------------------------|"))
    for (filename in names(figfiles)) {
      file <- figfiles[filename]
      report.text <- append(report.text, paste(sep = "", "| ", filename, " | ", ReportLink(source = rmd.report, target = file), " | "))
    }
  }

  ## ---- Session info ---------------------------------------------------------
  ## Print the complete list of libraries + versions used in this session
  report.text <- append(report.text, "\n\n## Session info\n")
  session.info <- sessionInfo()
  # print(session.info)
  report.text <- append(report.text, "\n```\n")
  report.text <- append(report.text, capture.output(session.info))
  report.text <- append(report.text, "\n```\n")



  ## ---- Generate the HTML index -----
  report.text <- append(report.text, paste(sep = "", "\n\nJob done: ", Sys.time()))

  writeLines(text = report.text, con = report.socket)
  close(report.socket)
  rmarkdown::render(rmd.report, output_format = "html_document")
  rmarkdown::render(rmd.report, output_format = "pdf_document")


  ## ---- job_done ------------------------------------------------------------
  message("Tables directory\t", dir.tables.samples)
  message("Figures directory\t", dir.figures.samples)
  message("Index of input/output files\t", index.file)
  message("Report in Rmd format\t", rmd.report)
  message("Report in HTML format\t", sub(pattern = '.Rmd', replacement  = ".html", x = rmd.report))
  message("Report in pdf format\t", sub(pattern = '.Rmd', replacement  = ".pdf", x = rmd.report))
  message("Job done")


  detach(dataset)
  return()

}
