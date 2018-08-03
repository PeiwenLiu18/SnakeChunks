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
#' @param count.table file containing a table with counts of reads per feature (row) in each sample (column).
#' @param configFile a yaml-formatted file defining the mandatory + some optional parameteres
#' @param main.dir=getwd()  directory from which the script runs. Paths are defined relative to this directory.
#' @param result.dir=file.path(dir.main,"results") directory where the results will be stored. Should be defined relative to main directory.
#' @param count.prefix=NULL basename to export the transformed count tables. If not specified, computed automatically from the name of the count table file.
#' @param edgeR.norm.methods="TMM" vector with one or more normalisation methods for edegR. Supported: TMM, RLE, upperquartile, quantiles, none
#' @param verbose=1 level of verbosity
#'
#' @export
RNAseqAnalysis <- function(count.table,
                           configFile,
                           main.dir = getwd(),
                           result.dir = file.path(dir.main, "results"),
                           edgeR.norm.methods = "TMM",
                           count.prefix = NULL,
                           verbose = 1) {


  ## ---- TO DO ----
  ##
  ## Check non-used variables: suppress or use them
  ## - cols.heatmap
  ## - DEG.selection.criterion
  ## - filtered.counts.log10
  ## - filtered.counts.log2
  ## - stdcounts.libsum
  ## - stdcounts.perc95
  ## - stdcounts.median
  ## - current.labels
  
  ## ---- Define the main parameters -----------------------------------------------------

  ## Prepare a list of input and output files.
  ## This will later serve to generate a report.
  message("\tInitializing indexes and main parameters")
  dirs <- vector()
  infiles <- vector()   ## Input files
  outfiles <- vector()  ## For tab-separated value files


  ## Figures in different format
  figure.formats <- c("png", "pdf")
  figure.files <- list()
  for (fig.format in figure.formats) {
    figure.files[[fig.format]] <- vector()
  }


  ## Define main parameters to generate this report
  #dirs["main"] <- "~/ko-rna-seq/" ## Main directory
  dirs["main"] <- main.dir
  setwd(dirs["main"])
  message("\tMain directory: ", dirs["main"])

  ## Define YAML configuration file
  #configFile <- "metadata/config_RNA-seq.yml"
  infiles["config"] <- configFile
  message("\tConfiguration file: ", configFile)

  ## Prefix for the count table
  #count.prefix <- "bowtie2_featureCounts_all"
  #message("\tPrefix for the count table: ", count.prefix)
  infiles["count_table"] <- count.table
  message("\tCount table: ", count.table)

  ## ---- Initialize the Rmd report (index of input/output file)
  index.Rmd <- "index.Rmd"
  index.socket <- file(index.Rmd)

  Rmd.header <- '---
title: "RNA-seq analysis report"
author:
name: "[AUTHOR NAME]"
email: "[AUTHOR EMAIL]"
date: Last update:`r format(Sys.time())`
output:
  html_document:
    fig_caption: yes
    highlight: zenburn
    self_contained: yes
    theme: cerulean
    toc: yes
    toc_depth: 3
    toc_float: yes
  pdf_document:
    fig_caption: yes
    highlight: zenburn
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
  fig.width = 7,
  fig.height = 5,
  fig.path = "figures/")
```
'

  index.text <- Rmd.header


  ## ---- knitr_setup, include=FALSE,  eval=TRUE, echo=FALSE, warning=FALSE----

  quick.test <- FALSE ## For debug

  ## To do: check if these options are taken into account at rendering
  knitr::opts_chunk$set(
    fig.path = "figures/",
    echo = FALSE,
    eval = TRUE,
    cache = FALSE,
    message = FALSE,
    warning = FALSE)

  ## Load required libraries
  required.cran.libraries <- c("knitr",
                               "yaml",
                               "pander",
                               # "xlsx",
                               "ascii",
                               "xtable",
                               "gplots",
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

  ## ---- Load configuration file (YAML-formatted) ----
  if (!exists("configFile")) {
    ## The prompt does not seem to work with the Rmd documents
    #   message("Choose the parameter file")
    #   parameter.file <- file.choose()
    stop("This report requires to specify a variable named configFile, containing the path to an YAML-formatted file describing the parameters for this analysis.")
  }
  parameters <- yaml.load_file(configFile)
  message("\tLoaded parameters from file ", configFile)

  if (is.null(parameters$edgeR$norm_method)) {
    edgeR.norm.methods <- c("TMM","RLE","upperquartile","none")
    #edgeR.norm.methods <- c("TMM","RLE","upperquartile","none")
  } else {
    edgeR.norm.methods <- parameters$edgeR$norm_method
  }


  ## Compute count prefix, i.e. the basename to export various transformations of the count tables
  if (is.null(parameters$dir$count_prefix)) {
    ## use the basename of the count table as prefix for output file
    count.prefix <- basename(path = count.table)

    ## Suppress usual extensions from the basename
    for (ext in c(".tsv", ".tab", ".txt", ".csv")) {
      count.prefix <- sub(pattern = ext, replacement = "", x = count.prefix)
    }
  } else {
    count.prefix <- parameters$dir$count_prefix
  }
  message("\tFile prefix for normalized count tables ", count.prefix)

  ## Directory to store differential expression results
  dirs["output"] <- result.dir ## Index dir for the report
  message("\tDirectory for differential expression: ", result.dir)
  dir.create(result.dir, showWarnings = FALSE, recursive = TRUE)

  ## Directory for Figures
  dir.figures.samples <- file.path(result.dir, "figures")
  dirs["sample_figures"] <- dir.figures.samples
  message("\tDirectory for the sample-related figures: ", dir.figures.samples)
  dir.create(dir.figures.samples, showWarnings = FALSE, recursive = TRUE)

  ## Directory to export result files in tba-separated value (tsv) format
  dir.tables.samples <- file.path(result.dir, "samples/tables")
  dirs["tsv"] <- dir.tables.samples ## Index directory for the report
  message("\tDirectory to export sample-related tables:\t", dir.tables.samples)
  dir.create(dir.tables.samples, showWarnings = FALSE, recursive = TRUE)

  ## ----default_parameters--------------------------------------------------
  ## In this chunk, we define a set of default parameters for the display and the analysis. These parameters can be modified but it is not necessary to adapt them to each project.
  if ((!exists("verbosity")) || (is.null(verbosity))) {
    verbosity <- 1
  }

  ## Color palette for heatmaps. I like this Red-Blue palette because
  ## - it suggests a subjective feeling of warm (high correlation)/cold (low correlation)
  ## - it can be seen by people suffering from red–green color blindness.
  if (!exists("cols.heatmap")) {
    cols.heatmap <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))
  }

  ## A trick: to enable log-scaled plots for 0 values, I add an epsilon increment
  if (is.null(parameters$DEG$epsilon)) {
    parameters$DEG$epsilon <- 0.1 # passed to file parameters.R, 2017-03-15
  }
  epsilon <- parameters$DEG$epsilon

  ## Default method for the selection of the final list of DEG
  ## Since DESeq2 seems more conservative, we use it by default
  if (is.null(parameters$DEG$selection_criterion)) {
    parameters$DEG$selection_criterion <- "DESeq2"
  }
  DEG.selection.criterion <- parameters$DEG$selection_criterion

  ## Sample description file
  if (is.null(parameters$metadata$samples)) {
    stop("The sample file must be defined in the metadata seection of the yaml config file: ", configFile)
  } else {
    infiles["sample descriptions"] <- parameters$metadata$samples
  }

  ## Design file
  if (is.null(parameters$metadata$design)) {
    stop("The design file must be defined in the metadata seection of the yaml config file: ", configFile)
  } else {
    infiles["design"] <- parameters$metadata$design
  }

  ## Count table
  # infiles["counts"] <- file.path(parameters$dir$diffexpr, paste(sep = "", count.prefix, ".tsv")
  infiles["counts"] <- count.table
  # infiles["counts"] <- file.path(dirs["main"], infiles["counts"])
  if (!file.exists(count.table)) {
    stop("Feature count table does not exist: ", count.table)
  } else {
    message("\tFeature count table: ", count.table)
  }



  ## ---- threshold_table -----------------------------------------------------
  if (is.null(parameters$DEG$thresholds)) {
    message("\tDEG thresholds were not defined in config file -> using default values")
    if (is.null(parameters$DEG)) {
      parameters$DEG <- list()
    }
    parameters$DEG$thresholds <- list(
      min.count = 1,
      mean.count = 5,
      padj = 0.05,
      FC = 1.2)
  }
  thresholds <- as.vector(parameters$DEG$thresholds)

  ## Print the threshold tables
  ## NOTE: I should evaluate what I do with the kable calls
  index.text <- append(index.text, "\n\n## Parameters\n")
  index.text <- append(
    index.text,
    kable(t(as.data.frame(thresholds)), col.names = "Threshold",
          caption = "Thresholds for the selection of differentially expressed genes. "))

  outfiles["thresholds"] <- file.path(dir.tables.samples, "thresholds.tsv")
  write.table(x = t(as.data.frame(thresholds)),
              file = outfiles["thresholds"],
              sep = "\t", row.names = TRUE, col.names = FALSE)
  # list.files(dir.tables.samples)
  # system(paste("open", dir.tables.samples))

  ## ----read_samples--------------------------------------------------------
  #setwd(dirs["main"]) ## !!!!! I don't understand why I have to reset the working directory at each chunk

  ## Read the sample description file, which indicates the
  ## condition associated to each sample ID.
  message("\tReading sample description file: ", infiles["sample descriptions"])
  if (!file.exists(infiles["sample descriptions"])) {
    stop("Sample description file does not exist\t", infiles["sample descriptions"])
  }
  sample.desc <- read.delim(
    file = infiles["sample descriptions"],
    sep = "\t",
    comment = ";", header = TRUE, row.names = 1)
  sample.ids <- row.names(sample.desc)
  message("\t\tNb of samples = ", length(sample.ids))
  # head(sample.desc)


  ## Experimental conditions
  sample.conditions <- as.vector(sample.desc[,1]) ## Condition associated to each sample
  names(sample.conditions) <- sample.ids
  # print(sample.conditions)

  ## Build sample labels by concatenating their ID and condition
  if (is.null(sample.desc$Label)) {
    sample.labels <- paste(sep = "_", sample.ids, sample.conditions)
    if (!is.null(sample.desc$Replicate)) {
      sample.labels <- paste(sep = "_", sample.labels, sample.desc$Replicate)
    }
    sample.desc$Label <- sample.labels
  } else {
    sample.labels <- as.vector(unlist(sample.desc$Label))
  }
  # print(sample.labels)


  ## Define a specific color for each distinct condition
  conditions <- unique(sample.conditions) ## Set of distinct conditions
  cols.conditions <- brewer.pal(max(3, length(conditions)),"Dark2")[1:length(conditions)]
  names(cols.conditions) <- conditions
  # print(cols.conditions)
  message("\t\tNb of conditions = ", length(conditions))
  message("\t\tConditions = ", paste(collapse = ", ", conditions))



  ## Define a color per sample according to its condition
  sample.desc$color <- cols.conditions[sample.conditions]
  # names(cols.samples) <- sample.ids
  # print(cols.samples)

  ## Print the sample descriptons
  index.text <- append(index.text, "\n\n## Sample descriptions\n")
  index.text <- append(index.text, kable(sample.desc, caption = "Sample description table"))


  samples.per.condition <- as.data.frame.table(table(sample.conditions))
  index.text <- append(index.text, "\n\n### Samples per condition\n")
  index.text <- append(index.text, kable(samples.per.condition, caption = "Samples per condition",
                                         col.names = c("Condition", "Nb samples")))
  if (verbose >= 2) { print(as.data.frame(samples.per.condition)) }
  
  ## ----read_design, warning=FALSE------------------------------------------
  # setwd(dirs["main"]) ## !!!!! I don't understand why I have to reset the working directory at each chunk

  ## Read the design file, which indicates the anlayses to be done.
  ## Each row specifies one differential expression analysis, which
  ## consists in comparing two conditions.
  message("\tReading design file: ", infiles["design"])
  design <- read.delim(file.path(dirs["main"], infiles["design"]), sep = "\t",
                       comment = c(";"), header = T, row.names = NULL)
  message("\t\tDesign file contains ", nrow(design), " comparisons. ")
  comparison.summary <- design ## Initialize a summary table for each DEG analysis
  comparison.summary$prefixes <- paste(sep = "_", design[,1], "vs", design[,2])

  ## Print out the design table (pairs of conditions to be compared)
  index.text <- append(index.text, "\n\n## Design\n")
  index.text <- append(index.text, kable(
    comparison.summary,
    row.names = TRUE,
    caption = "**Design**. Each row describes one comparison between two conditions."))



  ## ----load_count_table----------------------------------------------------
  message("\tLoading count table: ", infiles["counts"])
  ori.counts <- read.delim(infiles["counts"], row.names = 1, sep = "\t")
  # names(ori.counts)
  # dim(ori.counts)
  # View(ori.counts)

  ## Filter out the rows corresponding to non-assigned counts,
  ## e.g. __not_aligned, __ambiguous, __too_low_qAual, __not_aligned
  not.feature <- grep(rownames(ori.counts), pattern = "^__")
  if (length(not.feature) > 0) {
    all.counts <- ori.counts[-not.feature,]
  } else {
    all.counts <- ori.counts
  }
  # dim(all.counts)

  ## Just for quick test and  debug: select a random subset of features
  if (quick.test) {
    all.counts <- all.counts[sample(x = 1:nrow(all.counts), size = 1000, replace = FALSE),]
  }
  message("\t\tLoaded counts: ",
          nrow(all.counts), " features x ",
          ncol(all.counts), " samples")

  ## Check that the header of all.counts match the sample IDs
  ids.not.found <- setdiff(sample.ids, names(all.counts)) ## Identify sample IDs with no column in the count table
  if (length(ids.not.found) == length(sample.ids)) {
    colnames(all.counts) <- sample.ids
    ids.not.found <- setdiff(sample.ids, names(all.counts)) ## Identify
  } else if (length(ids.not.found) > 0) {
    stop(length(ids.not.found), " missing columns in count table\t", infiles["counts"],
         "\n\tMissing columns: ", paste(collapse = "; ", ids.not.found))
  }

  ## Restrict the count table to the sample IDs found in the sample description file
  all.counts <- all.counts[, sample.ids]
  # names(all.counts)
  # dim(all.counts)

  ##---- Feature filtering ----

  ## Load list of black-listed features
  if (is.null(parameters$DEG$blacklist)) {
    black.listed.features <- NULL
  } else {
    infiles["black_list"] <- parameters$DEG$blacklist
    black.list.table <- read.table(
      infiles["black_list"], sep = "\t",
      comment.char = "#",
      header = FALSE)
    black.listed.features <- as.vector(black.list.table[,1])
    black.listed.not.found <- setdiff(black.listed.features, row.names(all.counts))
    if (length(black.listed.not.found) > 0) {
      message("Warning: some IDs of the black list do not correspond to features of the count table")
      message("\tNot found IDs\t", paste(collapse = ", ", head(black.listed.not.found)))
    }
    black.listed.features <- intersect(black.listed.features, row.names(all.counts))
    if (length(black.listed.features) > 0) {
      message("Discarding ", length(black.listed.features), " black-listed features\t", parameters$DEG$blacklist)
    }
  }

  ## ---- Filter out features according to various user-specified criteria ----
  ## (parameters defined in the config files
  filtered.counts <- FilterCountTable(
    counts = all.counts,
    na.omit = TRUE,
    min.count = thresholds$min.count,
    mean.count = thresholds$mean.count,
    mean.per.condition = threshold$mean.per.condition,
    black.list = black.listed.features)

  ##---- Log-transform non-normalized counts ----

  ## Add an epsilon to 0 values only, in order to enable log-transform and display on logarithmic axes.
  message("\t\tTreating zero-values by adding epsilon = ", epsilon)
  filtered.counts.epsilon <- filtered.counts
  filtered.counts.epsilon[filtered.counts == 0] <- epsilon

  ## Log-transformed data for some plots.
  message("\t\tComputing log-transformed values")
  filtered.counts.log10 <- log10(filtered.counts.epsilon)
  filtered.counts.log2 <- log2(filtered.counts.epsilon)


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
  #stats.per.sample <- calc.stats.per.sample(sample.desc, all.counts)
  # View(stats.per.sample.all)
  # dim(all.counts)
  # dim(sample.desc)
  stats.per.sample.all <- cbind(
    sample.desc,
    ColStats(x = all.counts, verbose = verbose, selected.stats = selected.stats))
  #stats.per.sample.all$Mcounts <- stats.per.sample.all$sum / 1e6
  # View(stats.per.sample.all.all)
  outfiles["stats_per_sample_all_features"] <- file.path(dir.tables.samples, "stats_per_sample_all_features.tsv")
  message("\t\t", outfiles["stats_per_sample_all_features"])
  write.table(x = stats.per.sample.all, file = outfiles["stats_per_sample_all_features"], quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

  ## Compute statistics ommitting zero values
  message("\tComputing sample-wise statistics for non-zero counts")
  all.counts.nozero <- all.counts
  all.counts.nozero[all.counts.nozero == 0] <- NA
  stats.per.sample.nozero <- cbind(
    sample.desc,
    ColStats(all.counts.nozero, verbose = verbose, selected.stats = selected.stats))
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
    ColStats(filtered.counts, verbose = verbose, selected.stats = selected.stats))
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
  stdcounts.libsum <- cpm(filtered.counts.epsilon)    ## Counts per million reads, normalised by library sum
  stdcounts.perc75 <- cpm(filtered.counts.epsilon, lib.size = stats.per.sample.filtered$perc75)    ## Counts per million reads, normalised by 75th percentile
  stdcounts.perc95 <- cpm(filtered.counts.epsilon, lib.size = stats.per.sample.filtered$perc95)    ## Counts per million reads, normalised by 95th percentile
  stdcounts.median <- cpm(filtered.counts.epsilon, lib.size = stats.per.sample.filtered$median)    ## Counts per million reads, normalised by sample-wise median count

  ## Chose one of the standardization methods to get
  #stdcounts <- stdcounts.median ## Choose one normalization factor for the stdcounts used below
  stdcounts <- stdcounts.perc75 ## Choose one normalization factor for the stdcounts used below
  stdcounts.log10 <- log10(stdcounts) ## Log-10 transformed stdcounts, xwith the epsilon for 0 counts
  stdcounts.log2 <- log2(stdcounts) ## Log-10 transformed stdcounts, with the epsilon for 0 counts

  ## Export normalized counts (in log2-transformed counts per million reads)
  outfiles["stdcounts"] <- file.path(dir.tables.samples, paste(sep = "", count.prefix, "_stdcounts.tsv"))
  message("\tExporting standardized counts: ", outfiles["stdcounts"])
  write.table(x = stdcounts, row.names = TRUE, col.names = NA,
              file = outfiles["stdcounts"], sep = "\t", quote = FALSE)

  outfiles["log2stdcounts"] <- file.path(dir.tables.samples, paste(sep = "", count.prefix, "_stdcounts_log2.tsv"))
  message("\tExporting log2-transformed standardized counts: ", outfiles["log2stdcounts"])
  write.table(x = stdcounts.log2, row.names = TRUE, col.names = NA,
              file = outfiles["log2stdcounts"], sep = "\t", quote = FALSE)


  # ## Detect outliers, i.e. genes with a very high number of reads (hundreds of thousands), most of which result from problems with ribodepletion.
  # if (is.null(thresholds["max.log10.cpm"])) {
  #   outlier.threshold <- 8.5 ## Somewhat arbitrary threshold to discard
  # } else {
  #   outlier.threshold <- thresholds["max.log10.cpm"]
  # }
  # outliers <- (apply(stdcounts.log10, 1, max) > outlier.threshold)
  # message("\tDetected ", sum(outliers), " outliers with log10(stdcounts) higher than ", outlier.threshold)
  # rownames(stdcounts.log10[outliers,])
  # sum(outliers)

  ## Compute Trimmed Means of M Values (TMM): TO BE DONE
  stats.per.sample.filtered$cpm.mean <- apply(stdcounts, 2, mean)
  stats.per.sample.filtered$log2.cpm.mean <- apply(stdcounts.log2, 2, mean)
  stats.per.sample.filtered$log10.cpm.mean <- apply(stdcounts.log10, 2, mean)

  ## ---- Export stats per sample ----
  outfiles["stats_per_sample_filtered_features"] <- file.path(dir.tables.samples, "stats_per_sample_filtered_features.tsv")
  message("\t", "Exporting stats per sample for filtered features")
  message("\t\t", outfiles["stats_per_sample_filtered_features"])
  write.table(x = stats.per.sample.nozero, file = outfiles["stats_per_sample_filtered_features"], quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)


  ## ----print_sample_stats--------------------------------------------------
  ## Statistics per sample
  # setdiff(selected.stats, names(stats.per.sample))

  # names(stats.per.sample.all)
  index.text <- append(index.text, "\n\n## Sample-wise statistics\n")

  index.text <- append(index.text, "\n\n### All feature counts (zeros included)\n")
  index.text <- append(index.text, kable(stats.per.sample.all[,selected.stats], digits = 2,
                                         format.args = list(big.mark = ",", decimal.mark = "."),
                                         caption = "Sample-wise statistics for all features (zeros included)"))

  index.text <- append(index.text, "\n\n### All features, non-zero counts\n")
  index.text <- append(index.text, kable(stats.per.sample.nozero[,selected.stats], digits = 2,
                                         format.args = list(big.mark = ",", decimal.mark = "."),
                                         caption = "Sample-wise statistics for all features (zeros excluded)"))

  index.text <- append(index.text, "\n\n### Filtered features, non-zero counts\n")
  index.text <- append(index.text, kable(stats.per.sample.filtered[,selected.stats], digits = 2,
                                         format.args = list(big.mark = ",", decimal.mark = "."),
                                         caption = "Sample-wise statistics for filtered features"))


  ## ----library_sizes_barplot, fig.width=6, fig.height=6, fig.cap="**Barplot of assigned reads per sample. ** Bars indicate the sum of read counts assigned to features (genes) per sample (library)."----
  par.ori <- par(no.readonly = TRUE) # Store original parameters
  # par(c("mar", "mai"))
  # libsize.barplot(
  #   stats.per.sample,
  #   plot.file = NULL,
  #   main = "Assigned reads per sample (libsum)")
  # par(par.ori) # Restore original parameters
  #
  figname <- "libsize_barplot_all_features"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)
  for (fig.format in figure.formats) {
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 8, height = 8)

    LibsizeBarplot(counts = all.counts, sample.labels = sample.labels, sample.colors = sample.desc$color, main = "All features", cex.axis = 0.8)

    silence <- dev.off(); rm(silence)
    index.text <- index.figure(figname, figure.file, index.text)
    # system(paste("open", figure.file))
  }


  figname <- "libsize_barplot_filtered_features"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)
  for (fig.format in figure.formats) {
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 7, height = 8)

    LibsizeBarplot(counts = filtered.counts, sample.labels = sample.labels, sample.colors = sample.desc$color, main = "After filtering", cex.names = 0.8)

    silence <- dev.off(); rm(silence)
    index.text <- index.figure(figname, figure.file, index.text)
    # system(paste("open", figure.file))
  }



  ## ---- Between-sample normalisation -------------------------------------------------------

  if (is.null(parameters$DEG$norm_method)) {
    #norm.methods <- c("none", "mean", "median", "percentile", "TMM", "DESeq2", "quantiles")
    norm.methods <- c("none", "mean", "median", "percentile", "TMM", "DESeq2")
  } else {
    norm.methods <- parameters$DEG$norm_method
  }
  norm.comparison <- NormalizeCountTable(
    counts = filtered.counts, class.labels = sample.conditions, nozero = TRUE,
    method = norm.methods, percentile = 75, log2 = FALSE, epsilon = 0.1, detailed.sample.stats = TRUE,
    verbose = verbose)
  #names(norm.comparison)

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
  for (fig.format in figure.formats) {
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    # message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 7, height = 8)

    plot(size.factors, main = "Sample size factors", col = sample.desc$color)

    silence <- dev.off(); rm(silence)
    index.text <- index.figure(figname, figure.file, index.text)
    # system(paste("open", figure.file))
  }


  ## ----differential_expression_analysis, fig.width=8, fig.height=12--------
  # setwd(dirs["main"]) ## !!!!! I don't understand why I have to reset the working directory at each chunk

  # i <- 1  # for quick test
  for (i in 1:nrow(design)) {


    prefix <- list() ## list for output file prefixes

    deg.results <- list()

    ## Identify samples for the first condition
    cond1 <- as.vector(design[i,1])  ## First condition for the current comparison
    samples1 <- sample.ids[sample.conditions == cond1]
    if (length(samples1) < 2) {
      stop(paste("Cannot perform differential analysis. The count table contains less than 2 samples for condition", cond1))
    }

    ## Identify samples for the second condition
    cond2 <- as.vector(design[i,2])  ## Second condition for the current comparison
    samples2 <- sample.ids[sample.conditions == cond2]
    if (length(samples2) < 2) {
      stop(paste("Cannot perform differential analysis. The count table contains less than 2 samples for condition", cond2))
    }

    #  stop("HELLO", "\tprefix = ", prefix)

    message("\tDifferential analysis\t", i , "/", nrow(design), "\t", cond1, " vs ", cond2)

    ## Create a specific directory for the results of this comparison
    comparison.prefix <- comparison.summary$prefixes[i]
    dir.results <- file.path(result.dir, paste(sep = "", comparison.prefix))
    dirs[comparison.prefix] <- dir.results ## Index current result dir for the list of directories
    comparison.summary[i, "result.dir"] <- dir.results ## Include current result dir to the comparison summary table
    dir.create(path = file.path(dirs["main"], dir.results), showWarnings = FALSE, recursive = TRUE)
    prefix["comparison_file"] <- file.path(dir.results, comparison.prefix)
    message("\t\tresults:\t", dir.results)

    ## Create a specific directory for the figures of this comparison
    dir.figures.diffexpr <-  file.path(dirs[comparison.prefix], "figures")
    dirs[paste(sep = "_", comparison.prefix, "figures")] <- dir.figures.diffexpr ## Index current figures dir for the list of directories
    comparison.summary[i, "figures"] <- dir.figures.diffexpr ## Include current figures dir to the comparison summary table
    dir.create(path = file.path(dirs["main"], dir.figures.diffexpr), showWarnings = FALSE, recursive = TRUE)
    prefix["comparison_figure"] <- file.path(dir.figures.diffexpr, comparison.prefix)
    message("\t\tfigures:\t", dir.figures.diffexpr)
    #    paste(sep = "", comparison.prefix, "_",  suffix.deg))


    ## Select counts for the samples belonging to the two conditions
    current.samples <- c(samples1, samples2)
    current.counts <- data.frame(filtered.counts[,current.samples])
    # dim(current.counts)  ## For test
    # names(current.counts)

    if (sum(!names(current.counts) %in% sample.ids) > 0) {
      stop("Count table contains column names without ID in sample description file.")
    }

    ## Define conditions and labels for the samples of the current analysis
    current.conditions <- sample.conditions[current.samples]
    current.labels <- paste(current.conditions, names(current.counts), sep = "_")

    result.table <- init.deg.table(stdcounts, samples1, samples2, stats = FALSE)
    # View(result.table)
    # dim(result.table)

    ## ---- DESeq2 analysis ----
    message("\tDESeq2 analysis\t", comparison.prefix)
    deseq2.result <- deseq2.analysis(
      counts = current.counts,
      condition = current.conditions,
      comparison.prefix = comparison.prefix,
      ref.condition = cond2,
      thresholds = as.vector(thresholds),
      title = comparison.prefix,
      dir.figures = dir.figures.diffexpr, 
      verbose = verbose)
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
    # names(result.table)
    # dim(deseq2.result$result.table)
    # dim(deseq2.result$result.table)
    # names(deseq2.result$result.table)
    # View(deseq2.result$result.table)
    # View(result.table)


    ## Save the completed DESeq2 result table
    deseq2.result.file <- paste(sep = "", prefix["comparison_file"], "_DESeq2.tsv")
    comparison.summary[i,"deseq2"] <- deseq2.result.file
    message("\tExporting DESeq2 result table (tab): ", deseq2.result.file)
    write.table(
      x = deseq2.result$result.table, row.name    = FALSE,
      file = deseq2.result.file,
      sep = "\t", quote = FALSE)

    ## ---- edgeR analysis ----
    # norm.method <- "TMM" ## For quick test and debugging
    for (norm.method in edgeR.norm.methods) {

      edgeR.prefix <- paste(sep = "_", "edgeR", norm.method)

      edger.result <- edger.analysis(
        counts = current.counts,
        condition = current.conditions,
        test.condition = cond1,
        ref.condition = cond2,
        thresholds = as.vector(thresholds),
        comparison.prefix = comparison.prefix,
        norm.method = norm.method,
        title = paste(sep = "_", norm.method, comparison.prefix),
        dir.figures = dir.figures.diffexpr,
        verbose = verbose)
      deg.results[[edgeR.prefix]] <- edger.result

      ## A tricky way to add edgeR with normalisation in column names
      edger.to.bind <- edger.result$result.table[row.names(result.table),]
      colnames(edger.to.bind) <- paste(sep = "_", edgeR.prefix, colnames(edger.to.bind))
      # names(edger.to.bind)
      # View(x)
      # x <- rownames(result.table)
      # y <- rownames(edger.to.bind)
      # sum(x != y)
      # names (result.table)
      result.table <- cbind(
        result.table,
        edger.to.bind)
      # names (result.table)


      ## Export edgeR result table
      edger.result.file <- paste(sep = "", prefix["comparison_file"], "_", edgeR.prefix, ".tsv")
      comparison.summary[i,"edger"] <- edger.result.file
      message("\tExporting edgeR result table (tab): ", edger.result.file)
      write.table(x = edger.result$result.table,
                  file = edger.result.file,
                  row.names = FALSE,
                  sep = "\t", quote = FALSE)
    }

    ## Export full result table (DESeq2 + edgeR with different normalisation methods)
    ## in a tab-separated values (tsv) file
    result.file <- paste(sep = "",
                         prefix["comparison_file"],
                         "_diffexpr_DESeq2_and_edgeR.tsv")
    # comparison.summary[i,"result.table"] <- paste(sep=".", result.file, "tsv")
    message("\tExporting result table (tsv): ", result.file)
    write.table(x = result.table, row.names = FALSE,
                file = result.file, sep = "\t", quote = FALSE)


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
    # feature.scores <- log2(apply(filtered.counts.epsilon, 1, median))

    # hist(feature.scores, breaks = 100)

    # View(deg.compa$padj)
    ## compare DESeq2 and edgeR normalisatio results

    ## ---- Plot comparing DEGs obtained with DESeq2 and edgeR ----

    ## Comparison between adjusted p-values
    figname <- paste(sep = "", comparison.prefix, "_norm_compa_padj")
    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (fig.format in figure.formats) {
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
      index.text <- index.figure(figname, figure.file, index.text)
      # system(paste("open", figure.file))
    }

    figname <- paste(sep = "", comparison.prefix, "_norm_compa_log2FC")
    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (fig.format in figure.formats) {
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
      index.text <- index.figure(figname, figure.file, index.text)
      # system(paste("open", figure.file))
    }


    ## ---- Draw Volcano plots -----
    deg.names <- names(deg.results)
    nb.panels <- n2mfrow(length(deg.names))
    if (length(deg.names) <= 5) {
      nb.panels <- rev(nb.panels)
    }
    figname <- paste(sep = "", comparison.prefix, "_norm_compa_volcano_plots")

    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (fig.format in figure.formats) {
      figure.file <- paste(sep = "", file.prefix, ".", fig.format)
      figure.files[[fig.format]][figname] <- figure.file
      # message("\t\t\t", fig.format, " plot\t", figname)
      OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, height = 2 + nb.panels[1]*3, width = 2 + nb.panels[2]*3.5)

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
      index.text <- index.figure(figname, figure.file, index.text)
      # system(paste("open", figure.file))
    }


    ## system(paste("ls -ltr ", pdf.files[prefix]))



    ## ---- Comparison between p-value histograms -----
    ## Draw Volcano plots
    # deg.name <- "DESeq2"
    # deg.name <- "edgeR_TMM"
    figname <- paste(sep = "", comparison.prefix, "_norm_compa_pvalue_histograms")

    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (fig.format in figure.formats) {
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
      index.text <- index.figure(figname, figure.file, index.text)
      # system(paste("open", figure.file))
    }


    ## ---- Draw Venn diagram with number of genes declared significant ----
    ## according to the selection criteria (threshold fields).
    selection.fields <- c("padj", "FC")
    selection.thresholds <- thresholds[selection.fields]
    selection.columns <- paste(sep = "", selection.fields, "_", selection.thresholds)

    figname <- paste(sep = "", comparison.prefix, "_norm_compa_Venn_", paste(collapse = "_", selection.fields))

    file.prefix <- file.path(dir.figures.diffexpr, figname)
    message("\t\tGenerating figure\t", figname)
    for (fig.format in figure.formats) {
      figure.file <- paste(sep = "", file.prefix, ".", fig.format)
      figure.files[[fig.format]][figname] <- figure.file
      # message("\t\t\t", fig.format, " plot\t", figname)
      OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, height = 2 + nb.panels[1]*3, width = 2 + nb.panels[2]*3.5)

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
      index.text <- index.figure(figname, figure.file, index.text)
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

  index.text <- append(index.text, "\n\n## Directories and files\n")

  index.text <- append(index.text, paste(sep = "", "- Config file : ", "[", configFile, "](", configFile,")"))
  index.text <- append(index.text, paste(sep = "", "- Count table : ", "[", count.table, "](", count.table,")"))
  index.text <- append(index.text, paste(sep = "", "- Sample descriptions : ", "[", infiles["sample descriptions"], "](", infiles["sample descriptions"], ")"))
  index.text <- append(index.text, paste(sep = "", "- Design : ", "[", infiles["design"], "](", infiles["design"], ")"))

  index.text <- append(index.text, paste(sep = "", "\n\n### Directories\n"))
  index.text <- append(index.text, paste(sep = "", "| Content | Path |"))
  index.text <- append(index.text, paste(sep = "", "| ---------------------- | -------------------------------------------------- |"))
  for (dirname in names(dirs)) {
    dir <- dirs[dirname]
    index.text <- append(index.text, paste(sep = "", "| ", dirname, " | ", "[", dir, "](", dir, ") |"))
  }

  index.text <- append(index.text, paste(sep = "", "\n\n### Input files\n"))
  index.text <- append(index.text, paste(sep = "", "| Content | Path |"))
  index.text <- append(index.text, paste(sep = "", "| ---------------------- | -------------------------------------------------- |"))
  for (filename in names(infiles)) {
    file <- infiles[filename]
    index.text <- append(index.text, paste(sep = "", "| ", filename, " | ", "[", file, "](", file, ") |"))
  }

  index.text <- append(index.text, paste(sep = "", "\n\n### Output files"))
  index.text <- append(index.text, paste(sep = "", "| Content | Path |"))
  index.text <- append(index.text, paste(sep = "", "| ---------------------- | -------------------------------------------------- |"))
  for (filename in names(outfiles)) {
    file <- outfiles[filename]
    index.text <- append(index.text, paste(sep = "", "| ", filename, " | ", "[", file, "](", file, ") |"))
  }

  for (fig.format in names(figure.files)) {
    figfiles <- figure.files[fig.format]
    index.text <- append(index.text, paste(sep = "", "\n\n### Figures (", fig.format,")"))
    index.text <- append(index.text, paste(sep = "", "| Content | Path |"))
    index.text <- append(index.text, paste(sep = "", "| ---------------------- | -------------------------------------------------- |"))
    for (filename in names(figfiles)) {
      file <- figfiles[filename]
      index.text <- append(index.text, paste(sep = "", "| ", filename, " | ", "[", file, "](", file, ") |"))
    }
  }

  ## ---- Session info ---------------------------------------------------------
  ## Print the complete list of libraries + versions used in this session
  index.text <- append(index.text, "\n\n## Session info\n")
  session.info <- sessionInfo()
  # print(session.info)
  index.text <- append(index.text, "\n```\n")
  index.text <- append(index.text, capture.output(session.info))
  index.text <- append(index.text, "\n```\n")



  ## ---- Generate the HTML index -----
  index.text <- append(index.text, paste(sep = "", "\n\nJob done: ", Sys.time()))

  writeLines(text = index.text, con = index.socket)
  close(index.socket)
  rmarkdown::render(index.Rmd, output_format = "html_document")
  rmarkdown::render(index.Rmd, output_format = "pdf_document")


  ## ---- job_done ------------------------------------------------------------
  message("Tables directory\t", dir.tables.samples)
  message("Figures directory\t", dir.figures.samples)
  message("Index of input/output files\t", index.file)
  message("Job done")

}
