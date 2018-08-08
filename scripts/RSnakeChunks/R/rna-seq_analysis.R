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
#' @param result.dir="results" directory where the results will be stored. Should be defined relative to main directory.
#' @param count.prefix=NULL basename to export the transformed count tables. If not specified, computed automatically from the name of the count table file.
#' @param edgeR.norm.methods="TMM" vector with one or more normalisation methods for edegR. Supported: TMM, RLE, upperquartile, quantiles, none
#' @param verbose=1 level of verbosity
#'
#' @export
RNAseqAnalysis <- function(count.table,
                           configFile,
                           main.dir = getwd(),
                           result.dir = "results",
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
  ## Problem with relevel : deseq2 and edgeR compute log2-ratios in opposite directions !

  
  ## ---- Load required libraries ----
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
  infiles <- vector()   ## Input files
  outfiles <- vector()  ## For tab-separated value files
  
  
  ## ---- Load configuration file (YAML-formatted) ----
  infiles["config"] <- configFile
  message("\tConfiguration file: ", configFile)
  if (!exists("configFile")) {
    ## The prompt does not seem to work with the Rmd documents
    #   message("Choose the parameter file")
    #   parameter.file <- file.choose()
    stop("This report requires to specify a variable named configFile, containing the path to an YAML-formatted file describing the parameters for this analysis.")
  }
  parameters <- yaml.load_file(configFile)
  message("\tLoaded parameters from file ", configFile)
  # View(parameters)
  
  
  ## ---- Define the main parameters -----------------------------------------------------


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

  
  ## Directory to store differential expression results
  dirs["output"] <- result.dir ## Index dir for the report
  message("\tOutput directory: ", result.dir)
  dir.create(result.dir, showWarnings = FALSE, recursive = TRUE)
  
  ## Prefix for the count table
  #count.prefix <- "bowtie2_featureCounts_all"
  #message("\tPrefix for the count table: ", count.prefix)
  infiles["count_table"] <- count.table
  message("\tCount table: ", count.table)

 
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

  ## Normalisation method(s) for edgeR
  if (is.null(parameters$edgeR$norm_method)) {
    # If not defined in config file, use edgeR norm methods specified in the arguments
    parameters$edgeR$norm_method <- edgeR.norm.methods
  } else {
    edgeR.norm.methods <- parameters$edgeR$norm_method
  }

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

  ## ---- Set some default parameters ----
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
  
  ## A trick: to enable log-scaled plots for 0 values, I add an epsilon increment
  if (is.null(parameters$DEG$out_width)) {
    parameters$DEG$out_width <- "75%" # passed to file parameters.R, 2017-03-15
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
  infiles["counts"] <- count.table
  if (!file.exists(count.table)) {
    stop("Feature count table does not exist: ", count.table)
  } else {
    message("\tFeature count table: ", count.table)
  }


  ## ---- Define thresholds ----
  if (is.null(parameters$DEG$thresholds)) {
    message("\tDEG thresholds were not defined in config file -> using default values")
    if (is.null(parameters$DEG)) {
      parameters$DEG <- list()
    }
    parameters$DEG$thresholds <- list(
      min.count = 1,
      mean.count = 5,
      padj = 0.05,
      FC = 2)
  }
  thresholds <- as.vector(parameters$DEG$thresholds)

  ## ---- Initialize the Rmd report (index of input/output file) ----
  index.Rmd <- "index.Rmd"
  index.socket <- file(index.Rmd)

  
  ## Get some elements from the  config file (if defined)  
  if (is.null(parameters$title)) {
    parameters$title <- "RNA-seq analysis report"
  }
  if (is.null(parameters$description)) {
    parameters$description <- "Analysis workflow for RNA-seq data: sample normalisation, descriptive statistics, differential analysis. "
  }
  if (is.null(parameters$author)) {
    parameters$author <- "[AUTHORS] "
  }
  
  Rmd.header <- paste(
    sep = '', 
    '---
title:  "', parameters$title,'"
author: "', parameters$author,'"
date: Last update:`r format(Sys.time())`
output:
  html_document:
    fig_caption: yes
    highlight: kate
    self_contained: yes
    theme: cerulean
    toc: yes
    toc_depth: 3
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
  
  index.text <- Rmd.header
  
  ## Print the description
  if (!is.null(parameters$description)) {
    index.text <- append(index.text, "\n\n## Description\n")
    index.text <- append(index.text, parameters$description)    
  }
  
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

  ## ---- Read sample description table ----
  
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
  

  ## ---- Extract conditions from the sample table
  
  ## Identify if the sample description file contains a column with heading "condition"(case-insensitive)
  condition.column <- grep(pattern = "condition", 
            x = colnames(sample.desc), 
            ignore.case = TRUE)
  if (length(condition.column) != 1) {
    ## If not column contains "condition" in the sample descriptions, 
    ## use the first column
    condition.column <- 1
  }
  message("\t\tSample conditions in column ", condition.column, " of sample description table. ")
  sample.conditions <- as.vector(sample.desc[,condition.column]) ## Condition associated to each sample
  names(sample.conditions) <- sample.ids
  
  ## Unique conditions
  conditions <- unique(sample.conditions) ## Set of distinct conditions
  message("\t\tNb of conditions = ", length(conditions))
  message("\t\tConditions = ", paste(collapse = ", ", conditions))
  
  ## ---- Define condition-specific colors and apply them to each sample ----
  color.per.condition <- brewer.pal(max(3, length(conditions)),"Dark2")[1:length(conditions)]
  names(color.per.condition) <- conditions
  sample.desc$color <- color.per.condition[sample.conditions]

  # print(sample.conditions)

  ## ---- Build sample labels by concatenating their ID and condition ----
  labeL.column <- grep(pattern = "label", x = colnames(sample.desc), ignore.case = TRUE)
  if (length(labeL.column) == 1) {
    message("\t\tSample labels in column ", labeL.column, " of sample description table. ")
    sample.labels <- as.vector(unlist(sample.desc[, labeL.column]))
  } else {
    message("\t\tSample description file has no column with heading 'Label'. ")
    sample.labels <- paste(sep = "_", sample.ids, sample.conditions)
    if (is.null(sample.desc$Replicate)) {
      message("\t\tBuilding labels from sample IDs and conditions. ")
    } else {
      message("\t\tBuilding labels from sample IDs, conditions and replicate nb. ")
      sample.labels <- paste(sep = "_", sample.labels, sample.desc$Replicate)
    }
    sample.desc$Label <- sample.labels
  }
  # print(sample.labels)





  ## Print the sample descriptons
  index.text <- append(index.text, "\n\n## Sample descriptions\n")
  index.text <- append(index.text, kable(sample.desc, caption = "Sample description table"))


  samples.per.condition <- as.data.frame.table(table(sample.conditions))
  index.text <- append(index.text, "\n\n### Samples per condition\n")
  index.text <- append(index.text, kable(samples.per.condition, caption = "Samples per condition",
                                         col.names = c("Condition", "Nb samples")))
  if (verbose >= 2) { print(as.data.frame(samples.per.condition)) }
  
  ## ---- Read design table ----
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
  
  ## Identify reference column in the design file
  if ("reference" %in% (tolower(colnames(design)))) {
    reference.column <- which(tolower(colnames(design)) == "reference")
    if (length (reference.column) != 1) {
      stop("The design file contains several columns entitled 'reference'. ")
    }
  } else {
    message("The headers of the design table do not contain 'reference' column -> taking 1st column as reference")
    reference.column <- 1
  }

  ## Identify test column in the design file
  if ("test" %in% (tolower(colnames(design)))) {
    test.column <- which(tolower(colnames(design)) == "test")
    if (length (test.column) != 1) {
      stop("The design file contains several columns entitled 'test'. ")
    }
  } else {
    message("The headers of the design table do not contain 'test' column -> taking 1st column as test")
    test.column <- 2
  }
  
  ## Print out the design table (pairs of conditions to be compared)
  index.text <- append(index.text, "\n\n## Design\n")
  index.text <- append(index.text, kable(
    comparison.summary,
    row.names = TRUE,
    caption = "**Design**. Each row describes one comparison between two conditions."))



  ## ---- Load count table ----
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
  figname <- "libsize_barplot_all_features"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)
  for (f in 1:length(figure.formats)) {
    fig.format <- figure.formats[f]
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 8, height = 8)

    LibsizeBarplot(counts = all.counts, sample.labels = sample.labels, sample.colors = sample.desc$color, main = "All features", cex.axis = 0.8)

    silence <- dev.off(); rm(silence)
    if (f == 1) {
      index.text <- index.figure(figname, figure.file, index.text, out.width = parameters$DEG$out_width)
    }
    # system(paste("open", figure.file))
  }


  figname <- "libsize_barplot_filtered_features"
  file.prefix <- file.path(dir.figures.samples, figname)
  message("\t\tGenerating figure\t", figname)
  for (f in 1:length(figure.formats)) {
    fig.format <- figure.formats[f]
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 7, height = 8)

    LibsizeBarplot(counts = filtered.counts, sample.labels = sample.labels, sample.colors = sample.desc$color, main = "After filtering", cex.names = 0.8)

    silence <- dev.off(); rm(silence)
    if (f == 1) {
      index.text <- index.figure(figname, figure.file, index.text, out.width = parameters$DEG$out_width)
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
    plot.height <- max(2 + 0.4 * nrow(sample.desc), 6)
    
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
    par(mar=c(5.1, left.margin, 4.1, 1))
    
    x <- x[nrow(x):1,]
    indices <- sort(rep(1:nrow(x), times = 2))
    barplot.densities <- rep(c(50, -1), length.out = 2 * nrow(x))
    bplt <- barplot(
      as.matrix(t(x)), 
      horiz = TRUE, 
      beside = TRUE, 
      las=1,
      cex.names = 0.9,
      col = rev(sample.desc$color[indices]),
      density = barplot.densities, xlab = "Million reads", 
      main = "Library sizes\nbefore/after filtering",
      xlim = c(0, max(x) * 1.3))
    barplot.Mreads <- round(digits=1, as.vector(unlist(t(x))))
    text(x=pmax(barplot.Mreads, 3), labels =  barplot.Mreads, y = bplt, pos = 2, font = 2, cex=0.8)
    #sample.type <- rep(c("Filtered", "All"), length.out = 2 * nrow(x))
    #text(x = 0, labels = sample.type, y = bplt, pos = 4, font = 1, cex=0.5)
    legend("topright", legend = c("All genes", "Filtered"), col = "gray", density = c(-1, 50), cex = 0.7)
    
    par(par.ori)
    
    silence <- dev.off(); rm(silence)
    if (f == 1) {
      index.text <- index.figure(
        figname, figure.file, index.text, out.width = parameters$DEG$out_width)
      #, chunk.opt = paste(sep = "", ", out.height = ", plot.height))
    }
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
  for (f in 1:length(figure.formats)) {
    fig.format <- figure.formats[f]
    figure.file <- paste(sep = "", file.prefix, ".", fig.format)
    figure.files[[fig.format]][figname] <- figure.file
    # message("\t\t\t", fig.format, " plot\t", figname)
    OpenPlotDevice(file.prefix = file.prefix, fig.format = fig.format, width = 7, height = 8)

    plot(size.factors, main = "Sample size factors", col = sample.desc$color)

    silence <- dev.off(); rm(silence)
    if (f == 1) {
      index.text <- index.figure(figname, figure.file, index.text, out.width = parameters$DEG$out_width)
    }
    # system(paste("open", figure.file))
  }


  ## ----differential_expression_analysis, fig.width=8, fig.height=12--------
  # setwd(dirs["main"]) ## !!!!! I don't understand why I have to reset the working directory at each chunk

  index.text <- append(index.text, "\n\n## Differential analysis\n")
  
  
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


    index.text <- append(index.text, paste(sep= "", "\n\n### ", test.condition, " versus ", ref.condition, "\n\n"))
#    index.text <- append(index.text, paste(sep= "", "Test condition : ",  test.condition, "\n"))
#    index.text <- append(index.text, paste(sep= "", "Reference condition :",  ref.condition, "\n"))
    compa.table <- data.frame(
      type = c("Reference", "Test"),
      condition = c(ref.condition, test.condition),
      samples = c(length(ref.samples), length(test.samples))
    )
    index.text <- append(index.text, 
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
    current.counts <- data.frame(filtered.counts[,current.samples])
    # dim(current.counts)  ## For test
    # names(current.counts)

    if (sum(!names(current.counts) %in% sample.ids) > 0) {
      stop("Count table contains column names without ID in sample description file.")
    }

    ## Define conditions and labels for the samples of the current analysis
    current.conditions <- sample.conditions[current.samples]
    current.labels <- paste(current.conditions, names(current.counts), sep = "_")

    ## A big result table with all features and all statistics
    result.table <- init.deg.table(stdcounts, ref.samples, test.samples, stats = FALSE)
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
    outfiles[paste(sep="_", comparison.prefix, "DESeq2")]  <- deseq2.result.file
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
    # feature.scores <- log2(apply(filtered.counts.epsilon, 1, median))

    # hist(feature.scores, breaks = 100)

    # View(deg.compa$padj)
    ## compare DESeq2 and edgeR normalisatio results

    ## ---- Plot comparing DEGs obtained with DESeq2 and edgeR ----

    ## Comparison between adjusted p-values
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
        index.text <- index.figure(figname, figure.file, index.text, out.width = parameters$DEG$out_width)
      }
      # system(paste("open", figure.file))
    }

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
        index.text <- index.figure(figname, figure.file, index.text, out.width = parameters$DEG$out_width)
      }
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
        index.text <- index.figure(figname, figure.file, index.text, out.width = parameters$DEG$out_width)
      }
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
        index.text <- index.figure(figname, figure.file, index.text, out.width = parameters$DEG$out_width)
      }
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
        index.text <- index.figure(figname, figure.file, index.text, out.width = parameters$DEG$out_width)
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

  index.text <- append(index.text, "\n\n## Directories and files\n")

  index.text <- append(index.text, paste(sep = "", "- Config file : ", "[", configFile, "](", configFile,")"))
  index.text <- append(index.text, paste(sep = "", "- Count table : ", "[", count.table, "](", count.table,")"))
  index.text <- append(index.text, paste(sep = "", "- Sample descriptions : ", "[", infiles["sample descriptions"], "](", infiles["sample descriptions"], ")"))
  index.text <- append(index.text, paste(sep = "", "- Design : ", "[", infiles["design"], "](", infiles["design"], ")"))

  index.text <- append(index.text, paste(sep = "", "\n\n### Directories\n"))
  index.text <- append(index.text, paste(sep = "", "| Content | Path |"))
  index.text <- append(index.text, paste(sep = "", "|:----------------------|--------------------------------------------------|"))
  for (dirname in names(dirs)) {
    dir <- dirs[dirname]
    index.text <- append(index.text, paste(sep = "", "| ", dirname, " | ", "[", dir, "](", dir, ") |"))
  }

  index.text <- append(index.text, paste(sep = "", "\n\n### Input files\n"))
  index.text <- append(index.text, paste(sep = "", "| Content | Path |"))
  index.text <- append(index.text, paste(sep = "", "|:----------------------|--------------------------------------------------|"))
  for (filename in names(infiles)) {
    file <- infiles[filename]
    index.text <- append(index.text, paste(sep = "", "| ", filename, " | ", "[", file, "](", file, ") |"))
  }

  index.text <- append(index.text, paste(sep = "", "\n\n### Output files"))
  index.text <- append(index.text, paste(sep = "", "| Content | Path |"))
  index.text <- append(index.text, paste(sep = "", "|:----------------------|--------------------------------------------------|"))
  for (filename in names(outfiles)) {
    file <- outfiles[filename]
    index.text <- append(index.text, paste(sep = "", "| ", filename, " | ", "[", file, "](", file, ") |"))
  }

  f <- 1
  for (f in 1:length(figure.formats)) {
    fig.format <- figure.formats[f]
    figfiles <- unlist(figure.files[fig.format])
    index.text <- append(index.text, paste(sep = "", "\n\n### Figures (", fig.format,")"))
    index.text <- append(index.text, paste(sep = "", "| Content | Path |"))
    index.text <- append(index.text, paste(sep = "", "|----------------------|--------------------------------------------------|"))
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
