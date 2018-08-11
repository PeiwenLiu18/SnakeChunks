#' @title Load RNA-seq data, metadata and parameters
#' @author Jacques van Helden
#' @description Load an RNA-seq dataset (i.e. a table with raw counts per feature), 
#' together with sample descriptions, design for differential expression, and parameters
#' from  yaml formatted file. 
#' @param countFile a tab-separated-value file containing the counts of reads for each feature (row) in each sample (column). 
#' @param configFile a yaml-formatted file containing the parameters, including the path to metadata files (sample description table, design table). 
#' @param verbose=1 level or verbosity
#' 
#' @return a list with all the required data for further analysis
#' \itemize{
#' \item countTable the count table (data.frame) with one row per feature and one column per sample
#' \item sampleDescriptions a data.frame with one row per sample and one column per attribute.
#' \item design a data.frame with one row per analysis, indicating the reference and test conditions
#' \item parameters a list with the parameters read from the yaml-formatted configuration file
#' }
#' @export
LoadRNAseqDataset <- function(countFile, 
                              configFile,
                              verbose = 1) {
  
  
  
  ## ---- Load count table ----
  if (!file.exists(countFile)) {
    stop("Count table does not exist: ", countFile)
  }
  if (verbose >= 1) { message("\tLoading count table: ", countFile) }
  counts <- read.delim(countFile, row.names = 1, sep = "\t")
  if (verbose >= 1) { message("\t\tCount table dimensions: ", 
                              nrow(counts), " rows (features) x ", 
                              ncol(counts), "columns (samples). ") }
  
  
  
  ## ---- Load configuration file (YAML-formatted) ----
  if (!exists("configFile")) {
    ## The prompt does not seem to work with the Rmd documents
    #   message("Choose the parameter file")
    #   parameter.file <- file.choose()
    stop("This report requires to specify a variable named configFile, containing the path to an YAML-formatted file describing the parameters for this analysis.")
  }
  if (verbose >= 1) { message("\tLoading configuration file: ", configFile) }
  parameters <- yaml.load_file(configFile)
  
  
  ## --- Check parameters ----
  
  
  ## Sample description file
  if (is.null(parameters$metadata$samples)) {
    stop("The sample file must be defined in the metadata seection of the yaml config file: ", configFile)
  } 
  
  ## Design file
  if (is.null(parameters$metadata$design)) {
    stop("The design file must be defined in the metadata seection of the yaml config file: ", configFile)
  } 
  
  
  ## A trick: to enable log-scaled plots for 0 values, I add an epsilon increment
  if (is.null(parameters$DEG$epsilon)) {
    parameters$DEG$epsilon <- 0.1 # passed to file parameters.R, 2017-03-15
  }
  
  ## Figures in different format
  if (is.null(parameters$figure_formats)) {
    parameters$figure_formats <- c("png", "pdf")
  }  
  
  ## Output width parameter for figures in the Rmd report
  if (is.null(parameters$DEG$out_width)) {
    parameters$DEG$out_width <- "75%" # passed to file parameters.R, 2017-03-15
  }
  
  ## Default method for the selection of the final list of DEG
  ## Since DESeq2 seems more conservative, we use it by default
  if (is.null(parameters$DEG$selection_criterion)) {
    parameters$DEG$selection_criterion <- "DESeq2"
  }
  # DEG.selection.criterion <- parameters$DEG$selection_criterion
  
  
  
  ## ---- Normalisation methods ----
  if (is.null(parameters$DEG$norm_method)) {
    parameters$DEG$norm_method <- c("DESeq2", "TMM")
  }
  
  ## Normalisation method(s) for edgeR differential analysis.
  edgeR.supported.norm.methods <- c("TMM", "RLE", "upperquartile", "none")
  if (is.null(parameters$edgeR$norm_method)) {
    # If not defined in config file, use the norm methods specified in the arguments that are supporte by edgeR
    parameters$edgeR$norm_method <- intersect(norm.methods, edgeR.supported.norm.methods)
  } 
  ## Take edgeR normalisation methods from config file
  edgeR.norm.methods <- parameters$edgeR$norm_method
  ## Check that edgeR norm methods provided in config files are all supported by edgeR
  edgeR.not.supported <- setdiff(edgeR.norm.methods, edgeR.supported.norm.methods)
  if (length(edgeR.not.supported) > 0) {
    stop("Invalid normalisation method for edgeR: ",
         paste(collapse = ", ", edgeR.not.supported),
         ". Supported: ", paste(collapse = ",", edgeR.supported.norm.methods), ".")
  }

  ## ---- Normalisation percentiles ----
  if (is.null(parameters$DEG$norm_percentile)) {
    parameters$DEG$norm_percentile <- 0.75
  }
  
  ## Convert our percentile parameter into edgeR p parameter
  ## (comprized between 0 and 1, which is actually a quantile rather than a percentile).
  parameters$edgeR$norm.p <- parameters$DEG$norm_percentile/100
  
  ## ---- DEG selection thresholds ----
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
  
  ## Get some elements from the  config file (if defined)
  if (is.null(parameters$title)) {
    parameters$title <- "RNA-seq analysis report"
  }
  if (is.null(parameters$description)) {
    parameters$description <- "Analysis workflow for RNA-seq data: sample normalisation, descriptive statistics, differential analysis. "
  }
  if (is.null(parameters$author)) {
    parameters$author <- "[AUTHOR] "
  }
  if (is.null(parameters$author_email)) {
    parameters$author_email <- "[AUTHOR_EMAIL] "
  }
  
  
  ## ---- Load sample description table ----
  
  ## Read the sample description file, which indicates the
  ## condition associated to each sample ID.
  if (!file.exists(parameters$metadata$samples)) {
    stop("Sample description file does not exist\t", parameters$metadata$samples)
  }
  message("\tLoading sample description file: ", parameters$metadata$samples)
  sample.desc <- read.delim(
    file = parameters$metadata$samples,
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
  
  ## Create a table with attributes per condition
  condition.table <- data.frame(row.names = conditions,
                                color = brewer.pal(max(3, length(conditions)),"Dark2")[1:length(conditions)]
  )
  
  ## Define sample colors according to their condition
  color.per.condition <- condition.table$color
  names(color.per.condition) <- row.names(condition.table)
  sample.desc$color <- color.per.condition[sample.conditions]
  # print(sample.conditions)
  
  ## ---- Build sample labels by concatenating their ID and condition ----
  labeL.column <- which(tolower(colnames(sample.desc)) == "label")
  #    grep(pattern = "label", x = colnames(sample.desc), ignore.case = TRUE)
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
  
  ## ---- Read design table ----
  
  ## Read the design table, which indicates the anlayses to be done.
  ## Each row specifies one differential expression analysis, which
  ## consists in comparing two conditions.
  message("\tReading design file: ", parameters$metadata$design)
  design <- read.delim(parameters$metadata$design, sep = "\t",
                       comment = c(";"), header = T, row.names = NULL)
  message("\t\tDesign file contains ", nrow(design), " comparisons. ")
  comparison.summary <- design ## Initialize a summary table for each DEG analysis
  comparison.summary$prefixes <- paste(sep = "_", design[,1], "vs", design[,2])
  
  ## Identify reference column in the design file
  if ("reference" %in% (tolower(colnames(design)))) {
    reference.column <- which(tolower(colnames(design)) == "reference")
    if (length(reference.column) != 1) {
      stop("The design file contains several columns entitled 'reference'. ")
    }
  } else {
    message("The headers of the design table do not contain 'reference' column -> taking 1st column as reference")
    reference.column <- 1
  }
  
  ## Identify test column in the design file
  if ("test" %in% (tolower(colnames(design)))) {
    test.column <- which(tolower(colnames(design)) == "test")
    if (length(test.column) != 1) {
      stop("The design file contains several columns entitled 'test'. ")
    }
  } else {
    message("The headers of the design table do not contain 'test' column -> taking 1st column as test")
    test.column <- 2
  }
  
  
  
  ## Build the result object
  result <- list()
  result$counts <- counts
  result$parameters <- parameters
  result$sample.desc <- sample.desc
  result$sample.ids <- sample.ids
  result$sample.labels <- sample.labels
  result$sample.conditions <- sample.conditions
  result$conditions <- conditions
  result$condition.table <- condition.table
  result$design <- design
  result$comparison.summary <- comparison.summary
  result$reference.column <- reference.column
  result$test.column <- test.column
  
  return(result)
}
