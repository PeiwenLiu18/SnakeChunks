#' @title scaling RNA-seq count table with a variety of methods (including edgeR and DESeq2 methods)
#'
#' @author Jacques van Helden and Mustafa AbuElqumsan
#'
#' @description normalisation of RNA-seq count table.
#' More precisely this function runs a sample-wise scaling so that all the samples
#' have the same value for a user-defined scaling parameter.
#' By default, we use the percentile 75 as scaling factor.
#'
#' @param counts a data frame with counts per feature, with one feature (gene, transcript) per row and one sample per column.
#' @param class.labels the class labels associated each sample of the count table. Should be provided in order to adapt it in case samples would be suppressed because they have a scaling factor of 0.
#'
#' @param nozero=TRUE If TRUE, zero values are ommitted from the data table before computing the column-wise statistics.
#'
#' @param method="percentile" normalization method. Optionally, can be specified as a vector with several methods.
#'
#' Supported normalization methods: none, sum, mean, median, percentile, quantiles, TMM, RLE, DESeq2.
#'
#' The "raw" methods simply returns the raw counts (no normalization).
#'
#' Sum and mean are really not recommended because very sensitive to outliers.
#' They are implemented only for the sake of comparison.
#'
#' Percentile: Preferred percentile = 75. Percentile-base scaling is more robust to
#' outliers than the mean or sum. Note that for RNA-seq count tables, sample median counts may
#' be unsatisfying  because for some datasets even the median counts are very low. A common
#' approach is to use the third quartile (equal to percentile 75), which seems to provide a
#' good tradeoff between robustness and representativity (taking into account a
#' representative proportion of the total counts per samples).
#'
#' Median: actually runs percentile-based scaling factor with percentile=50.
#'
#' TMM: trimmed mean of M-values proposed by Robinson and Oshlack (2010), computed via edgeR::calcNormFactors().
#'
#' RLE: relative log expression proposed by Anders and Huber (2010), computed via edgeR::calcNormFactors().
#'
#' DESeq2: compute size factors via DESeq2::estimateSizeFactors()
#'
#' @param percentile=75 reference percentile used as scaling factor when percentile method is selected.
#'
#' @param log2=FALSE  Apply log2 transformation after sample size correction
#'
#' @param epsilon=0.1 value added to all counts before applying the log2 transformation
#' in order to avoid zero counts.
#'
#' @param detailed.sample.stats=FALSE compute detailed sample stats (takes some seconds)
#'
#' @return a list with the following elements.
#' \itemize{
#'   \item{counts} The original count table
#'   \item{class.labels} The original class labels
#'   \item{parameters}  Normalisation parameters
#'   \item{norm.result} Normalisation results. If several methods are provided, normcounts is a list with one entry per normalisation method.
#'}
#' Each normalisation results contains.
#' \itemize{
#'   \item{size.factor} Size factor for each sample.
#'   \item{scaling.factor} Scaling factor for each sample. Scaling.factor equals 1/size.factor
#'   \item{counts} Normalized counts
#' }
#'
#' @export
NormalizeCountTable <- function(counts,
                                class.labels,
                                nozero = TRUE,
                                method = "TMM",
                                percentile = 75,
                                log2 = FALSE,
                                epsilon = 0.1,
                                detailed.sample.stats = FALSE,
                                verbose = 1) {

  ## Initialise return list
  result <- list()
  result$raw.counts <- counts
  result$class.labels <- class.labels


  ## if required, discard the zero counts before computing size factors
  if (is.null(nozero)) {
    ## By default, nozero parameter is activated
    nozero <- TRUE
    # View(head(non.null.counts))
  }
  result$nozero <- nozero

  ## Include parameters in the result
  result$parameters <- list(
    method = method,
    nozero = nozero,
    percentile = percentile,
    log2 = log2,
    epsilon = epsilon
  )

  #### Measure size of input matrix ####
  result$raw.features <- nrow(counts)
  result$raw.samples <- ncol(counts)

  if (verbose >= 1) {
    message("\t", "Library size normalization\tMethod: ", paste(collapse=", ", method),
            "\tRaw counts: ", result$raw.features, " features x ", result$raw.samples, " samples. ")
  }

  #### Compute non-zero counts ####
  counts.nozero <- counts
  counts.nozero[counts  == 0] <- NA
  if (nozero) {
    if (verbose >= 2) {
      message("\t\tIgnoring zero values for normalization")
    }
    # View(head(non.null.counts))
    counts.to.norm <- counts.nozero ## Use only non-zero counts for normalisation
  } else {
    counts.to.norm <- counts ## Use all counts for normalization
  }


  #### Compute sample-wise statistics ####
  sampleStats <- data.frame(
    zeros = apply(counts == 0, 2, sum),
    na.values = apply(is.na(counts), 2, sum),
    na.values.nozero = apply(is.na(counts.nozero), 2, sum),
    sum.all = apply(counts, 2, sum, na.rm = TRUE),
    sum.nozero = apply(counts.nozero, 2, sum, na.rm = TRUE),
    mean.all = apply(counts, 2, mean, na.rm = TRUE),
    mean.nozero = apply(counts.nozero, 2, mean, na.rm = TRUE),
    sd.all = apply(counts, 2, sd, na.rm = TRUE),
    sd.nozero = apply(counts.nozero, 2, sd, na.rm = TRUE)
  )

  # head(sampleStats)
  ## Report number of NA and zero values
  if (verbose >= 2) {
    message("\t\tOrginal data table contains ",
            format(x = sum(sampleStats$zeros), big.mark = ","),
            " zeros ")
    if (sum(sampleStats$na.values) > 0) {
      message("\t\tRaw counts contain ",
              format(x = sum(sampleStats$na.values), big.mark = ","),
              " NA values. ")
    }
  }

  if (detailed.sample.stats) {
    if (verbose >= 2) {
      message("\t\t", "Computing sample-wise statistics with and without zeros")
    }
    sampleStats$all <- RowStats(counts)
    sampleStats$nozero <- RowStats(counts.nozero)
  }
  result$sampleStats <- sampleStats
  # head(self$sampleStats)

  #### Apply normalization method(s) ####
  method.names <- vector()
  for (m in method) {
    #### Define method name ####
    if (m == "percentile") {
      if (!exists("percentile")) {
        stop("NormalizeSamples()\tMissing required parameter: standardization percentile")
      }
      if (is.null(percentile)) {
        stop("NormalizeSamples()\tpercentile-based scaling requires a non-null percentile parameter.")
      }
      if (percentile == 25) {
        method.name <- "Q1"
      } else if (percentile == 75) {
        method.name <- "Q3"
      } else if (percentile == 50) {
        method.name <- "median"
      } else {
        method.name <- paste(sep = "", "q", percentile)
      }
    } else {
      method.name <- m
    }
    if (log2) {
      ## Append log2 suffix if log2-transformation is activated
      method.name <- paste(sep = "", method.name, "_log2")
    }
    method.names <- append(method.names, method.name)

    if (verbose >= 2) {
      message("\t\tNormalization method\t", method.name)
    }

    if (m %in% c("raw", "none")) {
      scaling.factor <- rep(x = 1, length.out = ncol(counts))
      size.factor <-  rep(x = 1, length.out = ncol(counts))

    } else if (m %in% c("percentile", "median")) {
        ## Median will be treated as percentile 0.5
      if (m == "median") {
        current.percentile <- 50
      } else {
        current.percentile <- percentile
      }
      current.quantile <- current.percentile/100

      percentile.method <- "custom" ## alternative: compute percentile-based scaling factors via edgeR
      if (percentile.method == "edgeR") {

        ## Compute quantile-based scaling factor via edgeR
        ## NOTE (2018-07-21) : with single-cell data containing MANY zeros, this returns Inf scaling factors for almost all the samples
        if (verbose >= 3) {
          message("\t\tNormalizing counts with edgeR::calcNormFactors(method=upperquartile, p=", current.quantile,")")
        }
        d <- DGEList(counts = counts, group = class.labels)
        # d$samples$group <- relevel(d$samples$group)
        d <- calcNormFactors(d, method = "upperquartile", p = current.quantile)                 ## Compute normalizing factors
        scaling.factor <- d$samples$norm.factors
        size.factor <- 1/scaling.factor

      } else {
        if (verbose >= 3) {
          message("\t\tScaling factor: sample percentile ", current.percentile)
        }
        sampleStats$norm.percentile <- apply(counts.to.norm, 2, quantile, na.rm = TRUE, probs = current.quantile)
        size.factor <-  sampleStats$norm.percentile
        scaling.factor <- 1 / sampleStats$norm.percentile
        # mean(scaling.factor[!is.infinite(scaling.factor)])
        # hist(scaling.factor, breaks = 1000)
      }
      null.scaling <- sum(size.factor == 0)
      # inf.scaling <- sum(is.infinite(scaling.factor))
      if (null.scaling > 1) {
        message("\t\tdiscarding ", null.scaling, " samples with null value for percentile ", percentile)
      }

    } else if (m %in% c("mean", "libsum", "TC", "sum")) {
      if (verbose >= 3) {
        message("\t\tScaling factor: library size (equivalent for sum, mean, total counts).  ")
      }
      ## Note: mean is very sensitive to outliers, which are very problematic with RNAseq data
      size.factor <- apply(counts.to.norm, 2, mean, na.rm = TRUE)
      scaling.factor <- 1 / size.factor
      scaling.factor <- scaling.factor / mean(scaling.factor)
      # hist(scaling.factor, breaks = 1000)
      # mean(scaling.factor)


    } else if (m == "TMM") {
      if (verbose >= 3) {
        message("\t\tRunning edgeR::calcNormFactors(method='TMM').")
      }
      d <- DGEList(counts = counts, group = class.labels)
      # d$samples$group <- relevel(d$samples$group) ## Ensure that condition 2 is considered as the reference
      d <- calcNormFactors(d, method = "TMM")                 ## Compute normalizing factors
      scaling.factor <- d$samples$norm.factors
      size.factor <- 1/scaling.factor

    } else if (m == "RLE") {
      ## Run edgeR to compute the relative log expression defined by Anders and Huber (2010).
      if (verbose >= 3) {
        message("\t\tRunning edgeR::calcNormFactors(method='RLE').")
      }
      d <- DGEList(counts = counts, group = class.labels)
      # d$samples$group <- relevel(d$samples$group) ## Ensure that condition 2 is considered as the reference
      d <- calcNormFactors(d, method = "RLE")                 ## Compute normalizing factors
      scaling.factor <- d$samples$norm.factors
      size.factor <- 1/scaling.factor

    } else if (m == "DESeq2") {
      if (verbose >= 3) {
        message("\t\tRunning DESeq2::estimateSizeFactors()")
      }
      ## Replace non-alphanumeric charactersi by "_" in class labels, to avoid message from DESeq2
      clean.class.labels <- gsub(pattern = "[^[:alnum:] ]", replacement = "_", class.labels)
      clean.class.labels <- gsub(pattern = " ", replacement = "_", clean.class.labels)
      # unique(clean.class.labels)
      dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = data.frame(classes = clean.class.labels), ~ classes)
      dds <- estimateSizeFactors(dds)
      size.factor <- sizeFactors(dds)
      scaling.factor <- 1/size.factor
      # plot(size.factor, sampleStats$sum) ## THE DIFFERENCE IS QUITE IMPRESSIVE

    } else if (m == "VSD") {
      # Compute variance stabilizing transformations (VST) via DESeq2 (Tibshirani 1988; Huber et al. 2003; Anders and Huber 2010)
      if (verbose >= 3) {
        message("\t\tRunning DESeq2::vsd()")
      }
      stop("NOT FINISHED THE IMPLEMENTATION YET")
      ## Create a DESeqDataset object from the count table
      dds <- DESeqDataSetFromMatrix(counts = counts, colData = data.frame(classes = class.labels), ~ classes  )
      vsd <- vst(dds, blind = FALSE)
      names(vsd)
      ## QUESTION: HOW DO I GET THE SIZE FACTORS FROM THE RESULTING OBJECT ?

      #    rld <- rlog(dds, blind = FALSE)

      # View(vsd)
      # head(assay(vsd), 3)
      # library("pheatmap")
      # library("vsn")
      # meanSdPlot(assay(vsd))

      ## Run  differential expression analysis with DESeq2
      #    dds <- DESeq(dds)


      #    rld <- rlog(dds, blind=FALSE)

    } else {
      stop(method, " is not a valid method for NormalizeSamples()")
    }

    #### Discarded samples ####
    ## Detect problems related to null scaling factors, which may happen in some datasets due to a very large number of zeros.
    discardedSamples <-
      (size.factor == 0) |
      is.infinite(scaling.factor) |
      is.na(size.factor)
    discaredSampleNames <- vector()
    if (sum(discardedSamples) > 0) {
      discaredSampleNames <- colnames(counts.to.norm[, discardedSamples])
      message("\t\tDiscarding ", sum(discardedSamples), " samples because their size factor is null or NA. ")
      if (verbose >= 3) {
        message("\t\tDiscarded samples: ", paste(collapse = "; ", discaredSampleNames))
      }
    }

    #### Compute normalised counts ####
    normTarget <- mean(scaling.factor[!discardedSamples]) ## Ensure library eize equality before and after standardization
    scaling.factor <- scaling.factor * normTarget
    normCounts <- t(t(counts.to.norm[, !discardedSamples]) * scaling.factor)

    #### log2 transformation (if required) ####
    if (log2) {
      normCounts <- log2(normCounts + epsilon)
    }
    norm.result <- list(
      method = m,
      method.names = method.names,
      size.factor = size.factor,
      scaling.factor = scaling.factor,
      normCounts = normCounts)

    if (length(method) == 1) {
      result$norm.result <- norm.result
    } else {
      result[[method.name]] <- norm.result
    }
  }

  ## Get the method name(s) in the result
  result$method.name <- method.names

  return(result)
}
