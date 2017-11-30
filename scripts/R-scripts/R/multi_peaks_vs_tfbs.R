#' @title Compare several sets of genomic regions (e.g. Chip_seq peaks) with a set of reference genomic sites (e.g. transcription factor binding sites). 
#'
#' @author Claire Rioualen (\email{rioualenclaire@gmail.com}) Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Evaluate and compare the accuracy of several peaksets by comparing each of 
#' them with a set of reference TFBSs. The comparison betwee a peak set and the TFBSs is done 
#' with PeaksVsSTBS, and this function recollects the results and presents them in a synthetic way. 
#'  
#' @details
#' First version: 2017-09
#' @param peakFiles a vector contianing a list of files. Each file is supposed to contain a set of genomic regions, in bed format.   
#' @param sites path to a bed-formatted file containing the reference TFBSs. Passed to 
#' @param main="Mutual coverage plot" main title for the plot. 
#' @param peakSetLabels a vector of the same size as peakFiles indicating a label for each peak set. 
#' @param verbose=FALSE if true, print messages during the execution
#' @param drawPlots=FALSE if true, draw diagnostic plots
#' @param ... Additional parameters are passed to plot()
#' 
#' @return
#' 
#' A table with one row per peak set, and one column per evaluation parameter (those reported by PeaksVsTFBS). 
#' 
#' @examples
#' ## Load libraries
#' library(GenomicRanges)
#' 
#' ## TEMPORARY: these two functions were borrowed from the roken package, still to be published
#' source("~/SnakeChunks/scripts/R-scripts/R/bed_to_granges.R")
#' source("~/SnakeChunks/scripts/R-scripts/R/bed_import.R")
#' source("~/SnakeChunks/scripts/R-scripts/R/peaks_vs_tfbs.R")
#' 
#' ## Define the root directory for this analysis
#' mainDir <- "~/RegulonHT_results/GSE93506"
#' setwd(mainDir)
#' 
#' ## File containing the reference TF binding sites
#' TF <- "OmpR"
#' bedFile <- paste(sep="",TF, "_unique_sites.bed")
#' TFBSFile <- file.path(mainDir, "RegulonDB", bedFile)
#' 
#' ## Define the directory containing the peaks
#' peakDir <- file.path(mainDir, "results", "peaks", paste(sep="", TF, "_vs_WT"))
#' 
#' ## Identify all the peaks in this directory (recursive search in all subfolders)
#' peakFiles <- list.files(path=peakDir, pattern=c(".*sickle_bowtie2.*\\.bed$"), recursive=TRUE, full.names = TRUE)
#' message("Peak directory contains ", length(peakFiles), " bed files")
#' 
#' ## Compare peaks and TFBSs
#' multiPeakVsRegDB <- MultiPeaksVsTFBS(
#'   peakFiles, main=TF,
#'   peakSetLabels = basename(sub(peakFiles, pattern="\\.bed$", replacement = "")), 
#'   TFBSFile, verbose=TRUE, 
#'   drawPlots=TRUE)
#' 
#' @export
MultiPeaksVsTFBS <- function(peakFiles,
                             TFBSFile,
                             PNGFile,
                             main="Mutual coverage plot",
                             peakSetLabels=NULL,
                             verbose = FALSE,
                             drawPlots = FALSE,
                             ...) {
  if (verbose) {
    message("Comparing ", length(peakFiles), " peak files with one TFBS file.")
  }

  ## Check if peakset-specitic labels were specified
  if (!is.null(peakSetLabels)) {
    ## Check it there are the sanem number of labels as peaks sets
    if (length(peakFiles) != length(peakSetLabels)) {
      stop("MultiPeaksVsTFBS() error: peakFiles and peakSetLabels have different lengths.")
    }
  }
  
  
  ## Prepare result list
  result <- list()
  result$nbPeakFiles <- length(peakFiles)
  result$peakFiles <- peakFiles
  result$TFBSFile <- TFBSFile
  result$peaksVsTFBStable <- data.frame()
  
  
  #peakFile <- peakFiles[1]
  for (i in 1:length(peakFiles)) {
    peakFile <- peakFiles[i]
    if (verbose) {
      message("MultiPeaksVsTFBS, peak file: ", peakFile)
    }
    
    ## Run the comparison between one peakset and the reference TFBS
    onePeakSetCompa <- PeaksVsTFBS(peakFile, TFBSFile, verbose, drawPlots)
    
    ## Collect the statistics
    result$peaksVsTFBStable <- rbind (
      result$peaksVsTFBStable,
      data.frame(
        peakNb = onePeakSetCompa$peakNb,
        coveredPeaks = onePeakSetCompa$coveredPeaks,
        missedPeaks = onePeakSetCompa$missedPeaks,
        peakCoverage = onePeakSetCompa$peakCoverage,
        siteNb = onePeakSetCompa$siteNb,
        coveredSites = onePeakSetCompa$coveredSites,
        missedSites = onePeakSetCompa$missedSites,
        siteCoverage = onePeakSetCompa$siteCoverage,
        peakFile = onePeakSetCompa$peakFile
      )
     )
    if (is.null(peakSetLabels)) {
      result$peakSetLabel <- 1:length(peakFiles)
    } else {
      result$peakSetLabel <- result$peakSetLabels
    }
  }
  
  ## Draw a summary plot
  if (drawPlots) {
    
    ## Coverage plot
    x <- 100*result$peaksVsTFBStable$siteCoverage
    y <- 100*result$peaksVsTFBStable$peakCoverage
    png(PNGFile)
    plot(x, y,
         main = main,
         xlim = c(0,100), 
         ylim = c(0, max(y) * 2), 
         xlab = "% TFBS",
         ylab = "% peaks",
         las=1,
         panel.first=grid(lty="solid", col="gray")
    )
    abline(v=seq(from=0, to=100, by=10), col="gray")
    text(x, y, labels = 1:length(peakFiles), pos = 4)
    if (!is.null(peakSetLabels)) {
      legend("topleft", 
             legend = paste(1:length(peakFiles), peakSetLabels), 
             cex=0.7)
    }
    dev.off()
  }
  
  return(result)
}
