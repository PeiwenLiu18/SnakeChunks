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
#' ## Define the directory containing the peaks
#' peakDir <- file.path(mainDir, "results", "peaks", "Nac_vs_WT", "macs2")
#' 
#' ## Identify all the peaks in this directory (recursive search in all subfolders)
#' peakFiles <- list.files(path=peakDir, pattern="*macs2.bed$", recursive=TRUE, full.names = TRUE)
#' message("Peak directory contains ", length(peakFiles), " bed files")
#' 
#' ## File containing the reference TF binding sites
#' TFBSFile <- file.path(mainDir, "RegulonDB", "Nac_unique_sites.bed")
#' 
#' ## Compare peaks and TFBSs
#' multiPeakVsRegDB <- MultiPeaksVsTFBS(peakFiles, TFBSFile, drawPlots=TRUE)
#' 
#' @export
MultiPeaksVsTFBS <- function(peakFiles,
                             TFBSFile,
                             verbose = FALSE,
                             drawPlots = FALSE,
                             ...) {
  if (verbose) {
    message("Comparing ", length(peakFiles), " peak files with one TFBS file.")
  }

  ## Prepare result list
  result <- list()
  result$nbPeakFiles <- length(peakFiles)
  result$peakFiles <- peakFiles
  result$TFBSFile <- TFBSFile
  result$peaksVsTFBStable <- data.frame()
  
  for (peakFile in peakFiles) {
    ## Run the comparison between one peakset and the reference TFBS
    onePeakSetCompa <- PeaksVsTFBS(peakFile, TFBSFile, verbose, drawPlots)
    
    ## Collect the statistics
    result$peaksVsTFBStable <- rbind (
      result$peaksVsTFBStable,
      data.frame(
        peakNb=onePeakSetCompa$peakNb,
        coveredPeaks=onePeakSetCompa$coveredPeaks,
        missedPeaks=onePeakSetCompa$missedPeaks,
        peakCoverage=onePeakSetCompa$peakCoverage,
        siteNb=onePeakSetCompa$siteNb,
        coveredSites=onePeakSetCompa$coveredSites,
        missedSites=onePeakSetCompa$missedSites,
        siteCoverage=onePeakSetCompa$siteCoverage
      )
    )
  }
  
  ## Draw a summary plot
  if (drawPlots) {
    plot(result$peaksVsTFBStable$siteCoverage,
         result$peaksVsTFBStable$peakCoverage,
         xlim=c(0,1), xlab = "TFBS coverage by peaks",
         ylim=c(0,1), ylab = "peak coverage by TFBSs"
    )
  }
  
  return(result)
}