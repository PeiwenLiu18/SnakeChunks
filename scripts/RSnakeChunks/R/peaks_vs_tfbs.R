#' @title Compare peaks with transcription factor binding sites (TFBSs)
#'
#' @author Claire Rioualen (\email{rioualenclaire@gmail.com}) Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Evaluate a set of ChIP-seq peaks (or any analogous type of genomic regions) by comparing it with some reference collection of transcription factor binding sites (TFBSs). 
#' 
#' @details
#' First version: 2017-09
#' @param peaks path to a bed-formatted file containing the peaks
#' @param sites path to a bed-formatted file containing the reference TFBSs
#' @param verbose=FALSE if true, print messages during the execution
#' @param drawPlots=FALSE if true, draw diagnostic plots
#' @param ... Additional parameters are passed to plot()
#' 
#' @return
#' A list with fields describing the different elements of comparison {TO BE DETAILED}
#' 
#' @examples
#' ## Load libraries
#' library(GenomicRanges)
#' 
#' ## TEMPORARY: these two functions were borrowed from the roken package, still to be published
#' source("~/SnakeChunks/scripts/R-scripts/R/bed_to_granges.R")
#' source("~/SnakeChunks/scripts/R-scripts/R/bed_import.R")
#' 
#' ## Define the root directory for this analysis
#' mainDir <- "~/RegulonHT_results/GSE93506"
#' setwd(mainDir)
#' 
#' peakFile <- file.path(mainDir, "results", "peaks", "Nac_vs_WT", "macs2", "Nac_vs_WT_sickle_bowtie2_macs2.bed")
#' TFBSFile <- file.path(mainDir, "RegulonDB", "Nac_unique_sites.bed")
#' 
#' ## Compare peaks and TFBSs
#' peakVsRegDB <- PeaksVsTFBS(peakFile, TFBSFile, drawPlots=TRUE)
#' 
#' @export


## !!!!!!
## JvH:  Il ne faut pas passer ces paramètres dans la fonction R elle-même, sinon elle ne peut être utilisée qu'à partir de snakemake. 
## Ces paramètres pourraient sans doute être spécifiés dans le bloc de code R de la règle. 
## A discuter; pour le moment je les commente pour pouvoir travailler hors snakemake. 
##
# peakFile <- snakemake@input[["peaks"]]
# TFBSFile <- snakemake@input[["tfbs"]]
##
## !!!!!!

PeaksVsTFBS <- function(peakFile,
                        TFBSFile,
                        verbose = FALSE,
                        drawPlots = FALSE,
                        ...) {
  require(GenomicRanges)
  
  ## Read peaks and sites from the bed files, and convert them to GenomicRanges
  if (verbose) {
    message("Loading peaks from bed file ", peakFile)
  }
  peaks <- BedToGranges(peakFile)
  peakNb <- length(attr(peaks, "seqnames"))

  if (verbose) {
    message("Loading TFBS from bed file ", TFBSFile)
  }
  sites <- BedToGranges(TFBSFile)
  siteNb <- length(attr(sites, "seqnames"))
  
  if (verbose) {
    message("Loaded ", peakNb, " peaks and ", siteNb, " TFBSs") 
  }
  
  
  ## We don't use this later, I comment it
  # overlaps <- GenomicRanges::findOverlaps(peaks, sites, type = "any")
  # intersections <- GenomicRanges::intersect(peaks, sites)
  # class(intersections) ## Check the classof the resulting object
  # attributes(intersections) ## Get the names of attibutes for this lass of objects
  # intersectNames <- intersections@seqnames ## Get the attribute "seqname" of the intersections object
  # class(intersectNames)

  ## Count the number of sites per peak
  sitesPerPeak <- GenomicRanges::countOverlaps(peaks, sites, type = "any")
  peaksPerSite <- GenomicRanges::countOverlaps(sites, peaks, type = "any")
  
  ## Peak statistics
  result <- list()
  result$peakNb <- peakNb
  result$sitesPerPeak <- sitesPerPeak ## Distribution of the number of sites per peak
  result$sitesPerPeakDistrib <- as.data.frame.table(table(sitesPerPeak))
  result$coveredPeaks <- sum(sitesPerPeak > 0)  ## Number of sites overlapped by at least one peak
  result$missedPeaks <- sum(sitesPerPeak == 0)  ## Number of sites overlapped by at least one peak
  if (peakNb == 0) {
    result$peakCoverage <- 0
  } else {
    result$peakCoverage <- result$coveredPeaks / peakNb
  }
  
  ## Site statistics
  result$siteNb <- siteNb
  result$peaksPerSite <- peaksPerSite
  result$peaksPerSiteDistrib <- as.data.frame.table(table(peaksPerSite))
  result$coveredSites <- sum(peaksPerSite > 0)  ## Number of sites overlapped by at least one peak
  result$missedSites <- sum(peaksPerSite == 0)  ## Number of sites overlapped by at least one peak
  if (siteNb == 0) {
    result$siteCoverage <- 0
  } else {
    result$siteCoverage <- result$coveredSites / siteNb
  }
  
  result$peakFile <- c(toString(peakFile))
  ################################################
  ## Draw some diagnostic plots if requested
  
  ## Define clor palette
  myColors <- c("peaks"="#88BBFF", 
                "sites"="#88FFBB")
  
  if (drawPlots) {
    par.ori <- par(no.readonly = TRUE) ## Store original plotting parameters
    par(mfrow=c(2,1))
    par(mar=c(5.1, 4.1, 4.1, 1.1))
    par(las=1)
    
    ## Peaks per site
    if (siteNb == 0) {
      plot(0, 0, main = "Not a single site", type="n", 
           xlab="Number of peaks", 
           ylab="Number of sites")
    } else {
      png(paste(peakFile, "_peaks_per_site.png", sep=""))
      hist(peaksPerSite, breaks = seq(from=min(peaksPerSite), to=max(peaksPerSite)+1)-0.5, 
           col=myColors["sites"], main="Peaks per sites", 
           xlab="Number of peaks", 
           ylab="Number of sites", ...)
      dev.off()
    }

    ## Sites per peak
    if (peakNb == 0) {
      plot(0, 0, main = "Not a single peak", type="n", 
           xlab="Number of sites", 
           ylab="Number of peaks")
    } else {
      png(paste(peakFile, "_sites_per_peak.png", sep=""))
      hist(sitesPerPeak, breaks = seq(from=min(sitesPerPeak), to=max(sitesPerPeak)+1)-0.5, 
         col=myColors["peaks"], main="Sites per peak", 
         xlab="Number of sites", 
         ylab="Number of peaks", ...)
      dev.off()
    }
    par(mfrow=c(1,1))
    par <- par.ori ## Restore original parameters
  }
    
  return(result)
  
  
}
