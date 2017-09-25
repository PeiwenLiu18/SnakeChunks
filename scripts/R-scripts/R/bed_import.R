#' @title Import bed files
#' @author Zacharie Menetrier
#' @description Imports bed files automatically in a data frame.
#' 
#' @param path The path of the file to import.
#' 
#' @return A data frame with the corresponding chromosomic regions 
#' of the bed file.
#' 
#' @usage BedImport(path)
#' 
#' @examples 
#' chromFile <- system.file("extdata", "hg19.genome", package = "roken")
#' chromSizes <- BedImport(chromFile)
#' 
#' @export
BedImport <- function(path) {
  # # Get the data frame from the file path.
  # regions <- as.data.frame(data.table::fread(path, 
  #                                            header = FALSE,
  #                                            comment.char="#",
  #                                            sep = "\t", 
  #                                            stringsAsFactors = FALSE, 
  #                                            quote = ""))

  # Get the data frame from the file path.
  regions <- tryCatch(read.table(path, 
                        header = FALSE,
                        comment.char="#",
                        sep = "\t", 
                        quote = ""),  error=function(e) NULL)

  # Renaming the data frame columns to fall in the bed standard.
  columnNames <- c("chrom", "chromStart", "chromEnd", "name", "score", 
                   "strand", "thickStart",
                   "thickEnd", "itemRgb", "blockCount", "blockSizes",
                   "blockStarts")
  
  
  if (is.null(regions)) {
    ## If bed file does not contain a single line, return a data.frame with the right number of columns but 0 lines
    message("Warning: no line in bed file ", path)
    regions <- as.data.frame(t(data.frame(row.names = columnNames)))
  } else {
    colnames(regions) <- columnNames[1:ncol(regions)]
  }
  
  return(regions)
}
