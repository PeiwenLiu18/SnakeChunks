#' @title compute the path of one file relative to another one
#' @author Jacques van Helden
#' @param source source, i.e. the file or directory from which the relative path starts (e.g. an index file)
#' @param target target, i.e. the file that will be pointed by the relative path
#' @param sourceIsDir=FALSE Boolean indicating whether the source is a file or a directory
#' @examples
#' 
#' ## From file to file with shared path
#' RelativePath(
#'   source = "RNA-seq/results/diffexpr/reports/myreport.Rmd", 
#'   target = "RNA-seq/results/diffexpr/a_vs_b/DEG/deg_table.tsv",
#'   sourceIsDir = FALSE)
#'
#' ## From file to file without shared path
#' RelativePath(
#'   source = "RNA-seq/results/diffexpr/reports/myreport.Rmd", 
#'   target = "metadata/sample_descriptions.tsv",
#'   sourceIsDir = FALSE)
#'
#' ## From file to file without source in current dir
#' RelativePath(
#'   source = "myreport.Rmd", 
#'   target = "metadata/sample_descriptions.tsv",
#'   sourceIsDir = FALSE)
#'
#' ## From file to dir
#' RelativePath(
#'   source = "RNA-seq/results/diffexpr/reports/myreport.Rmd", 
#'   target = "RNA-seq/results/samples/",
#'   sourceIsDir = FALSE)
#'   
#' ## From dir to dir: requires to specify sourceIsDIR=TRUE
#' RelativePath(
#'   source = "RNA-seq/results/diffexpr/reports", 
#'   target = "RNA-seq/results/samples",
#'   sourceIsDir = TRUE)
#'   
#' ## From file to upper dir
#' RelativePath(
#'   source = "RNA-seq/results/diffexpr/reports/myreport.Rmd", 
#'   target = "RNA-seq/results/",
#'   sourceIsDir = FALSE)
#'   
#' ## From file to its own dir
#' RelativePath(
#'   source = "RNA-seq/results/diffexpr/reports/myreport.Rmd", 
#'   target = "RNA-seq/results/diffexpr/reports/",
#'   sourceIsDir = FALSE)
#'   
#' @export
RelativePath <- function(source, target, sourceIsDir = FALSE, verbose = 0) {
  if (sourceIsDir) {
    sourcePath <- unlist(strsplit(x = source, split = "/"))
  } else {
    sourcePath <- unlist(strsplit(x = dirname(source), split = "/"))
  }
  if ((length(sourcePath) == 1) && (sourcePath == ".")) {
    sourcePath <- vector() ## replace dot by empty vector
  }
  targetPath <- unlist(strsplit(x = target, split = "/"))
  
  ## Compute match bewteen source and target paths
  pathMatch <- match(sourcePath, targetPath)
  if ((length(pathMatch) == 0) || is.na(pathMatch)) {
    common.depth <- 0
  } else {
    common.depth <- max(pathMatch, na.rm = TRUE)
  }  
  backPath <- paste(collapse = "/", rep(x = "..", length.out = length(sourcePath) - common.depth))
  
  
  if (length(targetPath) > common.depth) {
    fwdPath <- paste(collapse = "/", targetPath[(common.depth + 1):length(targetPath)])
    if (length(backPath) == 0) {
      relPath <- fwdPath
    } else {
      relPath <- file.path(backPath, fwdPath)
    }
  } else {
    relPath <- backPath
  }
  
  if (relPath == "") {
    relPath = "."
  }
  
  if (verbose >= 1) {
    message("\tRelativePath()")
    message("\t\tsource\t", source)
    message("\t\ttarget\t", target)
    message("\t\trelPath\t", relPath)
  }
  
  return(relPath)
}
