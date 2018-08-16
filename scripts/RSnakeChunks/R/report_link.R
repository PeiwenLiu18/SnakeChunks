#' @title IncludeLink
#' @author Jacques van Helden
#' @param source path of the report file, where the link will be inserted
#' @param target path of the file to link
#' @param displayPath=TRUE in the link text, display either the path or the basename of the target
#' @param report.text=NULL if not null, append the link to it
#' @examples 
#' 
#' ## Display the (relative) path  in the link
#' ReportLink(
#'   source = "RNA-seq/results/diffexpr/reports/myreport.Rmd", 
#'   target = "RNA-seq/results/diffexpr/a_vs_b/DEG/deg_table.tsv",
#'   displayPath = TRUE)
#' 
#' ## Display only the file name, not its path
#' ReportLink(
#'   source = "RNA-seq/results/diffexpr/reports/myreport.Rmd", 
#'   target = "RNA-seq/results/diffexpr/a_vs_b/DEG/deg_table.tsv",
#'   displayPath = FALSE)
#' 
#' @export
ReportLink <- function(source, target, displayPath=TRUE, report.text=NULL) {
  if (displayPath) {
    targetText <- target
  } else {
    targetText <- basename(target)
  }
  link.text <- paste(sep = "", 
                     "[", targetText, "]",
                     "(", RelativePath(source = source, target = target, sourceIsDir = FALSE), ")")
  if (is.null(report.text)) {
    return(link.text)
  } else {
    return(append(report.text, link.text))
  }
}
