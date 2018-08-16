#' @title generate a figure chunk of R code to insert in an Rmd report
#' @author Jacques van Helden
#' @description generate a figure chunk of R code to insert in an Rmd report.
#'
#' @param name name of the chunk
#' @param figureFile path of the figure file
#' @param reportFile path of the report file. The link will be computed relative to this file
#' @param out.width="75\%" figure width, relative to page width
#' @param report.text=NULL if specified, the chunk will be appended to it
#' @param chunk_opt=",eval=TRUE" a string specifying one or more options to append to the chunk declaration line. Must start with a comma.
#'
#' @export
ReportFigure <- function(name, 
                          figureFile,
                          reportFile,
                          report.text = NULL,
                          out.width="75%",
                          chunk.opt = ", eval=TRUE") {
  fig.chunk <-  paste(sep = "",
                      "
```{r fig='", name, "', out.width = '", out.width, "'", chunk.opt, " }
knitr::include_graphics(path = '", RelativePath(source = reportFile, target = figureFile), "', auto_pdf = TRUE)
```")
  
  if (is.null(report.text)) {
    return(fig.chunk)
  } else {
    return(append(report.text, fig.chunk))
  }
}
