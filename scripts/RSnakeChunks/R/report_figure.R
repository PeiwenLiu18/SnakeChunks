#' @title generate a figure chunk of R code to insert in an Rmd report
#' @author Jacques van Helden
#' @description generate a figure chunk of R code to insert in an Rmd report.
#'
#' @param name name of the chunk
#' @param figureFile path of the figure file
#' @param reportFile path of the report file. The link will be computed relative to this file
#' @param report.text=NULL if specified, the chunk will be appended to it
#' @param out.width="75\%" figure width, relative to page width
#' @param pdf.link=TRUE if TRUE, the chunk is inserted between square brackates followed with a link to the corresponding pdf figure (obtained by substituting the extension)
#' @param chunk_opt=",eval=TRUE" a string specifying one or more options to append to the chunk declaration line. Must start with a comma.
#'
#' @export
ReportFigure <- function(name,
                         figureFile,
                         reportFile,
                         report.text = NULL,
                         out.width="75%",
                         pdf.link = TRUE,
                         chunk.opt = ", eval=TRUE") {
  link.prefix <- "\n\n"
  link.suffix <- "\n\n"

  if (pdf.link) {
    figureFile.pdf <- sub(pattern = ".png", replacement = ".pdf", x = figureFile)
    figureFile.pdf <- sub(pattern = ".jpg", replacement = ".pdf", x = figureFile.pdf)
    figureFile.pdf <- sub(pattern = ".jpeg", replacement = ".pdf", x = figureFile.pdf)
    link.prefix <- paste(sep = "", link.prefix, "[")
    link.suffix <- paste(sep = "", "\n]", "(", RelativePath(source = reportFile, target = figureFile.pdf), ")", link.suffix)
  }

  fig.chunk <-  paste(
    sep = "",
    link.prefix,
    "
```{r fig='", name, "', out.width = '", out.width, "'", chunk.opt, " }
knitr::include_graphics(path = '", RelativePath(source = reportFile, target = figureFile), "', auto_pdf = TRUE)
```",
    link.suffix)

  if (is.null(report.text)) {
    return(fig.chunk)
  } else {
    return(append(report.text, fig.chunk))
  }
}
