#' @title generate a figure chunk of R code to insert in an Rmd report
#' @author Jacques van Helden
#' @param name name of the chunk
#' @param file file continaing the figure
#' @param out.width="75%" figure width, relative to page width
#' @param index.text=NULL if spcified, the chunk will be appended to it
#' @param chunk_opt=", eval=TRUE" a string specifying one or more options to append to the chunk declaration line. Must start with a comma.
#' @export
index.figure <- function(name, file,
                         index.text = NULL,
                         out.width="75%",
                         chunk.opt = ", eval=TRUE") {
  fig.chunk <-  paste(sep = "",
                      "
```{r fig='", name, "', out.width = '", out.width, "'", chunk.opt, " }
knitr::include_graphics(path = '", file, "', auto_pdf = TRUE)
```")

  if (is.null(index.text)) {
    return(fig.chunk)
  } else {
    return(append(index.text, fig.chunk))
  }
}
