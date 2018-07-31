#' @title Export the current device to one or several files with user-specified format(s).
#' @author Jacques van Helden
#' formats.
#' @examples
#' ExportPLot(file.prefix = "test",
#'    export.formats = c("postscript","jpg","png","bmp","pdf"))
#' @param  file.prefix="PlotExport"
#' @param export.formats="pdf" # supported: postscript, jpg, png, bmp, pdf
#' @param width=11 Width of output file (in inches)
#' @param height=8 Height of output file (in inches)
#' @param horizontal=TRUE parameter passed to postscript()
#' @param dpi=72 dots per inch. Automatically converted to width and height in pixels for bitmap formats. 
#' @param verbose = 0 level of verbosity
#' @param ... Additional parameters are passed to the export method
#' @export             
ExportPLot <- function(file.prefix = "PlotExport",
                        export.formats = "pdf", # supported: postscript, jpg, png, bmp, pdf
                        width = dev.size(units = "in")[1], # in inches
                        height = dev.size(units = "in")[1], # in inches
                        horizontal = TRUE,
                        dpi = 72,
                        verbose = 0, 
                        ... ## Additional parameters are passed to the export method
) {
  
  f <- export.formats[1]
  for (f in export.formats) {
    source.dev <- dev.cur();
    new.dev <- OpenPlotDevice(file.prefix, fig.format = f, width = width,height = height,horizontal = horizontal, dpi = dpi, verbose, ...)
#    new.dev <- dev.cur()
    dev.set(which = source.dev)
    dev.copy(which = new.dev)
    dev.set(which = new.dev)
    dev.off()
    dev.set(which = source.dev) ## This is required because dev.off() returns to the first, not the last, device
  }
}

#' @title open plot device with parameters adapted to the specificities of the different device functions. 
#' @author Jacques van Helden
#' @param file.prefix file prefix (the extension will be concatenated according to format)
#' @param fig.format="pdf" file firmat. Supported: postscript, jpg, png, bmp, pdf
#' @param width=8 Width of output file (in inches).
#' @param height=8 Height of output file (in inches)
#' @param horizontal=TRUE parameter passed to postscript()
#' @param dpi=72 dots per inch. Automatically converted to width and height in pixels for bitmap formats. 
#' @param verbose = 0 level of verbosity
#' @return returns the new device
#' @export
OpenPlotDevice <- function(file.prefix,
                             fig.format = 'pdf',
                             width = 8,
                             height = 8,
                             horizontal = TRUE,
                             dpi = 72, ## Points per Inch (screen resolution by default)
                             verbose = 0,
                             ...) {  
  
  fig.format <- tolower(fig.format)  
  file.ext <- c(
    x11 = "x11",
    postscript = "ps",
    pdf = "pdf",
    ps = "ps",
    eps = "eps",
    jpeg = "jpg",
    jpg = "jpg",
    bmp = "bmp",
    png = "png")
  file.name <- paste(file.prefix, file.ext[fig.format], sep = ".")
  if (verbose >= 0) { message("\t\tExporting plot to file ",file.name, " ", fig.format, " format", sep = "") }
  if ((fig.format   ==   "postscript") || (fig.format   ==   "ps")) {
    postscript(file.name,paper = "special",width = width,height = height,horizontal = horizontal, ...)
  } else if (fig.format   ==   "eps") {
    postscript(file.name,paper = "special",width = width,height = height,horizontal = horizontal,onefile = F, ...)
  } else if (fig.format   ==   "pdf") {
    pdf(file.name, paper = "special",width = width,height = height, ...)
  } else if ((fig.format   ==   "jpg") || (fig.format  ==  "jpeg")) {
    jpeg(file.name,width = width * dpi, height = height * dpi, quality = 100, ...)
  } else if (fig.format  ==  "png") {
    png(file.name,width = width * dpi, height = height * dpi, ...)
  } else if (fig.format  ==  "bmp") {
    bitmap(file.name,width = width*dpi,height = height*dpi, ...)
  } else if (fig.format  ==  "x11") {
    x11(width = width,height = height)
  } else {
    stop("OpenPlotDevice()\tInvalid format: ", fig.format, ". Supported: ps, eps, pdf, jpg, png, bmp, x11")
    return()
  }
  new.dev <- dev.cur()
  return(new.dev)
}
