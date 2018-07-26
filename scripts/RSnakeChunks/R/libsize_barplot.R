#' @title Draw a barplot with the number of reads per sample
#' @author Jacques van Helden
#' @param counts a count table with 1 row per feature and 1 column per sample
#' @param sample.labels=colnames(counts) sample labels (displayed as y legend)
#' @param sample.colors=NULL sample-specific colors (e.g. reflecting conditions)
#' @param main Main title for the plot
#' @param plot.file=NULL Path of a file to store the figure (pdf format).
LibsizeBarplot <- function(counts,
                           sample.labels=colnames(counts),
                           sample.colors=NULL,
                           main = "Library sizes",
                           plot.file = NULL,
                           ...) {

  ## Store original graphical parameters
  par.ori <- par(no.readonly = TRUE)

  ## Compute library sizes
  libsizes <- round(x = apply(counts, 2, sum, na.rm = TRUE) / 1e6, digits = 2)

  ## Sample colors
  if (is.null(sample.colors)) {
    #sample.colors <- colorRampPalette(c('blue','red'))(100)[as.numeric(cut(libsizes, breaks = 100))]
    sample.colors <- rep(c("#888888", "#DDDDDD"), length.out = length(libsizes))
  }

  ## Adapt boxplot size to the number of samples and label sizes
  boxplot.lmargin <- max(nchar(sample.labels))/3 + 5
  boxplot.height <- length(sample.ids)/3 + 2

  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file = plot.file, width = 8, height = boxplot.height)
  }

  par(mar = c(5, boxplot.lmargin, 4, 1)) ## adapt axes
  bplt <- barplot(libsizes,
                  names.arg = sample.labels,
                  main = main,
                  horiz = TRUE, las = 1,
                  xlab = "libsum (Million reads per sample)",
                  col = sample.colors, ...)
  grid(col = "white", lty = "solid", ny = 0)
  text(x = pmax(libsizes, 3),
       labels = libsizes, y = bplt, pos = 2, font = 2)
  if (!is.null(plot.file)) {
    silence <- dev.off(); rm(silence)
  }
  par(par.ori)
#  return(bplt)
}



#' #' @title Draw a barplot with the number of reads per sample
#' #' @author Jacques van Helden
#' #' @param stats.per.sample a table with the sample-wise statistics, which ca n be produced by RowStats().
#' #' @param main Main title for the plot
#' #' @plot.file=NULL Path of a file to store the figure (pdf format).
#' libsize.barplot <- function(stats.per.sample,
#'                             main = "Read library sizes (libsum per sample)",
#'                             plot.file = NULL,
#'                             ...) {
#'
#'   ## Store original graphical parameters
#'   par.ori <- par(no.readonly = TRUE)
#'
#'   ## Get sample IDs
#'   if (is.null(stats.per.sample$ID)) {
#'     sample.ids <- as.vector(rownames(stats.per.sample))
#'   } else {
#'     sample.ids <- as.vector(stats.per.sample$ID)
#'   }
#'
#'
#'   ## Check if sample descriptions contain label column. If not, use ID.
#'   if (is.null(stats.per.sample$label)) {
#'     sample.labels <- sample.ids
#'   } else {
#'     sample.labels <- stats.per.sample$label
#'   }
#'
#'   ## Adapt boxplot size to the number of samples and label sizes
#'   boxplot.lmargin <- max(nchar(sample.labels))/3 + 5
#'   boxplot.height <- length(sample.ids)/3 + 2
#'
#'   ## Sample-wise library sizes
#'   if (!is.null(plot.file)) {
#'     message("Generating plot", plot.file)
#'     pdf(file = plot.file, width = 8, height = boxplot.height)
#'   }
#'
#'   par(mar = c(5,boxplot.lmargin,4,1)) ## adapt axes
#'   bplt <- barplot(stats.per.sample$Mreads,
#'                   names.arg = sample.labels,
#'                   main = main,
#'                   horiz = TRUE, las = 1,
#'                   xlab = "libsum (Million reads per sample)",
#'                   col = stats.per.sample$color, ...)
#'   grid(col = "white", lty = "solid", ny = 0)
#'   text(x = pmax(stats.per.sample$Mreads, 3),
#'        labels = stats.per.sample$Mreads, y = bplt, pos = 2, font = 2)
#'   if (!is.null(plot.file)) {
#'     silence <- dev.off(); rm(silence)
#'   }
#'   par(par.ori)
#' }
