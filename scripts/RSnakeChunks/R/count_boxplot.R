#' @title Draw boxplots with read counts per genes for each sample
#' @author Jacques van Helden
#' @param counts a count table with one row per feature and one column per sample
#' @param sample.desc sample description table, one row per sample
#' @param sample.label.col=NULL column containing sample-associated colors
#' @param xlab="Raw counts" X axis label
#' @param main="Box plots per sample: raw counts" Main title to display on the plot
#' @param plot.file=NULL if specified, the graph will be stored in a pdf-formatted file.
count.boxplot <- function(counts,
                          sample.desc,
                          sample.label.col = NULL,
                          xlab="Raw counts",
                          main="Box plots per sample: raw counts",
                          plot.file=NULL) {

  ## Adapt boxplot size to the number of samples and label sizes
  boxplot.lmargin <- max(nchar(sample.desc$label))/3+5
  boxplot.height <- length(sample.ids)/3+2

  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file = plot.file, width = 8, height = boxplot.height)
  }

  par(mar=c(5,boxplot.lmargin,4,1)) ## adapt axes
  boxplot(counts, horizontal= TRUE, col=sample.desc$color,
          xlab=xlab, names=sample.desc[, sample.label.col],
          main=main, las=1)
  grid(col="grey", lty="solid",ny = 0)
  if (!is.null(plot.file)) {
    silence <- dev.off; rm(silence)
  }
}
