#' @title Draw boxplots with read counts per genes for each sample
#' @author Jacques van Helden
#' @description draw a box plot with one lane per sample (column) of an expression table. 
#' @param counts a count table with one row per feature and one column per sample
#' @param sample.desc sample description table, one row per sample. If it contains a column named "color", this colum will be used to define sample-specific colors. 
#' @param sample.label.col=NULL column containing sample-associated colors
#' @param xlab="Raw counts" X axis label
#' @param main="Box plots per sample: raw counts" Main title to display on the plot
#' @param plot.file=NULL if specified, the graph will be stored in a pdf-formatted file.
#' @export
BoxplotsPerSample <- function(counts,
                          sample.desc,
                          sample.label.col = NULL,
                          xlab="Raw counts",
                          main="Box plots per sample: raw counts",
                          plot.file=NULL) {
  
  ## Sample colors
  if(is.null(sample.desc$color)) {
    sample.colors <- rep(x = c("#BBBBBB", "#DDDDDD"), length.out = nrow(sample.desc))
  } else {
    sample.colors <- as.vector(unlist(sample.desc$color))
  }

  ## Sample labels
  if (is.null(sample.label.col)) {
    sample.labels <- rownames(sample.desc)
  } else {
    sample.labels <- as.vector(unlist(sample.desc[, sample.label.col]))
  }
  
  ## Adapt boxplot size to the number of samples and label sizes
  boxplot.lmargin <- ceiling(max(nchar(sample.labels))/3+5)
  boxplot.height <- ceiling(length(sample.ids)/3+2)

  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file = plot.file, width = 8, height = boxplot.height)
  }

  par(mar=c(5, boxplot.lmargin, 4, 1)) ## adapt axes
  boxplot(counts, horizontal= TRUE, col=sample.colors,
          xlab=xlab, names=sample.desc[, sample.label.col],
          main=main, las=1)
  grid(col="grey", lty="solid",ny = 0)
  if (!is.null(plot.file)) {
    silence <- dev.off; rm(silence)
  }
}
