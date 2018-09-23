#' @title Draw an MA plot
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Draw an MA plot from a table containing at least one column with the log2 fold change (displayed as ordinate) and one column with the geometric mean (displayed in abcsissa).
#' @param result.table A data frame containing one row per feature, and one column per statistics.
#' @param log2FC.col="log2FC" A column number or name, inidicating which column of the result table contains the log2 fold-change (log2FC)
#' @param log2mean.col="log2mean"  A column number or name, indicating which column of the result table contains the log2 of the mean signal per feature.
#' @param control.type="padj"  A column number or name, indicating which column of the result table contains the p-value or an equivalent indication of the significance of each feature (example: "padj", "p.value")
#' @param nb.tests=nrow(multitest.table) Total number of tests performed (by default, the number of rows in the table).
#' @param alpha=0.05    Alpha threshold for the control of false positives
#' @param log2FC.threshold=NULL Threshold on the log2FC
#' @param sort.by.pval=FALSE Sort row by p-value in order to plot significant elements on top of non-significant
#' @param xlab=paste(sep="","Effect_size_(",log2FC.col,")") Label for the X axis.
#' @param ylab=paste(sep="","-log10(",control.type,")") Label for the Y axis
#' @param xlim    Range of the X axis.
#' @param ylim    Range of the Y axis.
#' @param density.colors=FALSE Automatically set the colors according to feature status and local density in the Volcano space.
#' @param col.points='#888888' Color(s) for the points. can be either a single value (same color for all points), or a vector of the same length as the numer of genes (rows) of the input table.
#' @param col.positive='#4455DD' Color to highlight significant points (called positive). When NULL, positive points are not displayed.
#' @param col.lines='blue'  Color for the line highlighting the effect size thresholds.
#' @param col.alpha="darkred" Color for the line highlighting the significance threshold.
#' @param col.grid='#AAAAAA'   Grid color
#' @param lty.lines="solid" Line type for the lines.
#' @param lty.alpha="dashed" Line type for the horizontal line denoting the alpha theshold.
#' @param lty.grid="dotted" Line type for the grid.
#' @param pch.positive=3   Point shape for the positive tests (default: 3, i.e. + symbol)
#' @param pch.negative=1   Point shape for the negative tests (default: 1, i.e. o symbol)
#' @param legend.corner="top"  Corner for the legend. When NULL, no legend is displayed.
#' @param full.legend=TRUE Plot additional indications on the legend (number of elements passing the different adjusted p-values)
#' @param legend.cex=1   Font size for the legend.
#' @param cex=0.6     Point size, passed to plot() and lines()
#' @param plot.points=TRUE Plot one point per row with the significance (sig=-log10(p-value)) as a function of the effect size.
#' @param tick.size=0.05 Height of the vertical lines denoting the boundaries of confidence intervals.
#' @param ... Additional parameters are passed to plot()
#' @return no return object
#'
#' @export
MAplot.MultiTestTable <- function(
  multitest.table,
  log2FC.col="log2FC",
  log2mean.col="log2mean",
  control.type = "padj",
  nb.tests = nrow(multitest.table),
  alpha = 0.05,
  log2FC.threshold= NULL,
  sort.by.pval = TRUE,
  plot.points = TRUE,
  xlab = "A = log2(mean signal)",
  ylab = "M = log2(fold-change)",
  density.colors = FALSE,
  col.points = '#888888',
  col.positive = '#4455DD',
  col.grid = '#AAAAAA',
  col.lines = 'blue',
  col.alpha = "darkred",
  lty.lines = "solid",
  lty.alpha = "dashed",
  lty.grid = "dotted",
  pch.positive = 3,
  pch.negative = 1,
  xlim = NULL,
  ylim = NULL,
  cex = 0.6,
  legend.corner = "top",
  full.legend = FALSE,
  legend.cex = 1,
  tick.size = 0.05,
  #                         Y.score = "sig", ## Score to plot on the Y axis. Supported: sig (default), p.value e.value
  ... ## additional parameters are passed to the plot function
) {

  if (is.null(col.positive)) {
    col.positive = NA
  }

  ## Identify features declared positive according to the specified alpha
  #positive <-  (multitest.table[,control.type] <= alpha)
  multitest.table$is.positive <- !is.na(multitest.table[,control.type] ) & (multitest.table[,control.type] <= alpha)
  multitest.table$status <- "negative"
  multitest.table$status[multitest.table$is.positive] <- "positive"
  # summary(multitest.table$is.positive)

  ## Apply threshold on effect size if required
  if (!is.null(log2FC.threshold)) {
    multitest.table$is.positive[abs(multitest.table[,log2FC.col]) < log2FC.threshold] <- FALSE
  }
  # table(multitest.table$is.positive)

  ## Sort the table before plotting, to ensure that positive points appear visible on top of negative points
  if (sort.by.pval) {
    order <- order(multitest.table[,control.type], decreasing = TRUE)
    multitest.table.sorted <- multitest.table[order,]
    if (length(col.points) == nrow(multitest.table)) {
      col.points <- col.points[order]
    }
  } else {
    multitest.table.sorted <- multitest.table
  }

  ## Select Y values depending on the control type (p.value, e.value, fdr)
  y.values.ori <- multitest.table.sorted[, log2FC.col]
  # hist(multitest.table.sorted[, control.type], breaks=100)
  # summary(multitest.table.sorted[, control.type])

  ## Fix a problem with infinite Y values resulting from 0 values for the control (p-value or derived stat)
  y.values <- y.values.ori
  # y.value.max <- 320 ## Maximal value for Y corresponds to the precision of floating point computation for the p-values (~ 1e-320)
  # y.values[is.infinite(y.values)] <- y.value.max
  # # hist(y.values, breaks = 100)


  ## Identify features with NA values
  # table(is.na(y.values))
  # na.values <- is.na(y.values)
  # sum(na.values)

  ## Get X values
  x.values <- multitest.table.sorted[, log2mean.col]

  ## Define limits of X and Y scales
  if (is.null(xlim)) {
    xlim <- c(min(x.values), max(x.values) * 1.2) ## Keep place for the legend
  }
  if (is.null(ylim)) {
    ylim <- range(y.values)
  }

  ## Define point colors and shapes.
  ## Note: the attributes col.points and pch.points can either be either
  ## a single value (for all points) or a  vector with one user-specified
  ## color per point.
  if (density.colors) {
    ## Compute status and density-specific color
    multitest.table.sorted$color <- featureColors(
      multitest.table.sorted$status,
      positions = data.frame(y.values,y.values))
  }  else {
    multitest.table.sorted$color <- col.points
    if (length(col.points) != nrow(multitest.table.sorted)) {
      if (!is.na(col.positive) & !is.null(col.positive)) {
        multitest.table.sorted[multitest.table.sorted$is.positive, "color"] <- rep(col.positive, length.out = sum(positive, na.rm = TRUE))
      }
    }
    # table(multitest.table.sorted$color)
  }

  ## Point character for each feature
  multitest.table.sorted$pch <- pch.negative
  if (!is.null(pch.positive)) {
    multitest.table.sorted[multitest.table.sorted$is.positive, "pch"] <- pch.positive
  }
  # table(multitest.table.sorted$pch)

  ## Identify the points above Y limits, to denote them by a different symbol and color
  below.ylim <- y.values.ori < ylim[1]
  y.values[below.ylim] <- ylim[1]
  multitest.table.sorted$pch[below.ylim] <- 6
  multitest.table.sorted$color[below.ylim] <- "purple"
  # table(below.ylim)

  above.ylim <- y.values > ylim[2]
  y.values[above.ylim] <- ylim[2]
  multitest.table.sorted$pch[above.ylim] <- 17
  multitest.table.sorted$color[above.ylim] <- "purple"
  # table(above.ylim)

  ################################################################
  ## Draw the MA plot

  ## Draw the plot frame (to have the grid below the points)
  plot(x.values,
       y.values,
       xlab = xlab,
       ylab = ylab,
       cex = cex,
       xlim = xlim,
       ylim = ylim,
       type = "n",
       panel.first = grid(lty = lty.grid, col = col.grid),
       ...)


  ## Plot the points
  if (plot.points) {
    points(x = x.values,
           y = y.values,
           cex  = cex,
           #col = "red",
           col = multitest.table.sorted$color,
           pch = multitest.table.sorted$pch)
  }


  ## Vertical line to denote the null position (corresponding to no difference between groups)
  abline(h = 0, col = "black", lty = "solid", lwd = 1)

  ## Add horizontal lines to denote the threshold on effect size
  if (!is.null(log2FC.threshold)) {
    abline(h = c(-log2FC.threshold, log2FC.threshold),
           col = col.lines,
           lwd = 1, lty = lty.lines) ## Vertical line to denote the threshold on effect size
  }


  ## Plot the legend
  if (!is.null(legend.corner)) {
    ## Store legend prameters in a table
    legend.table <- data.frame()

    ## Legend for alpha threshold
    if (!is.null(alpha)) {
      legend.table <- rbind(
        legend.table,
        data.frame("name" = "alpha",
                   "legend" = paste(sep = "", control.type," < ", signif(digits = 4, alpha)),
                   "lwd" = 0, "col" = "white", "pch" = -1, "lty" = lty.alpha))
    }


    ## Legend for threshold on effect size
    if (!is.null(log2FC.threshold)) {
      legend.table <- rbind(
        legend.table,
        data.frame("name" = "effect",
                   "legend" = paste(sep = "", "abs(", log2FC.col,") >= ", signif(digits = 4, log2FC.threshold)),
                   "lwd" = 2, "col" = col.lines, "pch" = -1, "lty" = lty.lines))
    }

    ## Legend for the points, depending on their status + display options
    nb.positives <- sum(positive, na.rm = TRUE)
    nb.negatives <- nb.tests - nb.positives
    legend.table <- rbind(
      legend.table,
      data.frame("name"="positive",
                 "legend"=paste(sep="", nb.positives, " positives"),
                 "lwd"=1, "col"= col.positive[1], "pch"=-1, "lty"="blank"))
    legend.table <- rbind(legend.table,
                          data.frame("name"="negative",
                                     "legend"=paste(sep="", nb.negatives, " negatives"),
                                     "lwd"=1, "col"= col.points[1], "pch"=-1, "lty"="blank"))

    row.names(legend.table) <- legend.table[, "name"]
    if (plot.points) {
      legend.table$pch <- -1
      legend.table["positive", "pch"] <- pch.positive
      legend.table["negative", "pch"] <- pch.negative
    }


    ## Plot a legend with additional info on the number of positives for different multiple testing corrections
    if (full.legend) {
      legend.pch <- c(as.vector(legend.table$pch), -1,-1,-1,-1)
      legend(legend.corner,
             c(paste(sep="","N=",nrow(multitest.table)),
               paste(sep="", sum(positive), " positives (",control.type," <= ", alpha, ")"),
               paste(sum(!positive), "negatives"),
               paste(sep="", "P-value <= ", alpha, ": ", sum(multitest.table[,"p.value"] <= alpha)),
               paste(sep="", "FDR <= ", alpha, ": ", sum(multitest.table[,"fdr"] <= alpha)),
               paste(sep="", "E-value <=", alpha, ": ", sum(multitest.table[,"e.value"] <= alpha)),
               paste(sep="", "E-value <= 1: ", sum(multitest.table[,"e.value"] <= 1))
             ),
             pch=legend.pch,
             cex=legend.cex,
             lwd=as.vector(legend.table$lwd),
             lty=c(as.vector(legend.table$lty), "blank", "blank", "blank", "blank"),
             col=c(as.vector(legend.table$col), "white","white","white","white"),
             bty="o", bg="white")
    } else {
      legend(legend.corner,
             #              c(paste(sep="", control.type," = ", alpha),
             #                paste(sep="", sum(positive), " positives"),
             #                paste(sum(!positive), "negatives")
             #              ),
             legend=as.vector(legend.table$legend),
             pch=as.vector(legend.table$pch),
             cex=legend.cex,
             lwd=as.vector(legend.table$lwd),
             lty=as.vector(legend.table$lty),
             col=as.vector(legend.table$col),
             bty="o", bg="white")
    }
  }
  # return(multitest.table.sorted) ## For debugging
}
