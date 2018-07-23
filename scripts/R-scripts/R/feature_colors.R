#' @title define feature colors for scatter plots as a function of either their score or the local density
#' @author Jacques van Helden
#' @param type="2col" type of palette. Suppored: 2col (2-colors gradient), "gray" (gray scale), "dens" (proportional to local density). 
#' Palette types 2col and gray resuire to specify a score associated to each feature. 
#' Density palette requires to specify X and Y coordinates. 
#' @oaram score=NULL score associated to each feature
#' @param x=NULL feature abcsissas for density palettes
#' @param y=NULL feature ordinates for density palettes. Must have the same length as x.
#' @param levels=100 number of palette levels
#' @param end.colors=c('blue','red')  Colors and the two extemes of the palette.
FeatureColors <- function(palette.type = "gray",
                          scores=NULL,
                          x=NULL,
                          y=NULL,
                          levels = 100,
                          end.colors = c('blue', 'red')) {
  
  if (palette.type == "dens") {
    if (is.null(x) || is.null(y)) {
      stop("FeatureColors()\tpalette type 'dens' requires to specity x and y parameters. ")
    }
    feature.colors <- densCols(x = x, y = y)
  } else if (palette.type %in% c("gray", "grey", "2col")) {
    if (is.null(scores)) {
      stop("FeatureColors()\tpalette types 'gray' and '2col' require to specity scores. ")
    }
    
    if (palette.type %in% c("gray", "grey")) {
      feature.palette <- gray.colors(levels, start = 0.3, end = 1, gamma = 2.2, alpha = NULL)

    } else if (palette.type == "2col") {
      feature.palette <- colorRampPalette(c(end.colors[1], end.colors[2]))(levels)
    }
    feature.colors <- feature.palette[as.numeric(cut(scores, breaks = levels))]
  } else {
    stop("FeatureColors()\tInvalid palette type. Supported: 'gray', '2col', 'dens' ")
  }
  return(feature.colors)
}
