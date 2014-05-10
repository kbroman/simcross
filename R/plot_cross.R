## plot_cross.R

#  CCcolors
#
#' Collaborative Cross colors
#'
#' Get the vector of colors for the Collaborative Cross
#'
#' @return vector of eight colors
#'
#' @keywords color
#' @export
#' @examples
#' CCcolors()
CCcolors <-
function() {
  c("AJ"  =rgb(240,240,  0,maxColorValue=255),
    "B6"  =rgb(128,128,128,maxColorValue=255),
    "129" =rgb(240,128,128,maxColorValue=255),
    "NOD" =rgb( 16, 16,240,maxColorValue=255),
    "NZO" =rgb(  0,160,240,maxColorValue=255),
    "CAST"=rgb(  0,160,  0,maxColorValue=255),
    "PWK" =rgb(240,  0,  0,maxColorValue=255),
    "WSB" =rgb(144,  0,224,maxColorValue=255))
}  

#  plot_ind
#
#' Plot an individual
#'
#' Add an individual, as a pair of chromosomes, to a plot
#'
#' @param ind An individual object, as output by \code{\link{create.parent}}
#' or \code{\link{cross}}
#' @param center (x,y) vector for the center of the individual
#' @param length Length of chromosomes
#' @param width Width of chromsomes
#' @param gap Gap between chromsomes
#' @param col Vector of colors
#' @param border Color for border
#' @param lend Passed to \code{\link[graphics]{rect}}
#' @param ljoin Passed to \code{\link[graphics]{rect}}
#' @param allborders If TRUE, put borders around all segments
#' @param ... Additional arguments passed to rect()
#'
#' @return None.
#'
#' @importFrom graphics rect
#' @keywords hplot
#' @export
#' @examples
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 3:4)
#' kid <- cross(mom, dad)
#' plot(0,0, type="n", xlim=c(0, 100), ylim=c(0,100),
#'      xaxt="n", yaxt="n", xlab="", ylab="")
#' plot_ind(mom, c(25, 66), length=50)
#' plot_ind(dad, c(75, 66), length=50)
#' plot_ind(kid, c(50, 25), length=50)
plot_ind <-
function(ind, center, length=30, width=3, gap=3, col=CCcolors(),
         border="black", lend=1, ljoin=1, allborders=FALSE, ...)
{
  max_alleles <- max(sapply(ind, function(a) max(a[2,])))
  if(max_alleles > length(col))
    stop("Need more colors: length(col)=", length(col), " but max allele = ", max_alleles)

  for(i in 1:2) {
    sgn <- i*2-3
    chr <- ind[[i]]
    pos <- ind[[i]][1,]
    allele <- ind[[i]][2,]

    # rescale pos
    top <- center[2]+length/2
    bottom <- center[2]-length/2
    if(diff(par("usr")[3:4]) < 0) { # small values at top of figure
      tmp <- top
      top <- bottom
      bottom <- top
    }
    start <- pos[1]
    end <- pos[length(pos)]
    pos <- (pos-start)*(bottom-top)/(end-start) + top

    left <- center[1] + sgn*(gap/2+width)
    right <- center[1] + sgn*gap/2

    rect(left,  top, right, bottom,
         col=col[allele[1]],
         border=NA, lend=lend, ljoin=ljoin, ...)

    internalborder <- ifelse(allborders, border, NA)
    if(length(pos) > 2) {
      for(j in 2:length(pos))
        rect(rep(left, length(pos)-1),  pos[-length(pos)],
             rep(right, length(pos)-1), pos[-1],
             col=col[allele[-1]],
             border=internalborder, lend=lend, ljoin=ljoin, ...)
    }

    if(!is.na(border) && !is.null(border)) # draw border
      rect(center[1] + sgn*(gap/2+width), top,
           center[1] + sgn*gap/2,         bottom,
           col=NA, border=border, lend=lend, ljoin=ljoin, ...)
  }
}
