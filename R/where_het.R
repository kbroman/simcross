## where_het.R

# where_het
#
#' Find heterozygous regions
#'
#' Find regions of heterozygosity in an individual
#'
#' @param ind An individual object, as output be \code{\link{create_parent}}
#' or \code{\link{cross}}
#' @return A matrix with two columns; each row indicates the start
#' and end of a region where the individual is heterozygous
#'
#' @export
#' @examples
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 1:2)
#' child <- cross(mom, dad)
#' where_het(child)
where_het <-
function(ind)
{
  if(ncol(ind$mat)==ncol(ind$pat) && all(ind$mat == ind$pat)) {
    return(NULL)
  }
  u <- sort(unique(c(ind$mat[1,],ind$pat[1,])))
  het <- NULL
  for(i in 2:length(u)) {
    mat <- ind$mat[,ind$mat[1,] >= u[i],drop=FALSE]
    mat <- mat[2,1]

    pat <- ind$pat[,ind$pat[1,] >= u[i],drop=FALSE]
    pat <- pat[2,1]

    if(mat!=pat) { # heterozygous
      if(is.null(het)) het <- cbind(u[i-1],u[i])
      else het <- rbind(het,c(u[i-1],u[i]))
    }
  }

  # clean up
  if(nrow(het) > 1) {
    keep <- rep(TRUE,nrow(het))
    for(j in 2:nrow(het)) {
      if(het[j,1] == het[j-1,2]) {
        het[j,1] <- het[j-1,1]
        keep[j-1] <- FALSE
      }
    }
    het <- het[keep,,drop=FALSE]
  }
  het
}
