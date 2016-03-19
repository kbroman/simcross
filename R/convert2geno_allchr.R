#' Convert continuous allele information into marker genotypes for
#' multiple chromosomes
#'
#' Wrap up of convert2geno to adequate multiple chromosomes.
#'
#' @param xodat The sort of detailed genotype/crossover data generated
#' by \code{\link{sim_from_pedigree_allchr}}
#' @param map marker locations, a list with elements for each
#' chromosome
#' @param id ids for which individuals genotypes is desired
#' @param f.geno Optional matrix (size \code{n_markers} x 8) of
#' founder genotypes. If coded as 1/2, results are 1/2/3 genotypes. If
#' coded as A/T/G/C/N, results are A/T/G/C/N/H genotypes. If coded as
#' letters A-H for the 8 founders, results are two-letter genotypes
#' AA-HH with 36 possible values.
#' @param shift_map If TRUE, shift genetic map to start at 0
#' @param return.matrix If FALSE, the result is a list of length
#' \code{n_chrs}, otherwise it is converted into a matrix if size
#' \code{length(id)} x \code{n_markers}.
#'
#' @return If \code{founder_geno} is provided or there are just two
#' founders, the result is a numeric matrix of genotypes, individuals
#' x markers, with genotypes 1/2/3 codes for 11/12/22 genotypes. If
#' there are more than two founders and \code{founder_geno} are
#' letters, the result is a character matrix, too.
#'
#' If \code{founder_geno} is not provided and there are more than two
#' founders, the result is a 3-dimensional array, individuals x
#' markers x alleles, with the third dimensional corresponding to the
#' maternal and paternal allele.
#'
#' @export
#' @keywords utilities
#' @seealso \code{\link{convert2geno}}
#'
#' @examples
#' library(qtl)
#' # marker map
#' map <- sim.map(len=rep(100, 19), n.mar=10, include.x=FALSE)
#' # simulate AIL pedigree
#' tab <- sim_ail_pedigree(12, 30)
#' # simulate data from that pedigree
#' dat <- sim_from_pedigree_allchr(tab, map)
#' names(map) <- paste0("marker", seq(along=map))
#' # convert data to marker genotypes
#' id <- which(tab[, "gen"]==12)
#' geno <- convert2geno_allchr(dat, map, id)
#'
convert2geno_allchr <- function(xodat, map, id=NULL, f.geno=NULL,
                                  return.matrix=TRUE, shift_map=FALSE){

  stopifnot(length(map) == length(xodat))
  if(is.null(id)) id <- 1:length(xodat[[1]])

  if(is.null(f.geno)){
    geno <- NULL
    for(chr in 1:length(map)){
      geno[[chr]] <- convert2geno(xodat=xodat[[chr]][id], map=map[[chr]])
    }
    names(geno) <- names(map)
  }else{

    ## create a list of founder geno, one element for each chr.
    if(is.list(f.geno)){
      ## check that f.geno is a list of same length as map
      stopifnot(length(f.geno)==length(map))
    }else{  ## f.geno is a matrix.
      nm.chr <- sapply(map, length)
      nm.map <- sum(nm.chr)
      nm <- nrow(f.geno)
      if(nm.map!=nm) stop("f.geno should have ", nm.map, " columns but has ", nm)

      f.geno <- cbind(f.geno, f.geno)
      f.geno.list <- list()
      for(chr in 1:length(map)){
        if(chr==1)
            f.geno.list[[chr]] <- t(f.geno[1:nm.chr[chr], ])
        else
            f.geno.list[[chr]] <- t(f.geno[sum(nm.chr[1:(chr-1)]) + (1:nm.chr[chr]), ])
      }
      f.geno <- f.geno.list
    }
    geno <- NULL
    for(chr in 1:length(map)){
      geno[[chr]] <- convert2geno(xodat=xodat[[chr]][id], map=map[[chr]],
                                  founder_geno=f.geno[[chr]])
    }
    names(geno) <- names(map)
  }

  if(!return.matrix){ ## return the list directly
    return(geno)
  }else{ ## cbind all elememts to make a matrix.
    if(is.matrix(geno[[1]])){
      geno.mat <- NULL
      for(i in 1:length(geno)){
        geno.mat <- cbind(geno.mat, geno[[i]])
      }
    }else{ ## this could be a arrary when f.geno is missing...
      ## individuals x markers x alleles,
      geno.mat1 <- NULL
      geno.mat2 <- NULL
      for(i in 1:length(geno)){
        geno.mat1 <- cbind(geno.mat1, geno[[i]][,,1])
        geno.mat2 <- cbind(geno.mat2, geno[[i]][,,2])
      }
      geno.mat <- array(NA, dim=c(dim(geno.mat1), 2))
      geno.mat[,,1] <- geno.mat1
      geno.mat[,,2] <- geno.mat2
    }
    return(geno.mat)
  }
}
