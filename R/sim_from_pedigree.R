# sim_from_pedigree
#
#' Simulate a pedigree
#'
#' Simulate genotypes along one chromosome for a pedigree
#'
#' @param pedigree Matrix describing a pedigree, with first four
#' columns being individual ID, mom ID, dad ID, and sex (female as
#' \code{0}, male as \code{1}).
#' @param L Length of chromosome in cM.
#' @param xchr If TRUE, simulate X chromosome.
#' @param m Crossover interference parameter, for chi-square model
#' (m=0 corresponds to no interference).
#' @param obligate.chiasma If TRUE, 4-strand bundle has at least one chiasma.
#'
#' @return A list whose each component is the data for one individual,
#' as produced by the \code{\link{cross}} function.  Those results are
#' a list of two matrices, corresponding to the maternal and paternal
#' chromosomes. The chromosomes are represented as a matrix with two
#' rows: a vector of positions, and the alleles at each position (and
#' in interval to left).
#'
#' @export
#' @keywords datagen
#'
#' @examples
#' # simulate AIL pedigree
#' tab <- sim_ail_pedigree(12, 30)
#' # simulate data from that pedigree
#' dat <- sim_from_pedigree(tab)
sim_from_pedigree <-
function(pedigree, L=100, xchr=FALSE, m=10, obligate.chiasma=TRUE)
{
    if(length(unique(pedigree[,1])) != nrow(pedigree))
        stop("IDs must be unique")

    result <- vector("list", nrow(pedigree))
    names(result) <- as.character(pedigree[,1])

    for(i in 1:nrow(pedigree)) {
        if(pedigree[i,2]==0 || pedigree[i,3]==0) # founder
            result[[i]] <- create_parent(L, allele=pedigree[i,1])
        else {
            mom <- which(pedigree[,1]==pedigree[i,2])
            dad <- which(pedigree[,1]==pedigree[i,3])

            if(length(mom) < 1 || length(dad) < 1) {
                print(pedigree[i,,drop=FALSE])
                stop("parents not found")
            }
            if(mom >= i || dad >= i) {
                print(pedigree[i,,drop=FALSE])
                stop("Pedigree problem: parents follow individual")
            }

            result[[i]] <- cross(result[[mom]], result[[dad]],
                                 m=m, obligate.chiasma=obligate.chiasma,
                                 xchr=xchr, male=pedigree[i,4]==1)
        }
    }
    result
}
