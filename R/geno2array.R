# geno2array:
#
#' Convert multi-allele genotype matrix to an 3-dimensional array.
#'
#' Convert a multi-allele genotype matrix as output by
#' \code{\link{convert2geno}} to a 3-dimensional array, of dimension
#' \code{n_ind} x \code{n_markers} x 2, with the last dimension
#' corresponding to the two alleles; we can't recover phase
#' information, though.
#'
#' @param geno Matrix of genotypes, as output by
#' \code{\link{convert2geno}}
#' @param n_founders Number of founder genotypes; if not provided,
#' inferred from contents.
#'
#' @return A 3-dimensional array, individual x marker x allele
#'
#' @export
#' @keywords utilities
#' @seealso \code{\link{get_geno}}, \code{\link{convert2geno}}
#'
#' @examples
#' # simulate 8-way RIL pedigree
#' tab <- sim_ril_pedigree(ngen=20, parents=1:8)
#' # simulate data from that pedigree
#' dat <- sim_from_pedigree(tab, L=100)
#' # marker map
#' map <- seq(0, 100, by=10)
#' names(map) <- paste0("m", map)
#' # convert dat to integer genotypes
#' geno <- convert2geno(dat, map)
#' # convert to 3-d array
#' arr <- geno2array(geno, 8)
geno2array <-
function(geno, n_founders)
{
    if(missing(n_founders) || is.null(n_founders))
        n_founders <- ceiling(log2(max(geno)))

    codes <- outer(1:n_founders, 1:n_founders, function(a,b) 2^(a-1)+2^(b-1))
    codes_lowtri <- codes[!upper.tri(codes)]
    mat <- col(codes)[!upper.tri(codes)]
    pat <- row(codes)[!upper.tri(codes)]

    output <- array(dim=c(nrow(geno), ncol(geno), 2))
    dimnames(output) <- list(rownames(geno), colnames(geno), NULL)
    output[,,1] <- mat[match(geno, codes_lowtri)]
    output[,,2] <- pat[match(geno, codes_lowtri)]

    output
}
