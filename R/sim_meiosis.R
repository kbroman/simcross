## sim_meiosis.R

# create_parent
#
#' Create a parent object
#'
#' Create a parent object
#'
#' @param L chromosome length in cM
#' @param allele vector of integers for alleles, of length 1 or 2
#'
#' @return A list with two components, for the individual's two
#' chromosomes.  Each is a data frame with two columns: alleles in
#' chromosome intervals (as integers), and locations of the
#' right endpoints of those intervals.
#'
#' @keywords datagen
#' @export
#' @seealso \code{\link{cross}}, \code{\link{meiosis}}
#'
#' @examples
#' create_parent(100, 1)
#' create_parent(100, 1:2)
create_parent <-
function(L, allele=1L)
{
    if(length(allele) == 1) allele <- rep(allele,2)
    if(length(allele) != 2)
        stop("allele should be of length 1 or 2")
    if(!is.integer(allele)) {
        if(!is.numeric(allele))
            stop("allele should be a vector with 1 or 2 integers")
        allele <- as.integer(allele)
    }

    list(mat=data.frame(alleles=allele[1], locations=L),
         pat=data.frame(alleles=allele[2], locations=L))
}


# check that data for an individual conforms to expected format
check_individual <-
function(ind, tol=1e-12)
{
    # list with two components, named "mat" and "pat"
    if(!is.list(ind) || length(ind)!=2 ||
       !all(names(ind) == c("mat", "pat")))
        stop('ind should be list with "mat" and "pat"')

    # check each chromosome
    for(i in 1:2) {
        # data.frame with "alleles" and "locations" components
        if(!is.data.frame(ind[[i]]) || ncol(ind[[i]])!=2 ||
           !all(names(ind[[i]]) == c("alleles", "locations")))
            stop('chromosome should be data frame with "alleles" and "locations"')

        alleles <- ind[[i]]$alleles
        locations <- ind[[i]]$locations

        # locations numeric
        if(!is.numeric(locations))
            stop("locations should be numeric")

        # alleles integer
        if(!is.integer(alleles))
            stop("alleles should be integers")

        # locations non-decreasing
        if( length(locations) > 1 && min(diff(locations)) < 0 )
            stop("locations should be non-decreasing")
    }

    # start and end positions are the same
    starts <- vapply(ind, function(a) a$location[1], 0.0)
    ends <- vapply(ind, function(a) a$location[length(a$location)], 0.0)
    stopifnot(abs(diff(starts)) < tol, abs(diff(ends)) < tol)

    TRUE
}

# note: function sim_crossovers in C++

# meiosis
#
#' Simulate meiosis
#'
#' Output a random meiotic product from an input individual.
#'
#' @param parent An individual object, as output by
#' \code{\link{create_parent}} or \code{\link{cross}}
#' @param m interference parameter for chi-square model
#' @param p Proportion of chiasmata coming from no-interference process.
#'
#' @return A data frame with two columns: alleles in
#' chromosome intervals (as integers), and locations of the
#' right endpoints of those intervals.
#'
#' @keywords datagen
#' @export
#' @seealso \code{\link{create_parent}}, \code{\link{cross}}
#'
#' @examples
#' ind <- create_parent(100, 1:2)
#' prod <- meiosis(ind)
meiosis <-
function(parent, m=10, p=0)
{
    tol <- 1e-12
    L <- max(parent$mat$locations)
    if(abs(L - max(parent$pat$locations)) > tol)
        stop("parent's two chromosomes are not the same length")

    product <- sim_crossovers(L, m, p)
    cur_allele <- sample(1:2,1) # first allele

    if(length(product)==0) return(parent[[cur_allele]])
    else {
        product <- c(-1, product)
        loc <- alle <- NULL
        for(i in 2:length(product)) {
            interval <- which(parent[[cur_allele]]$locations >= product[i-1] &
                              parent[[cur_allele]]$locations < product[i])
            loc <- c(loc, parent[[cur_allele]]$locations[interval])
            alle <- c(alle, parent[[cur_allele]]$alleles[interval])

            toright <- parent[[cur_allele]]$locations>=product[i]
            next_allele <- parent[[cur_allele]]$alleles[toright][1]

            loc <- c(loc, product[i])
            alle <- c(alle, next_allele)

            cur_allele <- 3 - cur_allele # 1 -> 2 or 2 -> 1

        }
        toright <- which(parent[[cur_allele]]$locations >= max(product))
        loc <- c(loc, parent[[cur_allele]]$locations[toright])
        alle <- c(alle, parent[[cur_allele]]$alleles[toright])
    }

    clean_meiosis(loc, alle)
}

# clean up meiotic product
#  - no two adjacent alleles need be the same
clean_meiosis <-
function(loc, alle)
{
    if(length(loc) > 1) {
        keep <- rep(TRUE, length(loc))
        for(i in 1:(length(loc)-1))
            if(alle[i] == alle[i+1])
                keep[i] <- FALSE
        loc <- loc[keep]
        alle <- alle[keep]
    }

    data.frame(alleles=alle, locations=loc)
}


# cross
#
#' Cross two individuals
#'
#' Simulate the cross of two individuals to create a
#' single progeny
#'
#' @param mom An individual object, as produced by
#' \code{\link{create_parent}} or this function.
#' @param dad An individual object, as produced by
#' \code{\link{create_parent}} or this function.
#' @param m interference parameter for chi-square model
#' @param p proption of crossovers coming from no-interference process
#' @param xchr If TRUE, simulate X chromosome
#' @param male If TRUE, simulate a male (matters only if
#' \code{xchr=TRUE})
#'
#' @return A list with two components, for the individual's two
#' chromosomes.  Each is a data frame with two columns: alleles in
#' chromosome intervals (as integers), and locations of the
#' right endpoints of those intervals.
#'
#' @keywords datagen
#' @export
#' @seealso \code{\link{create_parent}}, \code{\link{meiosis}}
#'
#' @examples
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 1:2)
#' child <- cross(mom, dad)
cross <-
function(mom, dad, m=10, p=0, xchr=FALSE, male=FALSE)
{
    if(!xchr) {
        return(list(mat=meiosis(mom,m,p),
                    pat=meiosis(dad,m,p)))
    }
    else {
        if(male)
            return(list(mat=meiosis(mom,m,p),
                        pat=dad$pat))
        else
            return(list(mat=meiosis(mom,m,p),
                        pat=dad$mat))
    }
}
