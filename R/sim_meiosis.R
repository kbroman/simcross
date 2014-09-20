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
#' @seealso \code{\link{cross}}, \code{\link{sim_meiosis}}
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

# note: function sim_crossovers in src/sim_meiosis.cpp
# note: function sim_meiosis in src/sim_meiosis.cpp


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
#' @seealso \code{\link{create_parent}}, \code{\link{sim_meiosis}}
#'
#' @examples
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 1:2)
#' child <- cross(mom, dad)
cross <-
function(mom, dad, m=10, p=0, xchr=FALSE, male=FALSE)
{
    if(!xchr) {
        return(list(mat=sim_meiosis(mom,m,p),
                    pat=sim_meiosis(dad,m,p)))
    }
    else {
        if(male)
            return(list(mat=sim_meiosis(mom,m,p),
                        pat=dad$pat))
        else
            return(list(mat=sim_meiosis(mom,m,p),
                        pat=dad$mat))
    }
}
