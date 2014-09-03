# sim_do_pedigree
#
#' Simulate diversity outcross pedigree
#'
#' Simulate a diversity outcross pedigree (a table of individual, mom, dad, sex)
#'
#' @param ngen Number of generations of outbreeding
#' @param npairs Number of breeding pairs at each generation
#' @param ccgen Vector of same length as npairs, with the number of generations for each CC line
#' @param design How to choose crosses: either random but avoiding siblings, or completely at random
#'
#' @return A matrix with six columns: individual ID, mother ID, father
#' ID, sex, generation, and 1/0 indicator for whether DO or pre-DO.
#' Founders have \code{0} for mother and father ID. Sex is coded 0 for
#' female and 1 for male.
#'
#' @keywords datagen
#' @export
#'
#' @examples
#' tab <- sim_do_pedigree(8, 30, rep(6, 30))
sim_do_pedigree <-
function(ngen=12, npairs=30, ccgen=rep(0, npairs), nkids_per=5, design=c("nosib", "random"))
{
    if(length(ccgen)==1) ccgen <- rep(ccgen, npairs)
    stopifnot(length(ccgen) == npairs)
    stopifnot(all(ccgen >= 0))

    # initial generation : need to double each strain; 1-8 are the mothers; 9-16 are the corresponding dads
    id <- 1:16
    mom <- rep(0, 16)
    dad <- rep(0, 16)
    sex <- rep(0:1, each=8)
    gen <- rep(0, 16)
    result <- cbind(id=id, mom=mom, dad=dad, sex=sex, gen=gen)

    cur_nind <- 16
    lastgen <- NULL
    for(i in 1:npairs) {
        # random cross among the 8 strains
        parents <- sample(1:8)
        parents[seq(2, 8, by=2)] <- parents[seq(2, 8, by=2)] + 8

        tab <- sim_ril_pedigree(ccgen[i], selfing=FALSE, parents=parents, firstind=cur_nind+1)[-(1:8),]
        result <- rbind(result, tab)
        lastgen <- c(lastgen, tab[tab[,"gen"] == max(tab[,"gen"]),1])
        cur_nind <- max(tab[,1])
    }
    result <- cbind(result, do=rep(0, nrow(result)))

    moms <- lastgen[seq(1, length(lastgen), by=2)]
    dads <- lastgen[seq(2, length(lastgen), by=2)]
    for(i in 1:ngen) {
        dads <- sample(dads)
        while(design=="nosib" && any(dads - moms == 1)) { # sample until no sibs
            dads <- sample(dads)
        }

        if(i < ngen) {
            id <- 1:(npairs*2)+cur_nind

            mom <- c(mom, rep(moms, each=2))
            dad <- c(dad, rep(dads, each=2))
            sex <- c(sex, rep(c(0,1), npairs))
            gen <- c(gen, rep(i, npairs*2))

            moms <- id[seq(1, length(id), by=2)]
            dads <- id[seq(2, length(id), by=2)]
        }
        else { # last generation: expand
            id <- 1:(npairs*nkids_per)+cur_nind

            mom <- c(mom, rep(moms, each=nkids_per))
            dad <- c(dad, rep(dads, each=nkids_per))
            sex <- c(sex, rep(c(0,1), ceiling(npairs*nkids_per/2))[1:(npairs*nkids_per)])
            gen <- c(gen, rep(i, npairs*nkids_per))
        }
        cur_nind <- max(id)

        result <- rbind(result, cbind(id, mom, dad, sex, gen, rep(1, length(id))))
    }

    result
}
