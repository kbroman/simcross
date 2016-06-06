# sim_ril_pedigree
#
#' Generate a ril pedigree
#'
#' Generate a pedigree for multi-way recombinant inbred lines (a table
#' of individual, mom, dad, sex)
#'
#' @param ngen Number of generations of inbreeding
#' @param selfing If TRUE, use selfing
#' @param parents Vector of length 2, 4, or 8, indicating the parents'
#' IDs
#' @param firstind Positive integer to assign to the first child
#'
#' @return A data frame with five columns: individual ID, mother ID,
#' father ID, sex, and generation.  Founders have \code{0} for mother
#' and father ID. Sex is coded 0 for female and 1 for male.
#'
#' @export
#' @seealso \code{\link{sim_from_pedigree}},
#' \code{\link{sim_ail_pedigree}}, \code{\link{sim_do_pedigree}},
#' \code{\link{sim_4way_pedigree}}
#' @keywords datagen
#'
#' @examples
#' tab <- sim_ril_pedigree(7)
sim_ril_pedigree <-
    function(ngen=20, selfing=FALSE, parents=1:8, firstind=max(parents)+1)
{
    nparents <- length(parents)
    stopifnot(nparents == 2 || nparents==4 || nparents==8) # 2, 4, or 8 parents

    # pedigree up to the last outcross
    if(nparents==2) { # 2 parents
        id <- parents
        mom <- dad <- c(0, 0)
        sex <- c(0, 1)
        gen <- c(0, 0)
        nextid <- firstind
    }
    else if(nparents==4) { # 4 parents
        id <- c(parents, 0:1+firstind)
        mom <- c(0, 0, 0, 0, parents[1], parents[3])
        dad <- c(0, 0, 0, 0, parents[2], parents[4])
        sex <- rep(c(0,1), nparents/2 + 1)
        gen <- c(0, 0, 0, 0, 1, 1)
        nextid <-  firstind + 2
    }
    else { # 8 parents
        id <- c(parents, 0:5+firstind)
        mom <- c(rep(0, nparents), parents[seq(1, 8, by=2)], firstind,   firstind+2)
        dad <- c(rep(0, nparents), parents[seq(2, 8, by=2)], firstind+1, firstind+3)
        sex <- rep(c(0, 1), nparents/2 + 3)
        gen <- rep(0:2, c(8, 4, 2))
        nextid <- firstind + 6
    }

    if(selfing) {
        # just need one individual to start
        mom <- c(mom, id[length(id)-1])
        dad <- c(dad, id[length(id)])
        id <- c(id, nextid)
        sex <- c(sex, 0) # call the hermaphrodites female
        gen <- c(gen, max(gen)+1)
        nextid <- nextid + 1

        if(ngen > 0) {
            last_outcross <- max(gen)
            for(g in 1:ngen) {
                mom <- c(mom, id[length(id)])
                dad <- c(dad, id[length(id)])
                id <- c(id, nextid)
                sex <- c(sex, 0)
                gen <- c(gen, last_outcross+g)
                nextid <- nextid + 1
            }
        }
    }
    else {
        mom <- c(mom, rep(id[length(id)-1], 2))
        dad <- c(dad, rep(id[length(id)], 2))
        id <- c(id, 0:1 + nextid)
        sex <- c(sex, 0, 1)
        gen <- c(gen, rep(max(gen)+1, 2))
        nextid <- nextid + 2

        if(ngen > 0) {
            last_outcross <- max(gen)
            for(g in 1:ngen) {
                mom <- c(mom, rep(id[length(id)-1], 2))
                dad <- c(dad, rep(id[length(id)], 2))
                id <- c(id, nextid+0:1)
                sex <- c(sex, 0, 1)
                gen <- c(gen, rep(last_outcross+g, 2))
                nextid <- nextid + 2
            }
        }
    }

    data.frame(id=id, mom=mom, dad=dad, sex=sex, gen=gen)
}
