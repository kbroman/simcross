
context("geno2array")

test_that("geno2array works in multi-allele case", {

    set.seed(63730516)

    # sim 8-way RIL
    tab <- sim_ril_pedigree(ngen=20, parents=1:8)
    dat <- sim_from_pedigree(tab, L=100)

    # genotypes every 10 cM
    map <- seq(0, 100, by=10)
    names(map) <- paste0("m", map)
    geno <- convert2geno(dat, map)

    # array version
    arr <- geno2array(geno, 8)

    # compare to result of get_geno
    for(i in seq(along=map)) {
        g <- get_geno(dat, map[i])
        g <- t(apply(g, 1, sort))
        colnames(g) <- NULL
        expect_equal(g, arr[,i,])
    }

})
