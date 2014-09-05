
context("get_geno")

test_that("simple test of get_geno", {

    dat <- list("1"=list(mat=rbind(c(0, 9.22, 83.6, 100), c(1,1,2,1)),
                         pat=rbind(c(0, 48.7, 100), c(2, 2, 1))),
                "2"=list(mat=rbind(c(0, 18.64, 55.21, 100), c(2, 2, 1, 2)),
                         pat=rbind(c(0, 31.55, 95.42, 100), c(2, 2, 1, 2))),
                "3"=list(mat=rbind(c(0, 5.5, 39.4, 100), c(2, 2, 1, 2)),
                         pat=rbind(c(0, 89.1, 100), c(1, 1, 2))))

    ex_geno <- rbind(c(2, 2), c(1, 2), c(1, 1))
    geno <- get_geno(dat, 30)
    expect_equivalent(geno, ex_geno)

})

