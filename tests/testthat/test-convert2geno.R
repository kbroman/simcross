
context("convert2geno")

test_that("convert2geno works for 2-allele case", {

    # data as list with alleles in intervals
    #    (as produced by sim_from_pedigree)
    dat <- list("1" = list(mat = list(alleles = c(2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(1.63593688048422, 50.6488258950412,
                                                    54.441562271677, 72.7342922240496,
                                                    80.7659703074023, 100)),
                           pat = list(alleles = c(2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(8.94253242295235, 43.9252462703735,
                                                    82.6503853080794, 90.0170957203954,
                                                    91.0740879364312, 100))),
                "2" = list(mat = list(alleles = c(1L, 2L, 1L, 2L, 1L),
                                      locations = c(1.72073859721422, 2.41351707372814,
                                                    72.7342922240496, 80.7659703074023, 100)),
                           pat = list(alleles = c(2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(8.94253242295235, 43.9252462703735,
                                                    54.1366793680936, 75.1830249791965,
                                                    82.6503853080794, 90.0170957203954,
                                                    91.0740879364312, 100))),
                "3" = list(mat = list(alleles = c(1L, 2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(1.72073859721422, 2.41351707372814,
                                                    50.6488258950412, 54.441562271677,
                                                    89.9764276109636, 97.2906407434493, 100)),
                           pat = list(alleles = c(2L, 1L, 2L, 1L),
                                      locations = c(8.94253242295235, 43.9252462703735,
                                                    97.3880127770826, 100))),
                "4" = list(mat = list(alleles = c(2L, 1L, 2L, 1L),
                                      locations = c(1.63593688048422, 72.7342922240496,
                                                    80.7659703074023, 100)),
                           pat = list(alleles = c(2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(20.5045772017911, 38.2374913664535,
                                                    48.7418097676709, 75.1830249791965,
                                                    82.6503853080794, 88.3461387362331,
                                                    97.3880127770826, 100))))

    # marker map
    map <- seq(0, 100, by=5)
    names(map) <- paste0("marker", seq(along=map))

    # expected genotype matrix (force integers)
    expected <- rbind(c(3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1),
                      c(2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1),
                      c(2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1),
                      c(3, 2, 2, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 3, 1, 2, 2, 1))
    expected <- matrix(as.integer(expected), nrow=4)
    dimnames(expected) <- list(as.character(1:4), paste0("marker", 1:21))

    expect_equal(expected, convert2geno(dat, map))




})


test_that("convert2geno works for 8-allele case", {

    # data as list with alleles in intervals
    #    (as produced by sim_from_pedigree)
    dat <- list("1" = list(mat = list(alleles=c(8L, 2L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(8.20657967124134, 18.4184361249208,
                                                  55.8109006844461, 86.3054546294734,
                                                  88.7538079405203, 90.0979985482991, 100)),
                           pat = list(alleles = c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  32.0876344339922, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100))),
                "2" = list(mat = list(alleles=c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  32.7756949001923, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100)),
                           pat = list(alleles = c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  29.4709661742672, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100))),
                "3" = list(mat = list(alleles=c(8L, 2L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(8.20657967124134, 18.4184361249208,
                                                  55.8109006844461, 86.3054546294734,
                                                  88.7538079405203, 90.0979985482991, 100)),
                           pat = list(alleles = c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  32.7756949001923, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100))),
                "4" = list(mat = list(alleles=c(8L, 2L, 8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(8.20657967124134, 12.613788805902,
                                                  21.6286054579541, 27.2114308550954,
                                                  32.0876344339922, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100)),
                           pat = list(alleles = c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  29.4709661742672, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100))))


    # marker map
    map <- seq(0, 100, by=5)
    names(map) <- paste0("marker", seq(along=map))

    # expected genotype matrix (force integers)
    expected <- rbind(c(256, 256, 130, 130, 256, 130, 129, 256, 256, 256, 256, 256,
                        32, 32, 32, 32, 32, 32, 32, 8, 8),
                      c(256, 256, 256, 256, 256, 4, 129, 256, 256, 256, 256, 256, 32,
                        32, 32, 32, 32, 32, 32, 8, 8),
                      c(256, 256, 130, 130, 256, 130, 129, 256, 256, 256, 256, 256,
                        32, 32, 32, 32, 32, 32, 32, 8, 8),
                      c(256, 256, 130, 256, 256, 4, 129, 256, 256, 256, 256, 256, 32,
                        32, 32, 32, 32, 32, 32, 8, 8))
    expected <- matrix(as.integer(expected), nrow=4)
    dimnames(expected) <- list(as.character(1:4), paste0("marker", 1:21))

    expect_equal(expected, convert2geno(dat, map))

})
