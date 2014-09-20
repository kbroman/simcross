
context("sim_meiosis")

test_that("test of create_parent", {

    expected_p1 <- list(mat=data.frame(alleles=1L, locations=100.0),
                        pat=data.frame(alleles=1L, locations=100.0))

    expect_equal(create_parent(100, 1), expected_p1)

    expected_p2 <- list(mat=data.frame(alleles=2L, locations=100.0),
                        pat=data.frame(alleles=2L, locations=100.0))

    expect_equal(create_parent(100, 2), expected_p2)

    expected_f1 <- list(mat=data.frame(alleles=1L, locations=100.0),
                        pat=data.frame(alleles=2L, locations=100.0))

    expect_equal(create_parent(100, 1:2), expected_f1)

    # numeric alleles converted to integers
    expect_equal(create_parent(100, c(1.9, 2.9)), expected_f1)

    # non-numeric alleles not allowed
    expect_error(create_parent(100, "a"))

    # can't have more than 2 alleles
    expect_error(create_parent(100, 1:3))

    # must have at least 1 allele
    expect_error(create_parent(100, NULL))
    expect_error(create_parent(100, numeric(0)))

})


test_that("test of check_individual", {

    expected_p1 <- list(mat=data.frame(alleles=1L, locations=100.0),
                        pat=data.frame(alleles=1L, locations=100.0))

    expected_p2 <- list(mat=data.frame(alleles=2L, locations=100.0),
                        pat=data.frame(alleles=2L, locations=100.0))

    expected_f1 <- list(mat=data.frame(alleles=1L, locations=100.0),
                        pat=data.frame(alleles=2L, locations=100.0))

    # should be clean
    expect_true(check_individual(expected_p1))
    expect_true(check_individual(expected_p2))
    expect_true(check_individual(expected_f1))

    # no data
    expect_error(check_individual(NULL))
    expect_error(check_individual(list(mat=NULL, pat=NULL)))

    # alleles not integers
    z <- list(mat=data.frame(alleles=1, locations=100),
              pat=data.frame(alleles=2, locations=100))
    expect_error(check_individual(z))
    z <- list(mat=data.frame(alleles="a", locations=100),
              pat=data.frame(alleles="b", locations=100))
    expect_error(check_individual(z))

    # wrong orders
    z <- list(mat=data.frame(locations=100, alleles=1L),
              pat=data.frame(alleles=1L, locations=100))
    expect_error(check_individual(z))
    z <- list(mat=data.frame(alleles=1L, locations=100),
              pat=data.frame(locations=100, alleles=1L))
    expect_error(check_individual(z))
    z <- list(pat=data.frame(alleles=2L, locations=100),
              mat=data.frame(alleles=1L, locations=100))
    expect_error(check_individual(z))

    # locations not in order
    z <- list(mat=data.frame(alleles=c(2L, 1L), locations=c(5, 100)),
              pat=data.frame(alleles=c(1L, 2L), locations=c(100, 85)))
    expect_error(check_individual(z))
    z <- list(mat=data.frame(alleles=c(2L, 1L), locations=c(500, 100)),
              pat=data.frame(alleles=c(1L, 2L), locations=c(100, 105)))
    expect_error(check_individual(z))
})

