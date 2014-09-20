
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

