library(Rmmquant)
library(testthat)

context("Checking base function")

test_that("Running init", {
    expect_null(init())
})

test_that("Running simple files", {
    dir     <- system.file("extdata", package="Rmmquant", mustWork = TRUE)
    gtfFile <- file.path(dir, "test.gtf")
    samFile <- file.path(dir, "test.sam")
    init()
    setGtfFileName(gtfFile)
    addReadsFileName(samFile)
    table <- start()
    expect(is.matrix(table));
    m <- matrix(1, 1, 1)
    dimnames(m) <- list(c("geneA"), c("test"))
    expect_equal(table, m)
})