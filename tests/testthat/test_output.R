library(Rmmquant)
library(testthat)
library(GenomicRanges)

context("Output tests")
dir <- system.file("extdata", package="Rmmquant", mustWork = TRUE)

test_that("Running default test", {
    gtfFile     <- file.path(dir, "test.gtf")
    samFile     <- file.path(dir, "test.sam")
    object      <- RmmquantRun(gtfFile, samFile)
    table       <- counts(object)
    stats       <- stats(object)
    expect_is(table, "matrix")
    expect_is(stats, "data.frame")
})

test_that("Running empty output", {
    gtfFile <- file.path(dir, "test.gtf")
    samFile <- file.path(dir, "test.sam")
    object  <- RmmquantRun(annotationFile=gtfFile,
                           readsFiles    =samFile,
                           overlap       =10000)
    table   <- counts(object)
    expect(is.matrix(table));
    expect_equal(dim(table), c(0, 1))
})

test_that("Running error", {
    gtfFile <- file.path(dir, "test.gtf")
    expect_error(RmmquantRun(annotationFile=gtfFile))
})