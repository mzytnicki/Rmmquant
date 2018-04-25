library(Rmmquant)
library(testthat)
library(GenomicRanges)

context("Output tests")
dir <- system.file("extdata", package="Rmmquant", mustWork = TRUE)

test_that("Running default test", {
    gtfFile     <- file.path(dir, "test.gtf")
    samFile     <- file.path(dir, "test.sam")
    se          <- RmmquantRun(gtfFile, samFile)
    table       <- assays(se)[[1]]
    stats       <- colData(se)
    expect_is(se, "SummarizedExperiment")
    expect_is(table, "matrix")
    expect_is(stats, "DataFrame")
})

test_that("Running empty output", {
    gtfFile <- file.path(dir, "test.gtf")
    samFile <- file.path(dir, "test.sam")
    se      <- RmmquantRun(annotationFile=gtfFile,
                           readsFiles    =samFile,
                           overlap       =10000)
    table   <- assays(se)$counts
    expect(is.matrix(table));
    expect_equal(dim(table), c(0, 1))
})