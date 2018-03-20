library(Rmmquant)
library(testthat)
library(GenomicRanges)

context("Parameter tests")
dir <- system.file("extdata", package="Rmmquant", mustWork = TRUE)

test_that("Running full parameter test", {
    gtfFile     <- file.path(dir, "test.gtf")
    samFile     <- file.path(dir, "test.sam")
    object      <- RmmquantRun(annotationFile   =gtfFile,
                               genomicRanges    =GRanges(),
                               genomicRangesList=GRangesList(),
                               readsFiles       =c(samFile),
                               sampleNames      =c("test"),
                               overlap          =1,
                               strands          =c("U"),
                               sorts            =c(TRUE),
                               countThreshold   =1,
                               mergeThreshold   =0.0,
                               printGeneName    =FALSE,
                               quiet            =TRUE,
                               progress         =FALSE,
                               nThreads         =1,
                               formats          =c("SAM"),
                               nOverlapDiff     =30,
                               pcOverlapDiff    =0.5)
    table       <- counts(object)
    m           <- matrix(1, 1, 1)
    dimnames(m) <- list(c("geneA"), c("test"))
    expect_equal(table, m)
})

test_that("Running sample name test", {
    gtfFile     <- file.path(dir, "test.gtf")
    samFile     <- file.path(dir, "test.sam")
    object      <- RmmquantRun(annotationFile   =gtfFile,
                               genomicRanges    =GRanges(),
                               genomicRangesList=GRangesList(),
                               readsFiles       =c(samFile),
                               sampleNames      ="name")
    table       <- counts(object)
    m           <- matrix(1, 1, 1)
    dimnames(m) <- list(c("geneA"), c("name"))
    expect_equal(table, m)
})

test_that("Running unsorted reads test", {
    gtfFile     <- file.path(dir, "test.gtf")
    samFile     <- file.path(dir, "test.sam")
    object      <- RmmquantRun(annotationFile   =gtfFile,
                               genomicRanges    =GRanges(),
                               genomicRangesList=GRangesList(),
                               readsFiles       =c(samFile),
                               sorts            =c(FALSE))
    table       <- counts(object)
    m           <- matrix(1, 1, 1)
    dimnames(m) <- list(c("geneA"), c("test"))
    expect_equal(table, m)
})

test_that("Running two threads test", {
    gtfFile     <- file.path(dir, "test.gtf")
    samFile     <- file.path(dir, "test.sam")
    object      <- RmmquantRun(annotationFile   =gtfFile,
                               genomicRanges    =GRanges(),
                               genomicRangesList=GRangesList(),
                               readsFiles       =c(samFile, samFile),
                               nThreads         =2)
    table       <- counts(object)
    m           <- matrix(c(1, 1), 1, 2)
    dimnames(m) <- list(c("geneA"), c("test", "test"))
    expect_equal(table, m)
})