library(Rmmquant)
library(testthat)
library(GenomicRanges)

context("Input tests")
dir <- system.file("extdata", package="Rmmquant", mustWork = TRUE)

test_that("Running SAM file test", {
    gtfFile     <- file.path(dir, "test.gtf")
    samFile     <- file.path(dir, "test.sam")
    se          <- RmmquantRun(gtfFile, samFile)
    table       <- assays(se)$counts
    stats       <- colData(se)
    m           <- matrix(1, 1, 1)
    dimnames(m) <- list(c("geneA"), c("test"))
    expect_equal(table, m)
    expect_equal(stats$n.hits,                     2)
    expect_equal(stats$n.uniquely.mapped.reads,    1)
    expect_equal(stats$n.ambiguously.mapped.hits,  0)
    expect_equal(stats$n.non.uniquely.mapped.hits, 0)
    expect_equal(stats$n.unassigned.hits,          1)
})
    
test_that("Running two SAM files test", {
    gtfFile     <- file.path(dir, "test.gtf")
    samFile     <- file.path(dir, "test.sam")
    se          <- RmmquantRun(gtfFile, c(samFile, samFile))
    table       <- assays(se)$counts
    m           <- matrix(c(1, 1), 1, 2)
    dimnames(m) <- list(c("geneA"), c("test", "test.1"))
    expect_equal(table, m)
})
    
test_that("Running BAM file test", {
    gtfFile     <- file.path(dir, "test.gtf")
    samFile     <- file.path(dir, "test.bam")
    se          <- RmmquantRun(gtfFile, samFile)
    table       <- assays(se)$counts
    m           <- matrix(1, 1, 1)
    dimnames(m) <- list(c("geneA"), c("test"))
    expect_equal(table, m)
})
    
test_that("Running genomic ranges test", {
    samFile     <- file.path(dir, "test.sam")
    gr          <- GRanges(seqnames="chr1",
                           ranges=IRanges(1000, width=3000, names="geneA"),
                           strand="+")
    se          <- RmmquantRun(genomicRanges=gr, readsFiles=c(samFile))
    table       <- assays(se)$counts
    m           <- matrix(1, 1, 1)
    dimnames(m) <- list(c("geneA"), c("test"))
    expect_equal(table, m)
})

test_that("Running genomic ranges list test", {
    samFile     <- file.path(dir, "test.sam")
    gr1         <- GRanges(seqnames="chr1",
                           ranges=IRanges(c(1000, 3000),
                                          width=1001, names="geneA"),
                           strand="+")
    gr2         <- GRanges(seqnames="chr2",
                           ranges=IRanges(c(1000, 3000),
                                          width=1001, names="geneB"),
                           strand="+")
    grl         <- GRangesList("geneA"=gr1, "geneB"=gr2)
    se          <- RmmquantRun(genomicRangesList=grl, readsFiles=c(samFile))
    table       <- assays(se)$counts
    m           <- matrix(1, 1, 1)
    dimnames(m) <- list(c("geneA"), c("test"))
    expect_equal(table, m)
})