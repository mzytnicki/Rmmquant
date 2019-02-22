suppressWarnings(suppressMessages(library(Rmmquant)))
suppressWarnings(suppressMessages(library(testthat)))
suppressWarnings(suppressMessages(library(GenomicRanges)))

context("Case tests")
dir <- system.file("extdata", package="Rmmquant", mustWork = TRUE)

test_that("Running two chromosomes test", {
    gr1         <- GRanges(seqnames="chr1", ranges=IRanges(c(1000, 3000), width=1001, names="geneA"), strand="+")
    gr2         <- GRanges(seqnames="chr2", ranges=IRanges(c(5000, 8000), width=1001, names="geneB"), strand="+")
    grl         <- GRangesList("geneA"=gr1, "geneB"=gr2)
    samFile     <- file.path(dir, "test1.sam")
    se          <- RmmquantRun(genomicRangesList=grl, readsFiles=samFile)
    table       <- assays(se)$counts
    m           <- matrix(c(1, 1), 2, 1)
    dimnames(m) <- list(c("geneA", "geneB"), c("test1"))
    expect_equal(table, m)
})

test_that("Running multi-mapped read test", {
    gr1         <- GRanges(seqnames="chr1", ranges=IRanges(c(1000, 3000), width=1001, names="geneA"), strand="+")
    gr2         <- GRanges(seqnames="chr2", ranges=IRanges(c(5000, 8000), width=1001, names="geneB"), strand="+")
    grl         <- GRangesList("geneA"=gr1, "geneB"=gr2)
    samFile     <- file.path(dir, "test2.sam")
    se          <- RmmquantRun(genomicRangesList=grl, readsFiles=samFile)
    table       <- assays(se)$counts
    m           <- matrix(c(1), 1, 1)
    dimnames(m) <- list(c("geneA--geneB"), c("test2"))
    expect_equal(table, m)
})

test_that("Running antisense unknown test", {
    gr1         <- GRanges(seqnames="chr1", ranges=IRanges(c(1000, 3000), width=1001, names="geneA"), strand="+")
    gr2         <- GRanges(seqnames="chr2", ranges=IRanges(c(5000, 8000), width=1001, names="geneB"), strand="-")
    grl         <- GRangesList("geneA"=gr1, "geneB"=gr2)
    samFile     <- file.path(dir, "test1.sam")
    se          <- RmmquantRun(genomicRangesList=grl, readsFiles=samFile)
    table       <- assays(se)$counts
    m           <- matrix(c(1, 1), 2, 1)
    dimnames(m) <- list(c("geneA", "geneB"), c("test1"))
    expect_equal(table, m)
})

test_that("Running antisense forward test", {
    gr1         <- GRanges(seqnames="chr1", ranges=IRanges(c(1000, 3000), width=1001, names="geneA"), strand="+")
    gr2         <- GRanges(seqnames="chr2", ranges=IRanges(c(5000, 8000), width=1001, names="geneB"), strand="-")
    grl         <- GRangesList("geneA"=gr1, "geneB"=gr2)
    samFile     <- file.path(dir, "test1.sam")
    se          <- RmmquantRun(genomicRangesList=grl, readsFiles=samFile, strands="F")
    table       <- assays(se)$counts
    m           <- matrix(c(1), 1, 1)
    dimnames(m) <- list(c("geneA"), c("test1"))
    expect_equal(table, m)
})

test_that("Running antisense reverse test", {
    gr1         <- GRanges(seqnames="chr1", ranges=IRanges(c(1000, 3000), width=1001, names="geneA"), strand="+")
    gr2         <- GRanges(seqnames="chr2", ranges=IRanges(c(5000, 8000), width=1001, names="geneB"), strand="-")
    grl         <- GRangesList("geneA"=gr1, "geneB"=gr2)
    samFile     <- file.path(dir, "test1.sam")
    se          <- RmmquantRun(genomicRangesList=grl, readsFiles=samFile, strands="R")
    table       <- assays(se)$counts
    m           <- matrix(c(1), 1, 1)
    dimnames(m) <- list(c("geneB"), c("test1"))
    expect_equal(table, m)
})

test_that("Running antisense both test", {
    gr1         <- GRanges(seqnames="chr1", ranges=IRanges(c(1000, 3000), width=1001, names="geneA"), strand="-")
    gr2         <- GRanges(seqnames="chr2", ranges=IRanges(c(5000, 8000), width=1001, names="geneB"), strand="*")
    grl         <- GRangesList("geneA"=gr1, "geneB"=gr2)
    samFile     <- file.path(dir, "test1.sam")
    se          <- RmmquantRun(genomicRangesList=grl, readsFiles=samFile, strands="F")
    table       <- assays(se)$counts
    m           <- matrix(c(1), 1, 1)
    dimnames(m) <- list(c("geneB"), c("test1"))
    expect_equal(table, m)
})