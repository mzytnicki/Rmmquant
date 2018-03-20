#' Rmmquant object validation function
#'
#' @param object A \code{RmmquantClass} object.
#' @return       \code{TRUE}, if succeed, otherwise a \code{character}.
validateRmmquant <- function(object) {
    errors <- character()
    if ((! is.character(object@annotationFile)) |
        (! is.atomic(object@annotationFile))) {
        errors <- c(errors, "Annotation file should be one file name.")
    }
    if ((object@annotationFile == "") & (length(object@genomicRanges) == 0) &
        (length(object@genomicRangesList) == 0)) {
        errors <- c(errors, "Missing annotation.")
    }
    if (((object@annotationFile != "") &
         (length(object@genomicRanges) != 0)) |
        ((object@annotationFile != "") &
         (length(object@genomicRangesList) != 0)) |
        ((length(object@genomicRanges) != 0) &
         (length(object@genomicRangesList) != 0))) {
        errors <- c(errors, "At least annotation inputs given.")
    }
    return(if (length(errors) == 0) TRUE else errors)
}

#' An S4 class for Rmmquant.
#'
#' @slot annotationFile    The annotation file
#' @slot readsFiles        The reads files
#' @slot genomicRanges     The annotation, in a
#'                              \code{GenomicRanges} format.
#' @slot genomicRangesList The annotation, in a
#'                              \code{GenomicRangesList} format.
#' @slot sampleNames       The name of the samples
#' @slot overlap           The minimum number of overlapping base pairs to
#'                              declare a match.
#' @slot strands           Whether annotation of the same strand should be
#'                              considered.
#' @slot sorts             Whether the files are sorted.
#' @slot countThreshold    The reads files
#' @slot mergeThreshold    The reads files
#' @slot printGeneName     Whether the (vernacular) gene name is reported.
#' @slot quiet             Shut Rmmquant up.
#' @slot progress          Print the progress of the tool.
#' @slot nThreads          The number of threads.
#' @slot formats           The format of the reads files (SAM or BAM).
#' @slot nOverlapDiff      Difference of overlap between a primary map and a
#'                             secondary map.
#' @slot pcOverlapDiff     Ratio of overlap between a primary map and a
#'                             secondary map.
#' @slot counts            A matrix of the counts.
#' @slot stats             A data frame of the quantification statistics.
setClass("RmmquantClass",
         representation(
             annotationFile    ="character",
             readsFiles        ="character",
             genomicRanges     ="GRanges",
             genomicRangesList ="GRangesList",
             sampleNames       ="character",
             overlap           ="numeric",
             strands           ="character",
             sorts             ="logical",
             countThreshold    ="numeric",
             mergeThreshold    ="numeric",
             printGeneName     ="logical",
             quiet             ="logical",
             progress          ="logical",
             nThreads          ="numeric",
             formats           ="character",
             nOverlapDiff      ="numeric",
             pcOverlapDiff     ="numeric",
             counts            = "matrix",
             stats             = "data.frame"),
         prototype(counts=matrix(), stats=data.frame()),
         validity = validateRmmquant)

#' Main Rmmquant function, and constructor of the class.
#'
#' @param annotationFile    The annotation file
#' @param readsFiles        The reads files
#' @param genomicRanges     The annotation, in a
#'                              \code{GenomicRanges} format.
#' @param genomicRangesList The annotation, in a
#'                              \code{GenomicRangesList} format.
#' @param sampleNames       The name of the samples
#' @param overlap           The minimum number of overlapping base pairs to
#'                              declare a match.
#' @param strands           Whether annotation of the same strand should be
#'                              considered.
#' @param sorts             Whether the files are sorted.
#' @param countThreshold    The reads files
#' @param mergeThreshold    The reads files
#' @param printGeneName     Whether the (vernacular) gene name is reported.
#' @param quiet             Shut Rmmquant up.
#' @param progress          Print the progress of the tool.
#' @param nThreads          The number of threads.
#' @param formats           The format of the reads files (SAM or BAM).
#' @param nOverlapDiff      Difference of overlap between a primary map and a
#'                              secondary map.
#' @param pcOverlapDiff     Ratio of overlap between a primary map and a
#'                              secondary map.
#' @param lazyload          Usual for S4 functions.
#' @return                  A numeric matrix.
#' @examples
#' dir <- system.file("extdata", package="Rmmquant", mustWork = TRUE)
#' gtfFile <- file.path(dir, "test.gtf")
#' samFile <- file.path(dir, "test.sam")
#' table <- RmmquantRun(gtfFile, samFile)
#'
#' @export
RmmquantRun <- function(annotationFile    ="",
                        readsFiles,
                        genomicRanges     =GRanges(),
                        genomicRangesList =GRangesList(),
                        sampleNames       =character(0),
                        overlap           =NA_integer_,
                        strands           =character(0),
                        sorts             =logical(0),
                        countThreshold    =NA_integer_,
                        mergeThreshold    =NA_real_,
                        printGeneName     =FALSE,
                        quiet             =TRUE,
                        progress          =FALSE,
                        nThreads          =1,
                        formats           =character(0),
                        nOverlapDiff      =NA_integer_,
                        pcOverlapDiff     =NA_real_,
                        lazyload          =FALSE) {
    object <- new("RmmquantClass",
                  annotationFile   =annotationFile,
                  readsFiles       =readsFiles,
                  genomicRanges    =genomicRanges,
                  genomicRangesList=genomicRangesList,
                  sampleNames      =sampleNames,
                  overlap          =overlap,
                  strands          =strands,
                  sorts            =sorts,
                  countThreshold   =countThreshold,
                  mergeThreshold   =mergeThreshold,
                  printGeneName    =printGeneName,
                  quiet            =quiet,
                  progress         =progress,
                  nThreads         =nThreads,
                  nOverlapDiff     =nOverlapDiff,
                  pcOverlapDiff    =pcOverlapDiff)
    data <- rcpp_Rmmquant(annotationFile,
                          readsFiles,
                          genomicRanges,
                          genomicRangesList,
                          sampleNames,
                          overlap,
                          strands,
                          sorts,
                          countThreshold,
                          mergeThreshold,
                          printGeneName,
                          quiet,
                          progress,
                          nThreads,
                          formats,
                          nOverlapDiff,
                          pcOverlapDiff)
    counts <- data$counts
    if (is.null(colnames(counts))) {
        colnames(counts) <- if (length(sampleNames) == 0) seq(ncol(counts))
                            else sampleNames
    }
    stats <- as.data.frame(data$stats,
                           row.names=make.names(colnames(counts),unique=TRUE))
    object@counts <- counts
    object@stats  <- stats
    return(object)
}

#' Example constructor
#' @return An \code{RmmquantClass} object
#'
#' @examples
#' example <- RmmquantExample()
#'
#' @export
RmmquantExample <- function() {
    dir     <- system.file("extdata", package="Rmmquant", mustWork = TRUE)
    gtfFile <- file.path(dir, "test.gtf")
    samFile <- file.path(dir, "test.sam")
    return(RmmquantRun(gtfFile, samFile))
}

#' Overloading the show method
#' @param object An \code{RmmquantClass} object.
#' @return       A description of the object.
#'
#' @examples
#' example <- RmmquantExample()
#' example
#'
#' @export
setMethod("show",
          signature(object="RmmquantClass"),
          function(object) {
                cat("Object of class Rmmquant.\n",
                    "Counts (length = ",
                    length(object@counts),
                    "):\n", sep="")
                show(object@counts)
                cat("Stats:\n")
                show(object@stats)
            }
)

#' Get the counts
#' @rdname counts-method
#' @param  object  An \code{RmmquantClass} object.
#' @return         The count matrix
#'
#' @examples
#' example <- RmmquantExample()
#' counts(example)
#'
#' @export
setGeneric("counts", function(object) { standardGeneric("counts") } )

#' @rdname counts-method
#' @export
setMethod("counts",
          signature(object="RmmquantClass"),
          function(object) { object@counts })

#' Get the stats
#' @rdname stats-method
#' @param  object  An \code{RmmquantClass} object.
#' @return         The stats data frame
#'
#' @examples
#' example <- RmmquantExample()
#' stats(example)
#'
#' @export
setGeneric("stats", function(object) { standardGeneric("stats") } )

#' @rdname stats-method
#' @export
setMethod("stats",
          signature(object="RmmquantClass"),
          function(object) { object@stats })