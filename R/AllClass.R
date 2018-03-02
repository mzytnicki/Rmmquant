#' Function.
#'
#' @param annotationFile The annotation file
#' @param readsFiles     The reads files
#' @param sampleNames    The reads files
#' @param overlap        The reads files
#' @param strand         The reads files
#' @param sort           The reads files
#' @param countThreshold The reads files
#' @param mergeThreshold The reads files
#' @param useGeneName    The reads files
#' @param quiet          The reads files
#' @param progress       The reads files
#' @param nThreads       The reads files
#' @param format         The reads files
#' @param nOverlapDiff   The reads files
#' @param pcOverlapDiff  The reads files
#' @param lazyload       Usual for S4 functions.
#' @return             An \code{sRNADiff} object.
#' @examples
#' dir         <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' data        <- read.csv(file.path(dir, "data.csv"))
#' gtfFile     <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf.gz")
#' annotation  <- readWholeGenomeAnnotation(gtfFile)
#' bamFiles    <- file.path(dir, data$FileName)
#' replicates  <- data$SampleName
#' conditions  <- factor(data$Condition)
#' exp         <- sRNADiffExp(annotation, bamFiles, replicates, conditions)
#'
#' @export
RmmquantRun <- function(annotationFile,
                        readsFiles,
                        sampleNames   =NULL,
                        overlap       =NULL,
                        strands       =NULL,
                        sorts         =NULL,
                        countThreshold=NULL,
                        mergeThreshold=NULL,
                        useGeneName   =NULL,
                        quiet         =NULL,
                        progress      =NULL,
                        nThreads      =NULL,
                        format        =NULL,
                        nOverlapDiff  =NULL,
                        pcOverlapDiff =NULL,
                        lazyload      =FALSE) {
    parameters <- new(RmmquantParameters)
    if ((! is.character(annotationFile)) | (! is.atomic(annotationFile))) 
        stop(paste0("Annotation file should be a file name."))
    parameters.setGtfFileName(annotationFile)
    if (! is.character(readsFiles))
        stop(paste0("Read files should be a (set of) file name(s)."))
    if (is.atomic(readsFiles))
        parameters.addReadsFileName(readsFiles)
    else {
        for (readsFile in readsFiles) {
            parameters.addReadsFileName(readsFiles)
        }
    }
    if (! is.null(sampleNames)) {
        if (is.atomic(sampleNames))
            parameters.addName(sampleNames)
        else {
            for (sampleName in sampleNames) {
                parameters.addName(sampleName)
            }
        }
    }
    if (! is.null(overlap)) {
        parameters.setOverlap(overlap)
    }
    if (! is.null(strands)) {
        if (is.atomic(strands))
            parameters.addStrand(strands)
        else {
            for (strand in strands) {
                parameters.addStrand(strand)
            }
        }
    }
    if (! is.null(sorts)) {
        if (is.atomic(sorts))
            parameters.addSort(sorts)
        else {
            for (sort in sorts) {
                parameters.addSort(sort)
            }
        }
    }
    if (! is.null(countThreshold)) {
        parameters.setCountThreshold(countThreshold)
    }
    if (! is.null(mergeThreshold)) {
        parameters.setMergeThreshold(mergeThreshold)
    }
    if (! is.null(useGeneName)) {
        parameters.setGeneName(useGeneName)
    }
    if (! is.null(quiet)) {
        parameters.setQuiet(quiet)
    }
    if (! is.null(progress)) {
        parameters.setProgress(progress)
    }
    if (! is.null(nThreads)) {
        parameters.setNThreads(nThreads)
    }
    if (! is.null(format)) {
        parameters.setFormat(format)
    }
    if (! is.null(nOverlapDiff)) {
        parameters.setNOverlapDifference(nOverlapDiff)
    }
    if (! is.null(pcOverlapDiff)) {
        parameters.setPcOverlapDifference(pcOverlapDiff)
    }
    return(rStart(parameters))
}


#' Example constructor
#' @return An \code{srnadiff} object
#'
#' @examples
#' exp <- sRNADiffExample()
#'
#' @export
sRNADiffExample <- function() {
    dir        <- system.file("extdata", package="srnadiff", mustWork=TRUE)
    data       <- read.csv(file.path(dir, "data.csv"))
    gtfFile    <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf.gz")
    annotation <- readWholeGenomeAnnotation(gtfFile)
    bamFiles   <- file.path(dir, data$FileName)
    replicates <- data$SampleName
    conditions <- factor(data$Condition)
    object     <- sRNADiffExp(annotation, bamFiles, replicates, conditions)
    return(object)
}
