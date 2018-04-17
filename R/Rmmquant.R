#' Rmmquant: RNA-Seq multi-mapping Reads Quantification Tool
#'
#' Counts the number of reads per gene.
#'
#' @docType package
#' @name Rmmquant
#'
#' @author Matthias Zytnicki, \email{matthias.zytnicki@@inra.fr}
#'
#' @import S4Vectors
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import devtools
#' @import BiocStyle
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import org.Hs.eg.db
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 lfcShrink
#' @importFrom methods is
#' @importFrom methods new
#' @importFrom Rcpp evalCpp
#' @useDynLib Rmmquant
NULL
