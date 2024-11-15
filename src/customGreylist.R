# https://www.bioconductor.org/packages/release/bioc/vignettes/GreyListChIP/inst/doc/GreyList-demo.pdf
# https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
# https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html#looping-over-grangeslist-objects
#' @param args[1] path to input bam files
#' @param args[2] path to chromsizes file
#' @param args[3] outdir

args <- commandArgs(TRUE)

library(GreyListChIP)