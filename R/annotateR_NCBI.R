#' @title annotateR_NCBI
#' @description Changes the gene names of an NCBI BED file.
#'
#' @param a NCBI BED file containing 'XR/XM' nomenclature for genes.
#' @param b NCBI 'all fields' file containing the alternate gene names.
#'
#' @return Returns an annotated BED file containing both gene names.
#' @export
#'
#' @importFrom stats complete.cases
#'
#' @examples
#' \dontrun{}
#'

annotateR_NCBI <- function(a, b){
  a[] <- lapply(a, as.character)
  b[] <- lapply(b, as.character)
  refseq_transcripts <- a[,4]
  BED <- a
  new <- BED
  new[] <- lapply(BED, function(x) b$name2[match(x, b$name)])
  BED[,4] <- new[,4]
  BED[] <- lapply(BED, as.character)
  BED[,13] <- refseq_transcripts
  BED <- BED[complete.cases(BED[,4]),]
}
