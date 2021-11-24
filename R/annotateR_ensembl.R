
# a is ensembl BED file with the ENSCAFT names
# b is the ensemblToGeneName file with the alternate gene names

#' @title annotateR_ensembl
#' @description Changes the gene names of an ensembl BED file.
#'
#' @param a ensembl BED file containing 'ENSCAFTxxx' nomenclature for genes.
#' @param b ensembl 'all fields' file containing the alternate gene names.
#'
#' @importFrom stats complete.cases
#'
#' @return Returns an annotated BED file containing both gene names.
#' @export
#'
#' @examples
#' \donttest{}

annotateR_ensembl <- function(a, b){
  a[] <- lapply(a, as.character)
  b[] <- lapply(b, as.character)
  refseq_transcripts <- a[,4]
  BED <- a
  new <- BED
  new[] <- lapply(BED, function(x) b[,2][match(x, b[,1])])
  BED[,4] <- new[,4]
  BED[] <- lapply(BED, as.character)
  BED[,13] <- refseq_transcripts
  BED <- BED[stats::complete.cases(BED[,4]),]
}
