#' All Fields BED file from UCSC table browser.
#'
#' Canis familiares (CanFam3.1, Sep 2011) annotation data from
#' the Genes and Gene Predictions table displayed in a tab-separated format
#' suitable for import into spreadsheets and relational databases
#' downloaded from the UCSC table browser website.
#'
#' @docType data
#'
#' @usage data(dog_allfields_NCBI)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{bin}{Indexing field to speed chromosome range queries}
#'  \item{name}{Name of gene (usually transcript_id from GTF)}
#'  \item{chrom}{Reference sequence chromosome or scaffold}
#'  \item{strand}{+ or - for strand}
#'  \item{txStart}{Transcription start position (or end position for minus strand item)}
#'  \item{txEnd}{Transcription end position (or start position for minus strand item)}
#'  \item{cdsStart}{Coding region start (or end position for minus strand item)}
#'  \item{cdsEnd}{Coding region end (or start position for minus strand item)}
#'  \item{exonCount}{Number of exons}
#'  \item{exonStarts}{Exon start positions (or end positions for minus strand item)}
#'  \item{exonEnds}{Exon end positions (or start positions for minus strand item)}
#'  \item{score}{score}
#'  \item{name2}{Alternate name (e.g. gene_id from GTF)}
#'  \item{cdsStartStat}{Status of CDS start annotation (none, unknown, incomplete, or complete)}
#'  \item{cdsEndStat}{Status of CDS end annotation (none, unknown, incomplete, or complete)}
#'  \item{exonFrames}{xon frame {0,1,2}, or -1 if no frame for exon}
#' }
#' @references This data set was downloaded from the UCSC table browser website <https://genome.ucsc.edu/cgi-bin/hgTables>
#' @keywords datasets
#' @examples
#'
#' data(dog_allfields_NCBI)
#' head(dog_allfields_NCBI)
"dog_allfields_NCBI"
