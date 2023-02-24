#' BED file from UCSC table browser.
#'
#' Canis familiares (CanFam3.1, Sep 2011) annotation data from
#' the Genes and Gene Predictions table displayed Browser Extensible Data format
#' that provides a flexible way to define the data lines that are displayed in an annotation track.
#' The file was downloaded from the UCSC table browser website.
#'
#' @docType data
#'
#' @usage data(dog_BED_NCBI)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{chrom}{The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671)}
#'  \item{chromStart}{The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0}
#'  \item{chromEnd}{The ending position of the feature in the chromosome or scaffold}
#'  \item{name}{Defines the name of the BED line}
#'  \item{score}{A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray)}
#'  \item{strand}{Defines the strand. Either "." (=no strand) or "+" or "-"}
#'  \item{thickStart}{The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position}
#'  \item{thickEnd}{The ending position at which the feature is drawn thickly (for example the stop codon in gene displays)}
#'  \item{itemRgb}{An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line}
#'  \item{blockCount}{The number of blocks (exons) in the BED line}
#'  \item{blockSizes}{A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount}
#'  \item{blockStarts}{A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount}
#' }
#' @references This data set was downloaded from the UCSC table browser website <https://genome.ucsc.edu/cgi-bin/hgTables>
#' @keywords datasets
#' @examples
#'
#' data(dog_BED_NCBI)
#' head(dog_BED_NCBI)
"dog_BED_NCBI"
