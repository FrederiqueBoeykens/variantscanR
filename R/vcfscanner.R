#' @title vcfscanneR
#' @description Uploads a single sample VCF file (from a multi sample VCF file) into R.
#'
#' @param x Character string containing the name of the vcf file or file path.
#' @param sample Character string containing the name of the sample of interest. The sample name has to match with that of the VCF file.
#'
#' @return Returns a single sample vcf file ready for downstream analysis.
#' @export
#' @importFrom vcfR read.vcfR
#' @importFrom dplyr "%>%"
#'
#'@examples
#'vcf_file <- system.file("extdata", "SNPs.recode.subset.rename.vcf.gz")
#'\dontrun{
#'vcf <- vcfscanneR(vcf_file, "sample_9")
#'}


vcfscanneR <- function(x, sample){
  vcf <- vcfR::read.vcfR(x)
  fix <- vcf@fix
  fix <- as.data.frame(fix, stringsAsFactors = FALSE)
  gt <- vcf@gt
  gt <- as.data.frame(gt, stringsAsFactors = FALSE)
  fix <- cbind(fix, gt[,1])
  names(fix)[9] <- "FORMAT"
  gt <- gt %>% dplyr::select(dplyr::all_of(sample))
  vcf <- cbind.data.frame(fix,gt)
  vcf_names <- c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample")
  names(vcf) <- vcf_names
  return(vcf)
}
