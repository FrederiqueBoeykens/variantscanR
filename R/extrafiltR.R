#' @title extrafiltR
#' @description Filtering function that takes the single sample VCF file and filters for present variants in genes of
#' interest. These genes of interest are obtained from the user-defined variants file.
#'
#' @param vcf Single sample VCF file, uploaded with the vcfscanneR function.
#' @param variants_file File containing variants of interest, processed by the pRocess function.
#' @param BED_file_annot Annotated BED file created with the annotateR_NCBI/ensembl function.
#'
#' @return A dataframe will be obtained containing all the variants within the genes of interest.
#'
#' @export
#'
#' @importFrom svMisc progress
#' @importFrom R3port html_list html_combine
#' @importFrom stringr fixed str_replace_all
#'
#' @examples
#' \dontrun{}


extrafiltR <- function(vcf, variants_file, BED_file_annot){
  cicip <- colnames(variants_file)[5]
  if (cicip == "Check"){
    variants_file <- variants_file[,-c(4:6)]
  } else {variants_file <- variants_file}
  vcf$pos <- as.numeric(vcf$pos)
  variants_file$Start <- as.numeric(variants_file$Start)
  df <- data.frame()
  filter_list <- list() #Create an empty list to store the non-empty filter dataframes
  for (r in 1:nrow(variants_file)){ #for every row in the OMIA_file
    nrrow <- nrow(variants_file)
    x <- variants_file[r,1] #get the chromosome we're looking at  str(chr)
    gene <- variants_file[r,6] #get the gene within the chromosome  str(num)
    subsetgene <- BED_file_annot[BED_file_annot[,4] == gene,] #subset on only this gene
    # Check if the generated dataframe is empty
    if (nrow(subsetgene) == 0){
      # If it is, break out of the loop and move on to the next row
      next
    }
    a <- as.numeric(min(subsetgene[,2]))
    b <- as.numeric(max(subsetgene[,3]))
    subsetvcffile <- vcf[vcf$chrom == x,] #subset on only this chromosome
    subsetvcffile <- subsetvcffile[subsetvcffile$pos > a,]
    subsetvcffile <- subsetvcffile[subsetvcffile$pos > b,]
    filter <- unique(subsetvcffile)
    output1 <- data.frame()
    filter <- filter[,-3]
    filter <- filter[,-5]
    filter <- filter[,-6]
    samplename <- names(filter)[7]
    if (nrow(filter) > 0){
      filter_list[[r]] <- filter #Store the non-empty filter dataframe in the list
    }
  }
  filter_all <- do.call(rbind, filter_list)
  filter_all <- unique(filter_all)
  convert_row <- function(row) {
    format <- unlist(strsplit(as.character(row[6]), ":"))
    ex <- which(format == "GT")
    snp <- unlist(strsplit(as.character(row[7]), ":"))
    snp <- snp[ex]
    snp <- gsub("\\|", "/", snp)
    snp_sep <- unlist(strsplit(snp, "/"))
    snp1 <- snp_sep[1]
    snp2 <- snp_sep[2]
    if (snp1 == "." ){snp1 <- "A call cannot be made for this sample at this given locus"} else {
      snp1 <- as.numeric(snp1)
    }
    if (snp2 == "." ){snp2 <- "A call cannot be made for this sample at this given locus"} else {
      snp2 <- as.numeric(snp2)
    }
    zygosity <- if(snp1 == "A call cannot be made for this sample at this given locus" | snp2 == "A call cannot be made for this sample at this given locus" ) {"Zygosity could not be determined"} else if (snp1 == snp2){
      "Homozygous"} else {"Heterozygous"}
    allele_alt <- unlist(strsplit(as.character(row[4]), ","))
    if (is.numeric(snp1) == FALSE){allele1 <- snp1} else if(snp1 == 0){allele1 <- row[3]} else {allele1 <- allele_alt[snp1]}
    if (is.numeric(snp2) == FALSE){allele2 <- snp2} else if(snp2 == 0){allele2 <- row[3]} else {allele2 <- allele_alt[snp2]}
    output <- c(row[1:3], allele1, allele2, zygosity)
    return(output)
  }# Use apply() to apply the function to each row of filter
  output_list <- apply(filter_all, 1, convert_row)
  output_df <- as.data.frame(output_list)
  output_df <- t(output_df)
  output_df <- as.data.frame(output_df)
  colnames(output_df) <- c("Chromosome", "Location", "Wild_type", "First_allele", "Second_allele", "Zygosity")
  rownames(output_df) <- NULL
  result <- unique(output_df)
  result[] <- lapply(result, as.character)
  result$Location <- as.numeric(result$Location)
  rownames(result) <- NULL
  whatif <- result %>%
    filter(if_any(First_allele:Second_allele, ~ .x != Wild_type))
  nrrows <- nrow(whatif)
  generep <- as.data.frame(rep(c(gene), times = nrrows))
  colnames(generep) <- "gene"
  total <- cbind(whatif[,1], generep, whatif[,2:ncol(whatif)])
  df <- rbind.data.frame(df, total)
  colnames(df) <- c("Chromosome", "gene", "Location", "Wild_type", "First_allele", "Second_allele", "Zygosity")
  progress(r, nrrow)
  if ( r == nrrow) message("Filtering step 1 done...but it's not done yet! Next up: step 2 ")
  Sys.sleep(0.01)
  return(df)
}
