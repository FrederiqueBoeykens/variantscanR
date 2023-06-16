#' @title variantfiltR
#' @description Filtering function that takes the single sample VCF file and filters for present variants of
#' interest.
#'
#' @param VCF Single sample VCF file, uploaded with the vcfscanneR function.
#' @param variants_file File containing variants of interest, processed by the pRocess function.
#' @param BED_file_annot Annotated BED file created with the annotateR_NCBI/ensembl function.
#' @param breed Breed of interest (Capital letter insensitive).
#'
#' @return A html report will be obtained containing an interactive overview of the various tables acquired
#' after filtering. Next to this, breeding advice on the zygosity state of the variant is provided in the report.
#' @export
#'
#' @importFrom svMisc progress
#' @importFrom R3port html_list html_combine
#' @importFrom stringr fixed str_replace_all
#'
#' @examples
#' \dontrun{}


variantfiltR <- function(vcf, variants_file, BED_file_annot, breed){
  cicip <- colnames(variants_file)[5]
  if (cicip == "Check"){
    variants_file <- variants_file[,-c(5:6)]
  } else {variants_file <- variants_file}
  filter <- data.frame()
  vcf$pos <- as.numeric(vcf$pos)
  variants_file$Start <- as.numeric(variants_file$Start)
  for (r in 1:nrow(variants_file)){ #for every row in the variants_file
    nrrow <- nrow(variants_file)
    x <- variants_file[r,1] #get the chromosome we're looking at  str(chr)
    y <- as.numeric(variants_file[r,2]) #get the location whithin the chromosome  str(num)
    subsetvcffile <- vcf[vcf$chrom == x,] #subset on only this chromosome
    subsetvcffile <- subsetvcffile[subsetvcffile$pos == y,]  #next use location within chromosome to subset the subset
    filter2 <- data.frame(subsetvcffile)
    filter <- rbind.data.frame(filter, filter2)
    svMisc::progress(r, nrrow)
    if ( r == nrrow) message("Filtering step 1 done...but it's not done yet! Next up: step 2 ")
    Sys.sleep(0.01)
  }
  filter <- unique(filter)
  if (nrow(filter) == 0) { #changed after the reviewing process
    stop("It appears that none of the variants of interest are present in the provided VCF file.")
  }
  output1 <- data.frame()
  filter <- filter[,-3]
  filter <- filter[,-5:-7]
  samplename <- names(filter)[6]
  for (r in 1:nrow(filter)){
    nrrow <- nrow(filter)
    row <- filter[r,]
    format <- unlist(strsplit(as.character(row[,5]), ":"))
    ex <- which(format == "GT")  #depending on where GT is in the format string,
    snp <- unlist(strsplit(as.character(row[,6]), ":")) #snp in sample_file column
    snp <- snp[ex] #we choose the same number where GT was located
    snp <- gsub("\\|", "/", snp)
    snp_sep <- unlist(strsplit(snp, "/"))
    snp1 <- snp_sep[1] #homozyg 0/0 or 1/1 and heterozyg 0/1 or 1/0 so therefore as.numeric because
    #until now chr
    snp2 <- snp_sep[2]
    if (snp1 == "." ){snp1 <- "A call cannot be made for this sample at this given locus"} else {
      snp1 <- as.numeric(snp1)
    }
    if (snp2 == "." ){snp2 <- "A call cannot be made for this sample at this given locus"} else {
      snp2 <- as.numeric(snp2)
    }
    #if numerical, there are stored as numerical values for downstream analysis
    zygosity <- if(snp1 == "A call cannot be made for this sample at this given locus" | snp2 == "A call cannot be made for this sample at this given locus" ) {"Zygosity could not be determined"} else if (snp1 == snp2){
      "Homozygous"} else {"Heterozygous"}
    allele_alt <- unlist(strsplit(as.character(row[,4]), ","))
    if (is.numeric(snp1) == FALSE){allele1 <- snp1} else if(snp1 == 0){allele1 <- row[,3]} else {allele1 <- allele_alt[snp1]} #location in row of alternative
    #alleles for example alt: A,T,C and a 2 is given than 2 means T because 0 is ref_allele
    #and starting from 1 to.... give alternative allles, so we find that location
    if (is.numeric(snp2) == FALSE){allele2 <- snp2} else if(snp2 == 0){allele2 <- row[,3]} else {allele2 <- allele_alt[snp2]}
    output <- cbind(row[,1:3], allele1, allele2, zygosity)
    colnames(output) <- c("Chromosome", "Location", "Wild type", "First allele", "Second allele", "Zygosity")
    output2 <- data.frame(output)
    output1 <- rbind(output1, output2)
    svMisc::progress(r, nrrow)
    if ( r == nrrow) message("Filtering step2 done! Next up: annotation part1")
    Sys.sleep(0.01)
  }
  result <- unique(output1)
  result[] <- lapply(result, as.character)
  result$Location <- as.numeric(result$Location)
  #annotate filtered vcf file with the variants_file
  annot <- data.frame()
  for (r in 1:nrow(result)){ #for every row in the variants_file
    nrrow <- nrow(result)
    nccol <- ncol(variants_file)
    row <- result[r,]
    x <- result[r,1] #get the chromosome we're looking at
    y <- result[r,2]
    subset <- variants_file[variants_file$Chromosome == x,] #subset on only this chromosome
    subset <- subset[subset$Start == y,]  #next use location within chromosome
    annot1 <- cbind(row, subset[,5:nccol])
    annot2 <- data.frame(annot1)
    annot <- rbind(annot, annot2)
    svMisc::progress(r, nrrow)
    if ( r == nrrow) message("Annotation part1 done! Next up: annotation part2 ")
    Sys.sleep(0.01)
  }
  # now we want the essentials to report, which are
  # chr - location - exon - Refseq transcript - variant - nucleotide - variant protein?  - overerving - breeding advice
  rownames(annot) <- c(1:nrow(annot))
  #keep <- apply(annot[3:5], 1, function(x) length(unique(x[!is.na(x)])) != 1)
  #keept <- annot[keep, ]
  #kept <- rownames(keept)
  colnames(BED_file_annot) <- c("chrom", "start", "end", "gene", "score", "strand", "thickstart", "thickend", "rgb", "blockcount", "blocksizes", "blockstarts", "transcript")
  BED_file_annot$start <- as.numeric(BED_file_annot$start)
  BED_file_annot$end <- as.numeric(BED_file_annot$end)
  BED_file_annot$blockcount <- as.numeric(BED_file_annot$blockcount)
  annot2 <- data.frame()
  for (r in 1:nrow(annot)){ # last column is already the transcript that was originally in the BED file :p
    nnrow <- nrow(annot)
    row <- annot[r,]
    nncol <- ncol(BED_file_annot)
    nnncol <- ncol(annot)
    x <- annot[r,1] #chromosome
    y <- annot[r,2] #location
    gene <- annot[r,9] #gene
    subset <- BED_file_annot[BED_file_annot$gene == gene,] #we don't use location but gene
    #so when using duplicates it duplicates the amount of times PLUS the original one, so we have to do nrow(subset)-1 but what
    #can we do for the subsets with only 1 row, if you do -1 is becomes 0 and we lose the row and this can not happen.
    #we check the amount of rows
    #if is is more than 1, we can just use the duplicates manner
    #otherwise, duplicates becomes the row itself and doesn't need duplicating
    nr <- nrow(subset)
    if(nr > 1){
      nr <- nr-1
      duplicates <- rbind(row, row[rep(1,nr),])
    } else(duplicates <- row)
    BED <- cbind(duplicates[,1:2], duplicates[,9], duplicates[,3:6], duplicates[,8], duplicates[,10:11], subset[,5:nncol], duplicates[,11:nnncol])
    colnames(BED)[3] <- "Gene"
    colnames(BED)[8] <- "Inheritance.pattern"
    annot3 <- data.frame(BED)
    annot2 <- rbind(annot2, annot3)
    svMisc::progress(r, nnrow)
    if ( r == nnrow) message("Annotation part2 done! Next up: reporting ")
    Sys.sleep(0.01)
  }
  # We obtained a bigger dataframe now because every gene can have multiple transcripts
  # Hier heb je terug een loop nodig aangezien er meerdere (transcripts, niet genen) in de subset kunnen zitten.
  # because BED file (0-based) vs VCF (1-based) coordinate systems, we need to adjust one of the two.
  # the location is obtained from OMIA-file/VCF and thus, 1-based
  # I prefer adjusting the information of the BED file (so +1)
  annot2$Gene <- as.character(annot2$Gene)
  annot2$Location <- as.numeric(annot2$Location)  # stays the same
  annot2$thickstart <- as.numeric(annot2$thickstart)
  annot2$thickstart <- annot2$thickstart +1
  annot2$thickend <- as.numeric(annot2$thickend)
  annot2$thickend <- annot2$thickend +1
  annot2$blockcount <- as.numeric(annot2$blockcount) # the amount of exons/blocks stays the same
  report <- data.frame()
  rownames(annot2) <- c(1:nrow(annot2))
  for (r in 1:nrow(annot2)){
    row <- annot2[r,]
    nbrcol <- ncol(annot2)
    nrrow <- nrow(annot2)
    a <- annot2[r,13] #chromStart - thickstart (already +1)
    b <- annot2[r,14] #chromEnd - thickend (already +1)
    exons <- annot2[r,16] #the amount of exons
    exonstarts <- unlist(strsplit(as.character(annot2[r,18]), ","))
    exonstarts <- as.numeric(exonstarts)
    exonstarts <- exonstarts # DO NOT ADD 1 - otherwise in the end, you have added 2, the 'starts' depend on the original value 'a'
    exonstarts <- vapply(exonstarts, "+", a)
    exonlengths <- unlist(strsplit(as.character(annot2[r,17]), ","))
    exonlengths <- as.numeric(exonlengths)
    exonlengths <- exonstarts + exonlengths
    dfpos <- cbind(exonstarts, exonlengths)
    dfpos <- data.frame(dfpos)
    colnames(dfpos)[1] <- "start"
    colnames(dfpos)[2] <- "length"
    y <- annot2[r,2]
    posit <- with(dfpos, dfpos$start <= y & dfpos$length >= y)
    exonnr <- which(posit == TRUE)
    if (length(exonnr) == 0){
      exon_or_intron <- "Intronic"
    } else {
      exonnr <- as.character(exonnr)
      exon_or_intron <- "exon"
      exon_or_intron <- paste(exon_or_intron,exonnr, sep = "")
    }
    report1 <- cbind(row[,1:3], exon_or_intron, row[,4:nbrcol])
    report2 <- tibble(report1)
    report <- rbind(report,report2)
    svMisc::progress(r, nrrow)
    if ( r == nrrow) message("Annotation part2 done! Next up: reporting ")
    Sys.sleep(0.01)
  }
  finalreport <- data.frame()
  for (r in 1:nrow(report)){
    nrrow <- nrow(report)
    row <- report[r,]
    Chromosome <- row[,1]
    Location <- row[,2]
    Gene <- row[,3]
    Exon <- row[,4]
    WildType <- row[,5]
    Nucleotide1 <- row[,6]
    Nucleotide2 <- row[,7]
    zygosity <- row[,8]
    RefSeq_transcipt <- row[,20]
    Inheritance_pattern <- row[,9]
    breeds <- row[,11]
    Phenotype <- row[,10]
    finalreport1 <- cbind(Chromosome,Location,Gene,Exon,WildType,Nucleotide1,Nucleotide2,zygosity,RefSeq_transcipt,Inheritance_pattern,Phenotype,breeds)
    finalreport2 <- data.frame(finalreport1)
    finalreport <- rbind(finalreport, finalreport2)
    svMisc::progress(r, nrrow)
    if ( r == nrrow) message("Waiting for report")
    Sys.sleep(0.01)
  }
  colnames(finalreport) <- c("Chromosome", "Location", "Gene", "Exon or Intron", "Wild Type", "Allele 1", "Allele 2", "Zygosity", "Refseq Transcript", "Inheritance Pattern", "Variant Phenotype", "Breed(s)")
  keep <- apply(finalreport[5:7], 1, function(x) length(unique(x[!is.na(x)])) != 1)
  keept <- finalreport[keep, ]
  kept <- as.numeric(rownames(keept))
  high <- finalreport[kept,] #the more important rows
  finalreport <- finalreport[-kept,] #What if there are no homozygous WT? exceptional case, but still
  if (nrow(high) == 0){
    high <- high
  } else{
    rownames(high) <- c(1:nrow(high)) #Then this step gives problems
  }
  if (nrow(finalreport) == 0){
    finalreport <- finalreport
  } else{
    rownames(finalreport) <- c(1:nrow(finalreport))
  }
  #kept contains the 'name' of the rows that are homozygous variant OR heterozygous, so the more important rows that need to be reported with higher importance
  prefix <- "Sample:"
  sample_result <- paste(prefix, samplename)
  #breed argument
  if (missing(breed)){
    high <- data.frame(lapply(high, as.character), stringsAsFactors=FALSE)
    finalreport <- data.frame(lapply(finalreport, as.character), stringsAsFactors=FALSE)
    coln <- colnames(finalreport)
    title0 <- c("Variants", "present", "in", "sample", ":", " ", " ", " ", " ", " ", " ", " " )
    extrarow1 <- as.data.frame(matrix(title0, ncol = 12, byrow = T))
    colnames(extrarow1) <- coln
    title2 <- c("Variants", "of", "interest", "found", "in", "sample", "but", "homozygous", "wild", "type", ": ", "" )
    extrarow2 <- as.data.frame(matrix(title2, ncol = 12, byrow = T))
    colnames(extrarow2) <- coln
    empty <- c(" "," "," "," "," "," "," "," "," "," "," "," ")
    emptyrow <- as.data.frame(matrix(empty, ncol = 12, byrow = T))
    colnames(emptyrow) <- coln
    R3port::html_list(high,title="Variants of high importance found within sample",
                      footnote= sample_result,out="variants_maintaineda.html")
    R3port::html_list(finalreport,title="Variants present in sample but homozygous Wild Type",
                      footnote= sample_result,out="variants_maintainedb.html")
    df1 <- rbind(title0,emptyrow,high,emptyrow,extrarow2,emptyrow,finalreport)
  } else {
    vectori <- numeric(0)
    for (r in 1:nrow(high)){
      breed_low <- tolower(breed)
      breed_nospace <- stringr::str_replace_all(breed_low, stringr::fixed(" "), "")
      row <- high[r,]
      breed_compare <- tolower(as.character(row[12]))
      breeds_nospace <- str_replace_all(breed_compare, fixed(" "), "")
      breeds_split <- stringr::str_split(breeds_nospace, ",")
      compare <- grepl(breed_nospace, breeds_split)
      #if compare is TRUE, remember rowname and add it to vector
      #at the end make two separate dataframes
      #do this for high and finalreport dataframe
      #at the end end, paste 4 dataframes to each other with in between an empty row containing a title
      #if argument breed was not given, just ad high and finalreport together
      # OR finalreport just is a separate raw.html file?
      if (compare == TRUE){
        new_row <- r
        vectori <- c(vectori, new_row)
      } else{
        high <- high
      }
    }
    vectori <- as.numeric(vectori)
    highhigh <- high[vectori,]
    if (length(vectori) == 0){
      highhigh <- highhigh
      high <- high
    } else {
      highhigh <- high[vectori,]
      high <- high[-vectori,]
      rownames(highhigh) <- NULL #changed after reviewing process
      rownames(high) <- NULL #changed after reviewin process.
    }
    R3port::html_list(highhigh,title="Variants of high importance found within breed of interest",
                      footnote= sample_result,out="variants_maintained1.html")
    R3port::html_list(high,title="Variants of high importance but not found within breed of interest",
                      footnote= sample_result,out="variants_maintained2.html")
    R3port::html_list(finalreport,title="Variants present in sample but homozygous Wild Type",
                      footnote= sample_result,out="variants_maintained3.html")
    #add all 3 dataframes together
    highhigh <- data.frame(lapply(highhigh, as.character), stringsAsFactors=FALSE)
    high <- data.frame(lapply(high, as.character), stringsAsFactors=FALSE)
    finalreport <- data.frame(lapply(finalreport, as.character), stringsAsFactors=FALSE)
    coln <- colnames(finalreport)
    title0 <- c("Variants", "present", "in", "sample", "and", "found", "in", "breed", "of", "interest", "! ", ":" )
    extrarow0 <- as.data.frame(matrix(title0, ncol = 12, byrow = T))
    colnames(extrarow0) <- coln
    title1 <- c("Variants", "present", "in", "sample", "but", "not", "found", "in", "breed", "of ", "interest ", ":" )
    extrarow1 <- as.data.frame(matrix(title1, ncol = 12, byrow = T))
    colnames(extrarow1) <- coln
    title2 <- c("Variants", "of", "interest", "found", "in", "sample", "but", "homozygous", "wild", "type", ": ", "" )
    extrarow2 <- as.data.frame(matrix(title2, ncol = 12, byrow = T))
    colnames(extrarow2) <- coln
    empty <- c(" "," "," "," "," "," "," "," "," "," "," "," ")
    emptyrow <- as.data.frame(matrix(empty, ncol = 12, byrow = T))
    colnames(emptyrow) <- coln
    df1 <- rbind(title0,emptyrow,highhigh,emptyrow,extrarow1,emptyrow,high,emptyrow,extrarow2,emptyrow,finalreport)
  }
  #HOMOZYGOUS DATAFRAME
  homozygousdf <- data.frame()
  homozygousdf[1,1] <- "autosomal recessive"
  homozygousdf[2,1] <- "autosomal dominant"
  homozygousdf[3,1] <- "X linked dominant (M)"
  homozygousdf[4,1] <- "X linked dominant (F)"
  homozygousdf[5,1] <- "X linked recessive (M)"
  homozygousdf[6,1] <- "X linked recessive (F)"
  homozygousdf[7,1] <- "mitochondrial"
  homozygousdf[8,1] <- "Y linked"
  homozygousdf[9,1] <- "NA"
  homozygousdf[1,2] <- ":...ONLY combine with wild type animal! Offspring will be carrier and, on its turn, can ONLY be combined with wild type animal.....................***"
  homozygousdf[2,2] <- ":...Animal can NOT be used for breeding purposes. Offspring would be carrier or homozygous and could also develop symptoms........................."
  homozygousdf[3,2] <- ":...If male: affected animal can NOT be used for breeding purposes as all female offspring will be carrier and thus will be affected....................."
  homozygousdf[4,2] <- ":...If female: Do NOT use animal for breeding purposes because all male offspring will inherit the defective X chromosome..............................."
  homozygousdf[5,2] <- ":...If male: affected animal can be used ONLY if combined with a wild type female animal. Female offspring: carrier.  Male offspring free of disease."
  homozygousdf[6,2] <- ":...If female: Do NOT use animal for breeding purposes as all male offspring will inherit defective X chromosome..........................................."
  homozygousdf[7,2] <- ":...If female: animal can NOT be used for breeding purposes. Mitochondria are inherited maternally. If male: animal can be used......................."
  homozygousdf[8,2] <- ":...If male: animal can NOT be used for breeding purposes as every male offspring will inherit the defective Y chromosome..............................."
  homozygousdf[9,2] <- ":...Not able to provide breeding advice because inheritance pattern was not made available...................................................................."
  colnames(homozygousdf) <- c("Inheritance pattern", "Breeding advice")
  R3port::html_list(homozygousdf,title="Zygosity: Homozygous",
                    footnote="***!IMPORTANT SIDE NOTE!: Do NOT use animal if signs or symptoms are already showing!
          Animal must be capable of carrying out pregnancy.",out="homozygousIP.html")
  #HETEROZYGOUS DATAFRAME
  heterozygousdf <- data.frame()
  heterozygousdf[1,1] <- "autosomal recessive"
  heterozygousdf[2,1] <- "autosomal dominant"
  heterozygousdf[3,1] <- "X linked dominant (M)"
  heterozygousdf[4,1] <- "X linked dominant (F)"
  heterozygousdf[5,1] <- "X linked recessive (M)"
  heterozygousdf[6,1] <- "X linked recessive (F)"
  heterozygousdf[7,1] <- "mitochondrial"
  heterozygousdf[8,1] <- "Y linked"
  heterozygousdf[9,1] <- "NA"
  heterozygousdf[1,2] <- ":...Animal can be used for breeding purposes ONLY if combined with a wild type animal. Offspring will be 50% carrier and 50% wild type.............."
  heterozygousdf[2,2] <- ":...Animal can NOT be used for breeding purposes......................................................................................................................."
  heterozygousdf[3,2] <- ":...If male: affected animal can NOT be used for breeding purposes as all female offspring will be carrier and thus will be affected....................."
  heterozygousdf[4,2] <- ":...If female: Do NOT use animal for breeding purposes as all male offspring will inherit defective X chromosome..........................................."
  heterozygousdf[5,2] <- ":...If male: affected animal can be used ONLY if combined with a wild type female animal. Female offspring: carrier.  Male offspring free of disease."
  heterozygousdf[6,2] <- ":...If female: animal can NOT be used for breeding purposes as 50% of male offspring will inherit defective X chromosome................................"
  heterozygousdf[7,2] <- ":...If female: animal can NOT be used for breeding purposes. Mitochondria are inherited maternally. If male: animal can be used......................."
  heterozygousdf[8,2] <- ":...If male: animal can NOT be used for breeding purposes as every male offspring will inherit the defective Y chromosome..............................."
  heterozygousdf[9,2] <- ":...Not able to provide breeding advice because inheritance pattern was not made available...................................................................."
  colnames(heterozygousdf) <- c("Inheritance pattern", "Breeding advice")
  R3port::html_list(heterozygousdf,title="Zygosity: Heterozygous",out="heterozygousIP.html")
  R3port::html_combine(out="report_variantscanR.html",toctheme=TRUE,
                       template=paste0(system.file(package="R3port"),"/bootstrap.html"))
  return(df1)
}



