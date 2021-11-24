#' @title pRocess
#' @description Optional function that processes the user-defined file with the variants of interest. If desired, this function performes some quality control steps. This quality control might include providing a reference genome assembly (refseq argument, e.g. CanFam3.1) and removes all variants from other assemblies.
#'
#' @param variants_file File containing variants of interest. Important: Specific format needed, check vignette for more information
#' @param BED_file Annotated BED file processed by annotateR_NCBI() or annotateR_ensembl() function
#' @param refseq Reference genome assembly (OPTIONAL) e.g. Canfam3.1
#' @param QC Quality Control (OPTIONAL). If TRUE (default) the function checks whether the nucleotide in the variants_file is the same as that of the reference genome
#' @param organism Optional argument but mandatory if QC = TRUE. This should be a standard BSgenome organism code (list with available organisms and codes is provided in the DATA folder)
#' @param OMIA Optional argument which is by default FALSE. If TRUE, the variants_file will be handled differently. OMIA should be TRUE if the organism file is directly downloaded from the OMIA website.
#'
#' @return Returns a processed file containing the variants of interest, ready for downstream analysis.
#' @export
#'
#' @importFrom stats complete.cases
#' @import rebus
#' @importMethodsFrom Biostrings getSeq
#' @importMethodsFrom Biostrings reverseComplement
#' @importClassesFrom  Biostrings DNAString
#'
#' @examples
#' \dontrun{}


pRocess <- function(variants_file, BED_file, refseq, QC = TRUE, organism = NULL, OMIA = FALSE){
  if (OMIA == TRUE){
    variants_file <- variants_file[stats::complete.cases(variants_file[1]), ]
    colnames(variants_file)[1] <- "Chr."
    variants_file <- variants_file[stats::complete.cases(variants_file[2]), ]
    colnames(variants_file)[2] <- "g. or m"
    colnames(variants_file)[3] <- "Reference Sequence"
  } else {
    variants_file <- variants_file[stats::complete.cases(variants_file[1]), ]
    colnames(variants_file)[1] <- "Chromosome"
    variants_file <- variants_file[stats::complete.cases(variants_file[2]), ]
    colnames(variants_file)[2] <- "Start"
    colnames(variants_file)[3] <- "End"
    colnames(variants_file)[4] <- "Reference"
    colnames(variants_file)[5] <- "Reference Sequence"
  }
  for (r in 1:nrow(variants_file)){
    if (OMIA == TRUE){
      id <- variants_file[r,3]
      id <- gsub(" ", "", id) #remove potential spaces between CanFam and 3.1 so we can filter them out in next step
    } else {
      id <- variants_file[r,5]
      id <- gsub(" ", "", id) #remove potential spaces between CanFam and 3.1 so we can filter them out in next step
    }
  }
  if (missing(refseq)){
    variants_file <- variants_file
  } else {
    variants_file <- subset(variants_file, variants_file$'Reference Sequence'== refseq )
  }
  if (OMIA == TRUE){
    between <- variants_file
    for (r in 1:nrow(between)){
      original <- between[r,2]
      original <- gsub(",","", original)
      original <- gsub(">", "", original)
      original <- gsub("_", "-", original)
      between[r,2] <- original
    }
    final <- data.frame()
    for (r in 1:nrow(between)){
      output = NULL
      row <- between[r,]
      nrc <- ncol(between)
      id <- between[r,2]
      id <- stringr::str_split(id, pattern = ("-"), n = 2, simplify = TRUE)
      start_pos <- id[1,1]
      end_pos <- id[1,2]
      output <- cbind(row[,1], start_pos, end_pos, row[,3:nrc])
      output_temp <- data.frame(output)
      final <- rbind(final, output_temp)
    }
    colnames(final)[1:4] <- c("Chromosome", "Start", "End", "Reference Sequence")
    final <- data.frame(lapply(final, as.character), stringsAsFactors = FALSE)
    for(r in 1:nrow(final)){
      pattern_loc_start <- rebus::optional(or("g","m")) %R% rebus::optional(DOT) %R% rebus::capture(one_or_more(DGT)) %R% rebus::optional(WRD)
      pattern_loc_end <- rebus::capture(one_or_more(DGT)) %R% rebus::optional(WRD) %R% rebus::optional(one_or_more(DGT))
      start <- stringr::str_match(final[r,2], pattern = pattern_loc_start)
      end <- stringr::str_match(final[r,3], pattern = pattern_loc_end)
      final[r,2] <- start[,2]
      final[r,3] <- end[,2]
    }
    for (r in 1:nrow(variants_file)){
      original <- as.character(variants_file[r,2])
      original <- gsub(",","", original)
      original <- gsub("_", "", original)
      original <- gsub("g.", "", original)
      original <- gsub("m.", "", original)
      original <- gsub("[0-9]", "", original)
      original <- gsub("del", "", original)
      original <- gsub("dup", "", original)
      original <- gsub("inv", "", original)
      original <- gsub("ins", "", original)
      if (grepl(">", original, fixed = TRUE)){
        pattern <- rebus::capture(or("T", "A", "C", "G")) %R% ">" %R% or("T", "A", "C", "G")
        nt <- stringr::str_match(original, pattern = pattern)
        original <- nt[1,2]
      } else {
        original <- original
      }
      variants_file[r,2] <- original
    }
    final$Chromosome <- paste0("chr", final$Chromosome)
    nrccc <- ncol(final)
    final <- cbind(final[,1:3], variants_file[,2], final[,4:nrccc])
    colnames(final)[4] <- "Reference"
  } else {
    final <- variants_file
    final$Chromosome <- paste0("chr", final$Chromosome)
  }
  if (QC == TRUE){ #until now everything was the same for QC = TRUE or FALSE
    random1 <- data.frame()
    for (r in 1:nrow(final)){
      row <- final[r,]
      nrcc <- ncol(final)
      a <- final[r,1]
      b <- as.numeric(final[r,2])
      c <- as.numeric(final[r,3])
      c[is.na(c)] <- b
      Check <- Biostrings::getSeq(organism, a, start = b, end = c)
      Check <- as.character(Check)
      random <- cbind(row[,1:4], Check, row[,5:nrcc], row.names = NULL)
      random2 <- data.frame(random)
      random1 <- rbind(random1, random2)
    }
    test <- unique(random1)
    proc_variants_file <- data.frame()
    for (r in 1:nrow(test)){
      end <- ncol(test)
      row <- test[r,]
      gene <- test[r,8]
      index <- which(BED_file$V4 == gene) #gene name
      orientation <- BED_file[index,6]    #orientatio + or -
      orientation <- as.character(orientation)
      if (length(unique(orientation)) == 1){
        orientation <- orientation[1]
      } else {
        orientation <- "Information lacking"
      }
      proc_variants_file1 <- cbind(row[,1:5], orientation, row[,6:end], row.names = NULL)
      proc_variants_file2 <- data.frame(proc_variants_file1)
      proc_variants_file <- rbind(proc_variants_file, proc_variants_file2)
    }
    proc_variants_file <- data.frame(lapply(proc_variants_file, as.character), stringsAsFactors = FALSE)
    for (r in 1:nrow(proc_variants_file)){
      if (proc_variants_file[r,6] == "-") {
        dna <- proc_variants_file[r,5]
        dna <- Biostrings::DNAString(dna)
        dna <- Biostrings::reverseComplement(dna)
        dna <- as.character(dna)
        proc_variants_file[r,5] <- dna
      } else if (proc_variants_file[r,6] == "+"){
        dna <- proc_variants_file[r,5]
        proc_variants_file[r,5] <- dna
      } else {
        dna <- "Check manually"
        proc_variants_file[r,5] <- dna
      }
    }
    proc_variants_file[is.na(proc_variants_file)] <- "NA"
    for (r in 1:nrow(proc_variants_file)){
      if(proc_variants_file[r,4] == proc_variants_file[r,5]){
        proc_variants_file[r,5] <- proc_variants_file[r,5]
      } else {
        proc_variants_file[r,5] <- "Check manually"
      }
    }
  } else {proc_variants_file <- final}
  return(proc_variants_file)
}
