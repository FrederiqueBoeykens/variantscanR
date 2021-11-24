#' @title diveRsity
#' @description Provides a measure of diversity based on the percentage of heterozygosity per sample.
#'
#' @param VCF Multisample VCF file
#' @param breeds An Excel file containing the breed of every sample.
#' @param sample_name A string containing the name of the sample of interest that needs to be analysed.
#'
#' @return A dataframe containing the level of heterozygosity of every sample, together with a graph, representing the same information.
#' @export
#'
#' @import vcfR
#' @import dplyr
#' @importFrom  ggplot2 ggplot
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom R3port html_plot html_list html_combine
#'
#'
#' @examples
#' \dontrun{}

diveRsity <- function(VCF, breeds, sample_name){
  ID <- Heterozygosity <- Population <- NULL
  vcf <- vcfR::read.vcfR(VCF)
  example <- vcfR::gt.to.popsum(vcf) # population parameters, also n, we need all the rows where n = max!
  vcf_gt <- vcf@gt # becomes a matrix with only the gt info of the
  population_size <- as.numeric(ncol(vcf_gt)) -1 #the format column is included so: ncol - 1 = n -> this n value is the value we want to subset our example df on
  subsetvcf <- example[example$n == population_size,]
  rowindices <- rownames(subsetvcf) #characterstring chr [1:743614] "30" "31" "46" "48" ... these rowindices are the rows where n = 27
  rowindices <- as.vector(rowindices)
  gtinfosubset <- vcfR::extract.gt(vcf)
  gtinfosubset <- data.frame(gtinfosubset)
  row.names(gtinfosubset) <- 1:nrow(gtinfosubset)
  newgtinfosubset <- gtinfosubset[rownames(gtinfosubset) %in% rowindices,]
  newgtinfosubset <- data.frame(lapply(newgtinfosubset, as.character), stringsAsFactors = FALSE)
  try <- apply(newgtinfosubset, c(1,2), function(x) if (x == "1/1" | x == "0/0"){x <- 0} else {x <- 1})
  try2 <- as.data.frame(colSums(try != 0)) #you get a list of named, numeric columns of unequal size
  colnames(try2) <- "Heterozygosity"
  nrrows <- nrow(try)
  pop <- breeds$Breed #extract pop vector from file that contains breeds of patients
  try2$`Heterozygosity`<- as.numeric(try2$`Heterozygosity`) / nrrows
  try2$Population <- pop
  try2$Population <- as.factor(try2$Population)
  try2 <- cbind(ID = rownames(try2), try2)
  rownames(try2) <- 1:nrow(try2)
  path <- getwd()
  add_folder <- "/Diversity_"
  sample_folder <- paste(add_folder, sample_name, sep = "")
  new_path <- paste(path, sample_folder, sep = "")
  dir.create(new_path) #created a new folder with the sample as its name
  #now go into that folder to store the created raw html files.
  setwd(new_path)
  p <- ggplot2::ggplot(try2, ggplot2::aes(y=Heterozygosity, x=Population)) + ggplot2::geom_dotplot(binaxis = 'y', stackdir = 'center')
  q <- ggplot2::ggplot(try2, ggplot2::aes(y=Heterozygosity, x=Population, label=ID)) + ggplot2::geom_dotplot(binaxis = 'y', stackdir = 'center') + ggrepel::geom_label_repel(size=3, box.padding = ggplot2::unit(0.5, "lines")) + ggplot2::theme_classic()
  subset_max <- try2 %>% group_by(Population) %>% filter(Heterozygosity == max(Heterozygosity))
  subset_max <- data.frame(subset_max)
  subset_max$ID <- as.character(subset_max$ID)
  try3 <- try2
  try3$ID <- as.character(try3$ID)
  try3$name <- ""
  for(r in 1:nrow(subset_max)){
    name <- subset_max[r,1]
    index <- which(try3$ID == name)
    sample <- try3[index,1]
    try3[index,4] <- sample
  }
  pq <- ggplot2::ggplot(try3, ggplot2::aes(x=Population, y=Heterozygosity, label=name)) + ggrepel::geom_text_repel() + ggplot2::geom_point(color = ifelse(try3$name == "", "grey50", "orange"))
  try4 <- try2
  try4$ID <- as.character(try4$ID)
  try4$name <- ""
  for (r in 1:nrow(try4)){
    index <- which(try4$ID == sample_name)
    sample <- try4[index,1]
    try4[index,4] <- sample
  }
  pqr <- ggplot2::ggplot(try4, ggplot2::aes(x=Population, y=Heterozygosity, label=name)) + ggrepel::geom_text_repel() + ggplot2::geom_point(color = ifelse(try4$name == "", "grey50", "orange"))
  html <- "_diversity.html"
  sample_plot <- "Diversity:"
  html_full <- paste(sample_name, html, sep = "")
  sample_plot_full <- paste(sample_plot, sample_name, sep = " ")
  R3port::html_plot(pqr,out= html_full,
            title= sample_plot_full)
  R3port::html_plot(pq,out="highlight_plot.html",
            title="Diversity: highlights")
  R3port::html_plot(p,out="basic_plot.html",
            title="Diversity: basic overview")
  R3port::html_plot(q,out="annotated_plot.html",
            title="Diversity: fully annotated")
  R3port::html_list(try2,title="Diversity", out="Diversity_table.html")
  R3port::html_combine(out="report_variantscanR.html",toctheme=TRUE,
               template=paste0(system.file(package="R3port"),"/bootstrap.html"))
  #At the end, so here, we have to go back to the original folder. Otherwise we would get folder inception
  setwd(path)
  return(list(try2,p,q,pq,pqr))
}
