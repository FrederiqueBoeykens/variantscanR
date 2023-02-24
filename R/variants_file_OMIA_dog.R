#' File containing the variants of interest.
#'
#' File containing the variants of interest for the dog.
#' This file was downloaded from the OMIA.org website and manually curated.
#' (CanFam3.1, Sep 2011)
#'
#'
#' @docType data
#'
#' @usage data(variants_file_OMIA_dog)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{Chr.}{The name of chromosome}
#'  \item{g. or m.}{Genomic or mitochondrial location}
#'  \item{Reference Sequence}{Reference sequence chromosome or scaffold}
#'  \item{Inheritance pattern}{Inheritance pattern}
#'  \item{Gene}{Gene in which the variant is located}
#'  \item{Variant Phenotype}{Variant phenotype}
#'  \item{Breed}{Breeds in which the variant has been found}
#'  \item{OMIA ID(s)}{OMIA variant ID}
#'  \item{Species Name}{Name of species in which the variant was described}
#'  \item{Allele}{Allele}
#'  \item{Type of Variant}{Type of variant}
#'  \item{Deleterious?}{Is the individual's susceptibility or predisposition incread to a certain disease or disorder?}
#'  \item{c. or n.}{Coding or non-codion region}
#'  \item{p.}{Protein}
#'  \item{EVA ID}{European Variation Archive ID}
#'  \item{Year published}{Year in which the variant was first published}
#'  \item{PubMed ID(s)}{PubMed ID(s) of publication(s) referring to the variant}
#' }
#' @references This data set was downloaded from <https://omia.org/home/> and manually curated.
#' @keywords datasets
#' @examples
#'
#' data(variants_file_OMIA_dog)
#' head(variants_file_OMIA_dog)
"variants_file_OMIA_dog"
