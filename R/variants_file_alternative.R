#' File containing the variants of interest in an alternative format.
#'
#' File containing the variants of interest for the dog.
#' This file was manually created for demonstrative purposes.
#'
#'
#'
#' @docType data
#'
#' @usage data(variants_file_alternative)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{Chromosome}{The name of chromosome}
#'  \item{Start}{Position start of the variant}
#'  \item{End}{Position end of the varians (in case of insertion/deletion)}
#'  \item{Reference}{The name of chromosome}
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
#' @references This data set was manually created for demonstrative purposes for the variantscanR package.
#' @keywords datasets
#' @examples
#'
#' data(variants_file_alternative)
#' head(variants_file_alternative)
"variants_file_alternative"
