#' File.csv containing the variants of interest.
#'
#' File containing the variants of interest for the dog.
#' This file was downloaded from the OMIA.org website and manually curated.
#' (CanFam3.1, Sep 2011)
#'
#'
#' @docType data
#'
#' @usage data(omia_raw)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{OMIA.Variant.ID}{OMIA ID given to the variant}
#'  \item{OMIA.Phene.Species.ID.s.}{OMIA phenotype species ID}
#'  \item{Species.Name}{Name of species in which the variant was described}
#'  \item{Breed.s.}{Breeds in which the variant has been found}
#'  \item{Variant.Phenotype}{Variant phenotype}
#'  \item{Gene}{Gene in which the variant is located}
#'  \item{Allele}{Allele}
#'  \item{Type.of.Variant}{Type of variant}
#'  \item{Source.of.Genetic.Variant}{Source of the genetic variant}
#'  \item{Deleterious.}{Is the individual's susceptibility or predisposition incread to a certain disease or disorder?}
#'  \item{Reference.Sequence}{Reference sequence chromosome or scaffold}
#'  \item{Chr.}{The name of chromosome}
#'  \item{g..or.m.}{Genomic or mitochondrial location}
#'  \item{c..or.n.}{Coding or non-codion region}
#'  \item{p.}{Protein}
#'  \item{Verbal.Description}{Verbal description}
#'  \item{EVA.ID}{European Variation Archive ID}
#'  \item{Inferred.EVA.ID}{Inferred EVA ID}
#'  \item{Year.Published}{Year in which the variant was first published}
#'  \item{PubMed.ID.s.}{PubMed ID(s) of publication(s) referring to the variant}
#'  \item{Acknowledgements}{Acknowledgements}
#' }
#' @references This data set was downloaded from <https://omia.org/home/> and manually curated.
#' @keywords datasets
#' @examples
#'
#' data(omia_raw)
#' head(omia_raw)
"omia_raw"
