% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variantfiltR.R
\name{variantfiltR}
\alias{variantfiltR}
\title{variantfiltR}
\usage{
variantfiltR(VCF, variants_file, BED_file_annot, breed)
}
\arguments{
\item{VCF}{Single sample VCF file, uploaded with the vcfscanneR function.}

\item{variants_file}{File containing variants of interest, processed by the pRocess function.}

\item{BED_file_annot}{Annotated BED file created with the annotateR_NCBI/ensembl function.}

\item{breed}{Breed of interest (Capital letter insensitive).}
}
\value{
A html report will be obtained containing an interactive overview of the various tables acquired
after filtering. Next to this, breeding advice on the zygosity state of the variant is provided in the report.
}
\description{
Filtering function that takes the single sample VCF file and filters for present variants of
interest.
}
\examples{
\dontrun{}
}
