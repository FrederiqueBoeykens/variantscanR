% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateR_ensembl.R
\name{annotateR_ensembl}
\alias{annotateR_ensembl}
\title{annotateR_ensembl}
\usage{
annotateR_ensembl(a, b)
}
\arguments{
\item{a}{ensembl BED file containing 'ENSCAFTxxx' nomenclature for genes.}

\item{b}{ensembl 'all fields' file containing the alternate gene names.}
}
\value{
Returns an annotated BED file containing both gene names.
}
\description{
Changes the gene names of an ensembl BED file.
}
\examples{
\donttest{}
}
