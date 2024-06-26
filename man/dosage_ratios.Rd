% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dosage_ratios.R
\name{dosage_ratios}
\alias{dosage_ratios}
\title{Calculate the Percentage of Each Dosage Value}
\usage{
dosage_ratios(data, ploidy)
}
\arguments{
\item{data}{Genotype matrix or data.frame}

\item{ploidy}{The ploidy of the species being analyzed}
}
\value{
A data.frame with percentages of dosage values in the genotype matrix
}
\description{
This function calculates the percentage of each dosage value within a genotype matrix.
It assumes that the samples are the columns, and the genomic markers are in rows. Missing data should
be set as NA, which will then be ignored for the calculations. All samples must have the same ploidy.
}
