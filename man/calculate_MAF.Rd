% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_MAF.R
\name{calculate_MAF}
\alias{calculate_MAF}
\title{Calculate Minor Allele Frequency from a Genotype Matrix}
\usage{
calculate_MAF(df, ploidy)
}
\arguments{
\item{df}{Genotype matrix or data.frame}

\item{ploidy}{The ploidy of the species being analyzed}
}
\value{
A dataframe of AF and MAF values for each marker
}
\description{
This function calculates the allele frequency and minor allele frequency from a genotype matrix.
It assumes that the Samples are the columns, and the genomic markers are in rows. Missing data should
be set as NA, which will then be ignored for the calculations. All samples must have the same ploidy.
}
\examples{
# example input for a diploid
geno <- data.frame(
            Sample1 = c(0, 1, 2, NA, 0),
            Sample2 = c(1, 1, 2, 0, NA),
            Sample3 = c(0, 1, 1, 0, 2),
            Sample4 = c(0, 0, 1, 1, NA)
           )
row.names(geno) <- c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")

ploidy <- 2

# calculate allele frequency
result <- calculate_MAF(geno, ploidy)

print(result)

}
