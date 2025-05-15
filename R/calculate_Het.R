#' Calculate Observed Heterozygosity from a Genotype Matrix
#'
#' This function calculates the observed heterozygosity from a genotype matrix.
#' It assumes that the samples are the columns, and the genomic markers are in rows. Missing data should
#' be set as NA, which will then be ignored for the calculations. All samples must have the same ploidy.
#'
#' @param geno Genotype matrix or data.frame
#' @param ploidy The ploidy of the species being analyzed
#' @return A dataframe of observed heterozygosity values for each sample
#' @examples
#' # example input for a diploid
#' geno <- data.frame(
#'             Sample1 = c(0, 1, 2, NA, 0),
#'             Sample2 = c(1, 1, 2, 0, NA),
#'             Sample3 = c(0, 1, 1, 0, 2),
#'             Sample4 = c(0, 0, 1, 1, NA)
#'            )
#' row.names(geno) <- c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")
#'
#' ploidy <- 2
#'
#' # calculate observed heterozygosity
#' result <- calculate_Het(geno, ploidy)
#'
#' print(result)
#'
#' @export
calculate_Het <- function(geno, ploidy) {
  # Determine the heterozygous values based on ploidy
  heterozygous_values <- seq(1, ploidy - 1)

  # Create a logical matrix where TRUE represents heterozygous loci
  is_heterozygous <- sapply(geno, function(x) x %in% heterozygous_values)

  # Count the number of heterozygous loci per sample, ignoring NAs
  heterozygosity_counts <- colSums(is_heterozygous, na.rm = TRUE)

  # Calculate the total number of non-NA loci per sample
  total_non_na_loci <- colSums(!is.na(geno))

  # Compute the proportion of heterozygous loci
  heterozygosity_proportion <- heterozygosity_counts / total_non_na_loci

  # Create a dataframe with Sample ID and Observed Heterozygosity
  result_df <- data.frame(
    SampleID = colnames(geno),
    ObservedHeterozygosity = heterozygosity_proportion
  )

  return(result_df)
}
