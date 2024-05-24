#' Switch Dosage Values from a Genotype Matrix
#'
#' This function converts the dosage count values to the opposite value. This is primarily used when converting
#' dosage values from reference based (0 = homozygous reference) to alternate count based (0 = homozygous alternate).
#' It assumes that the Samples are the columns, and the genomic markers are in rows. Missing data should
#' be set as NA, which will then be ignored for the calculations. All samples must have the same ploidy.
#'
#' @param df Genotype matrix or data.frame
#' @param ploidy The ploidy of the species being analyzed
#' @return A genotype matrix
#' @export
flip_dosage <- function(df, ploidy, is_reference = TRUE) {
  if (is_reference) {
    # Convert from reference to alternate alleles
    return(abs(df - ploidy))
  } else {
    # Data already represents alternate alleles
    return(df)
  }
}
