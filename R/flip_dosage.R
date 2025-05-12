#' Switch Dosage Values from a Genotype Matrix
#'
#' This function converts the dosage count values to the opposite value. This is primarily used when converting
#' dosage values from reference based (0 = homozygous reference) to alternate count based (0 = homozygous alternate).
#' It assumes that the Samples are the columns, and the genomic markers are in rows. Missing data should
#' be set as NA, which will then be ignored for the calculations. All samples must have the same ploidy.
#'
#' @param df Genotype matrix or data.frame
#' @param ploidy The ploidy of the species being analyzed
#' @param is.reference The dosage calls value is based on the count of reference alleles (TRUE/FALSE)
#' @return A genotype matrix
#' @examples
#' # example code
#'
#' # example numeric genotype matrix for a tetraploid
#' n_ind <- 5
#' n_snps <- 10
#'
#' geno <- matrix(as.numeric(sample(0:4, n_ind * n_snps, replace = TRUE)), nrow = n_snps, ncol = n_ind)
#' colnames(geno) <- paste0("Ind", 1:n_ind)
#' rownames(geno) <- paste0("SNP", 1:n_snps)
#' ploidy <- 4
#'
#' # Output matrix with the allele count reversed
#' results <- flip_dosage(geno, ploidy, is.reference = TRUE)
#'
#' print(results)
#'
#' @export
flip_dosage <- function(df, ploidy, is.reference = TRUE) {
  if (is.reference) {
    # Convert from reference to alternate alleles
    return(abs(df - ploidy))
  } else {
    # Data already represents alternate alleles
    return(df)
  }
}
