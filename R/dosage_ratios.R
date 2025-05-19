#' Calculate the Percentage of Each Dosage Value
#'
#' This function calculates the percentage of each dosage value within a genotype matrix.
#' It assumes that the samples are the columns, and the genomic markers are in rows. Missing data should
#' be set as NA, which will then be ignored for the calculations. All samples must have the same ploidy.
#'
#' @param data Genotype matrix or data.frame
#' @param ploidy The ploidy of the species being analyzed
#' @return A data.frame with percentages of dosage values in the genotype matrix
#' @examples
#' # example numeric genotype matrix for a tetraploid
#' n_ind <- 5
#' n_snps <- 10
#'
#' geno <- matrix(as.numeric(sample(0:4, n_ind * n_snps, replace = TRUE)), nrow = n_snps, ncol = n_ind)
#' colnames(geno) <- paste0("Ind", 1:n_ind)
#' rownames(geno) <- paste0("SNP", 1:n_snps)
#' ploidy <- 4
#'
#' # ratio of dosage value (numeric genotypes) across samples in dataset
#' result <- dosage_ratios(geno, ploidy)
#'
#' print(result)
#'
#' @export
dosage_ratios <- function(data, ploidy) {
  percentages <- apply(data, 2, function(col) {
    counts <- table(col)
    prop <- prop.table(counts) * 100
    #max_val <- max(as.numeric(names(counts)))  # Find the maximum value in the column
    prop[as.character(0:ploidy)]  # Adjust the range based on the max value (consider entering the ploidy value explicitly for max_val)
  })

  # Calculate percentages for genotype matrices
  percentages_df <- as.data.frame(t(percentages))
  expected_colnames <- as.character(0:ploidy)
  if(ncol(percentages_df) == length(expected_colnames)) { # Basic check
    colnames(percentages_df) <- expected_colnames
  }
  percentages_df$Data <- "Dosages"
  melted_data <- percentages_df %>%
    pivot_longer(cols = -(Data),names_to = "Dosage", values_to = "Percentage")

  return(melted_data)
}

