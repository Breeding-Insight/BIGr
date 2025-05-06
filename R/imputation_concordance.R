#' Calculate Concordance between Imputed and Reference Genotypes
#'
#' This function calculates the concordance between imputed and reference genotypes. It assumes that samples are rows and markers are columns.
#' It is recommended to use allele dosages (0, 1, 2) but will work with other formats. Missing data in reference or imputed genotypes
#' will not be considered for concordance if the `missing_code` argument is used. If a specific subset of markers should be excluded,
#' it can be provided using the `snps_2_exclude` argument.
#'
#' @param reference_genos A data frame containing reference genotype data, with rows as samples and columns as markers. Dosage format (0, 1, 2) is recommended.
#' @param imputed_genos A data frame containing imputed genotype data, with rows as samples and columns as markers. Dosage format (0, 1, 2) is recommended.
#' @param missing_code An optional value to specify missing data. If provided, loci with this value in either dataset will be excluded from the concordance calculation.
#' @param snps_2_exclude An optional vector of marker IDs to exclude from the concordance calculation.
#' @param verbose A logical value indicating whether to print a summary of the concordance results. Default is FALSE.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{result_df}: A data frame with sample IDs and their concordance percentages.
#'     \item \code{summary_concordance}: A summary of concordance percentages, including minimum, maximum, mean, and quartiles.
#'   }
#'
#' @details The function identifies common samples and markers between the reference and imputed genotype datasets. It calculates the percentage of matching genotypes for each sample, excluding missing data and specified markers. The concordance is reported as a percentage for each sample, along with a summary of the overall concordance distribution.
#'
#' @import dplyr
#' @export
#'
imputation_concordance <- function(reference_genos,
                                   imputed_genos,
                                   missing_code = NULL,
                                   snps_2_exclude = NULL,
                                   verbose = FALSE) {

  # Find common IDs
  common_ids <- intersect(imputed_genos$ID, reference_genos$ID)

  imputed_genos <- imputed_genos %>% filter(ID %in% common_ids) %>% arrange(ID)
  reference_genos <- reference_genos %>% filter(ID %in% common_ids) %>% arrange(ID)

  # Find common SNPs, excluding those in snps_2_exclude if provided
  common_snps <- setdiff(
    intersect(colnames(imputed_genos), colnames(reference_genos)),
    as.vector(unlist(snps_2_exclude))
  )

  # Subset and convert to matrices for faster computation
  imputed_matrix <- as.matrix(imputed_genos[, common_snps])
  reference_matrix <- as.matrix(reference_genos[, common_snps])

  # Identify valid SNPs that are not missing in either dataset
  if (!is.null(missing_code)) {
    valid_snps <- (imputed_matrix != missing_code) & (reference_matrix != missing_code)
  } else {
    valid_snps <- matrix(TRUE, nrow = nrow(imputed_matrix), ncol = ncol(imputed_matrix))
  }

  # Compute concordance (row-wise percentage of matching SNPs)
  matches <- (imputed_matrix == reference_matrix) & valid_snps
  percentage_match <- rowSums(matches, na.rm = TRUE) / rowSums(valid_snps, na.rm = TRUE)

  # Create output data frame
  result_df <- data.frame(
    ID = common_ids,
    Concordance = paste0(round(percentage_match * 100, 2), "%")
  )

  # Print mean concordance
  summary_concordance <- summary(percentage_match, na.rm = TRUE) * 100
  names(summary_concordance) <- c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")

  if (verbose) {
    message("Concordance Summary:\n")
    for (name in names(summary_concordance)) {
      cat(name, ":", round(summary_concordance[name], 2), "%\n")
    }
  }

  return(result_df)
}

