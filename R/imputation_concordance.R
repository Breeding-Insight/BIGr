#' Calculate Concordance between Imputed and Reference Genotypes
#'
#' This calculates the concordance between imputed and reference genotypes. It assumes that samples are rows and markers are columns.
#' It is recommended to use allele dosages (0,1,2) but will work with other formats. Missing data in reference or imputed genotypes
#' will not be considered for concordance if argument missing_code used. If a specific subset of markers should it can be provided as argument snps_2_exclude.
#'
#' @param reference_genos Genotype data.frame with rows as samples and columns as markers. Dosage recommended.
#' @param imputed_genos Genotype data.frame with rows as samples and columns as markers. Dosage recommended.
#' @param missing_code Optional input to consider missing data to exclude in concordance calculation.
#' @param snps_2_exclude Optional input to exclude specific markers from concordance calculation. Single column of marker ids.
#' @param output_df_name Optional input to assign the output dataframe to a specific variable name. Default is "imputation_concordance"
#' @import dplyr
#' @return 2 outputs: 1) A data frame with sample IDs and concordance percentages. 2) A summary of concordance percentages.
#' @export
#'
imputation_concordance <- function(reference_genos, imputed_genos, missing_code = NULL, snps_2_exclude = NULL, output_df_name = "imputation_concordance") {

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

  # Assign the result dataframe to the output variable
  assign(output_df_name, result_df, envir = .GlobalEnv)

  # Print mean concordance
  summary_concordance <- summary(percentage_match, na.rm = TRUE) * 100
  names(summary_concordance) <- c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")

  cat("Concordance Summary:\n")
  for (name in names(summary_concordance)) {
    cat(name, ":", round(summary_concordance[name], 2), "%\n")
  }
}

