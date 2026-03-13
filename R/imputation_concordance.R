#' Calculate Concordance between Imputed and Reference Genotypes
#'
#' This function calculates the concordance between imputed and reference
#' genotypes. It assumes that samples are rows and markers are columns.
#' Allele dosages (0, 1, 2) are recommended but other numeric formats are supported.
#' Missing data in either dataset can be excluded from the concordance calculation
#' using the `missing_code` argument. Specific markers can be excluded using
#' the `snps_2_exclude` argument.
#'
#' @param reference_genos A data frame containing reference genotype data,
#' with rows as samples and columns as markers. Must include a column named `ID`.
#'
#' @param imputed_genos A data frame containing imputed genotype data,
#' with rows as samples and columns as markers. Must include a column named `ID`.
#'
#' @param missing_code Optional value specifying missing data. If provided,
#' loci with this value in either dataset will be excluded from the concordance calculation.
#'
#' @param snps_2_exclude Optional vector of marker IDs to exclude from the concordance calculation.
#'
#' @param verbose Logical. If `TRUE`, prints summary statistics (minimum, quartiles,
#' median, mean, maximum) of concordance percentages.
#'
#' @param plot Logical. If `TRUE`, produces a bar plot of concordance percentage
#' by sample.
#'
#' @param print_result Logical. If `TRUE` (default), prints the concordance
#' results data frame to the console. If `FALSE`, results are returned invisibly.
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{ID}: Sample identifiers shared between the datasets.
#'   \item \code{Concordance}: Percentage of matching genotypes per sample.
#' }
#' If \code{print_result = FALSE}, the data frame is returned invisibly.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Identifies common samples and markers between the datasets.
#'   \item Optionally excludes specified SNPs.
#'   \item Removes loci with missing data (if \code{missing_code} is provided).
#'   \item Computes per-sample concordance as the percentage of matching genotypes.
#' }
#'
#' When \code{plot = TRUE}, a bar plot showing concordance percentage per sample
#' is generated using \pkg{ggplot2}.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#' result <- imputation_concordance(
#'   reference_genos = ref,
#'   imputed_genos = test,
#'   snps_2_exclude = snps,
#'   missing_code = 5,
#'   verbose = TRUE,
#'   plot = TRUE
#' )
#'
#' @export
imputation_concordance <- function(reference_genos,
                                   imputed_genos,
                                   missing_code = NULL,
                                   snps_2_exclude = NULL,
                                   verbose = FALSE,
                                   plot = FALSE,
                                   print_result = TRUE) {

  # Find common IDs
  common_ids <- intersect(imputed_genos$ID, reference_genos$ID)

  imputed_genos <- imputed_genos %>%
    filter(ID %in% common_ids) %>%
    arrange(ID)

  reference_genos <- reference_genos %>%
    filter(ID %in% common_ids) %>%
    arrange(ID)

  # Find common SNPs
  common_snps <- setdiff(
    intersect(colnames(imputed_genos), colnames(reference_genos)),
    as.vector(unlist(snps_2_exclude))
  )

  # Remove ID column if present
  common_snps <- setdiff(common_snps, "ID")

  # Convert to matrices
  imputed_matrix <- as.matrix(imputed_genos[, common_snps])
  reference_matrix <- as.matrix(reference_genos[, common_snps])

  # Identify valid SNPs
  if (!is.null(missing_code)) {
    valid_snps <- (imputed_matrix != missing_code) &
      (reference_matrix != missing_code)
  } else {
    valid_snps <- matrix(TRUE,
                         nrow = nrow(imputed_matrix),
                         ncol = ncol(imputed_matrix))
  }

  # Compute concordance
  matches <- (imputed_matrix == reference_matrix) & valid_snps
  percentage_match <- rowSums(matches, na.rm = TRUE) /
    rowSums(valid_snps, na.rm = TRUE)

  percentage_match[is.nan(percentage_match)] <- NA

  # Output data frame (original structure preserved)
  result_df <- data.frame(
    ID = imputed_genos$ID,
    Concordance = paste0(round(percentage_match * 100, 2), "%")
  )

  # Summary statistics
  summary_concordance <- summary(percentage_match, na.rm = TRUE) * 100
  names(summary_concordance) <- c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")

  if (verbose) {
    cat("Concordance Summary:\n")
    for (name in names(summary_concordance)) {
      cat(name, ":", round(summary_concordance[name], 2), "%\n")
    }
  }

  # Print results to console (NEW OPTION)
  if (print_result) {
    print(result_df)
  }

  # Optional plot
  if (plot) {

    plot_df <- data.frame(
      ID = imputed_genos$ID,
      Concordance = percentage_match * 100
    )

    concordance_plot <- ggplot(plot_df,
                               aes(x = reorder(ID, Concordance),
                                   y = Concordance)) +
      geom_bar(stat = "identity") +
      labs(title = "Imputation Concordance by Sample",
           x = "Sample ID",
           y = "Concordance (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    print(concordance_plot)
  }

  if (print_result) {
    return(result_df)
  } else {
    invisible(result_df)
  }
}
