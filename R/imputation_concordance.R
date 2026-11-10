#' Calculate concordance between imputed and reference genotypes
#'
#' Compares imputed and reference genotype datasets sample by sample, returning
#' the percentage of matching genotypes per sample. Expects samples as rows and
#' markers as columns, with allele dosages (0, 1, 2) recommended. Other numeric
#' formats are supported.
#'
#' @param reference_genos A data frame of reference genotypes (samples × markers)
#'   with a column named `ID`.
#' @param imputed_genos A data frame of imputed genotypes (samples × markers)
#'   with a column named `ID`.
#' @param missing_code Optional value indicating missing data. Loci carrying this
#'   value in either dataset are excluded from the concordance calculation.
#' @param snps_2_exclude Optional character vector of marker names to exclude.
#' @param verbose Logical. Print a five-number summary of concordance? Default `FALSE`.
#' @param plot Logical. Produce a bar plot of concordance by sample? Default `FALSE`.
#' @param print_result Logical. Print the results data frame to the console?
#'   Default `TRUE`. If `FALSE`, results are returned invisibly.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{ID}{Sample identifiers shared between the two datasets.}
#'   \item{Concordance}{Percentage of matching genotypes per sample.}
#' }
#'
#' @details
#' The function identifies common samples and markers between the two datasets,
#' optionally removes specified SNPs and loci with missing data, then computes
#' per-sample concordance as the percentage of matching genotypes over valid loci.
#' When `plot = TRUE`, a bar plot is produced with \pkg{ggplot2}.
#'
#' @examples
#' ref <- data.frame(
#'   ID   = c("S1", "S2", "S3"),
#'   SNP1 = c(0, 1, 2),
#'   SNP2 = c(1, 1, 0),
#'   SNP3 = c(2, 5, 1)
#' )
#' imp <- data.frame(
#'   ID   = c("S1", "S2", "S3"),
#'   SNP1 = c(0, 0, 2),
#'   SNP2 = c(1, 1, 1),
#'   SNP3 = c(2, 5, 0)
#' )
#' imputation_concordance(
#'   reference_genos = ref,
#'   imputed_genos   = imp,
#'   snps_2_exclude  = "SNP2",
#'   missing_code    = 5,
#'   print_result    = FALSE
#' )
#'
#' @author Josué Chinchilla-Vargas
#'
#' @importFrom dplyr %>% filter arrange
#' @importFrom stats reorder
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal theme element_text
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
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required when plot = TRUE.", call. = FALSE)
    }

    plot_df <- data.frame(
      ID = imputed_genos$ID,
      Concordance = percentage_match * 100
    )

    concordance_plot <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = reorder(ID, Concordance), y = Concordance)
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(title = "Imputation Concordance by Sample",
                    x = "Sample ID",
                    y = "Concordance (%)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

    print(concordance_plot)
  }

  if (print_result) {
    return(result_df)
  } else {
    invisible(result_df)
  }
}
