#' Thin a dataframe of SNPs based on genomic position
#'
#' This function groups SNPs by chromosome, sorts them by physical position,
#' and then iteratively selects SNPs such that no two selected SNPs within
#' the same chromosome are closer than a specified minimum distance.
#'
#' @param df The input dataframe.
#' @param chrom_col_name A string specifying the name of the chromosome column.
#' @param pos_col_name A string specifying the name of the physical position column.
#' @param min_distance A numeric value for the minimum distance between selected SNPs.
#'   The unit of this distance should match the unit of the `pos_col_name` column (e.g., base pairs).
#'
#' @import dplyr
#' @import rlang
#' @return A thinned dataframe with the same columns as the input.
#'
#' @examples
#' # Create sample SNP data
#' set.seed(123)
#' n_snps <- 20
#' snp_data <- data.frame(
#'   MarkerID = paste0("SNP", 1:n_snps),
#'   Chrom = sample(c("chr1", "chr2"), n_snps, replace = TRUE),
#'   ChromPosPhysical = c(
#'     sort(sample(1:1000, 5)), # SNPs on chr1
#'     sort(sample(1:1000, 5)) + 500, # More SNPs on chr1
#'     sort(sample(1:2000, 10))      # SNPs on chr2
#'   ),
#'   Allele = sample(c("A/T", "G/C"), n_snps, replace = TRUE)
#' )
#' # Ensure it's sorted by Chrom and ChromPosPhysical for clarity in example
#' snp_data <- snp_data[order(snp_data$Chrom, snp_data$ChromPosPhysical), ]
#' rownames(snp_data) <- NULL
#'
#' print("Original SNP data:")
#' print(snp_data)
#'
#' # Thin the SNPs, keeping a minimum distance of 100 units (e.g., bp)
#' thinned_snps <- thinSNP(
#'   df = snp_data,
#'   chrom_col_name = "Chrom",
#'   pos_col_name = "ChromPosPhysical",
#'   min_distance = 100
#' )
#'
#' print("Thinned SNP data (min_distance = 100):")
#' print(thinned_snps)
#'
#' # Thin with a larger distance
#' thinned_snps_large_dist <- thinSNP(
#'   df = snp_data,
#'   chrom_col_name = "Chrom",
#'   pos_col_name = "ChromPosPhysical",
#'   min_distance = 500
#' )
#' print("Thinned SNP data (min_distance = 500):")
#' print(thinned_snps_large_dist)
#' @export
thinSNP <- function(df, chrom_col_name, pos_col_name, min_distance) {
  # Convert column name strings to symbols for dplyr
  chrom_sym <- rlang::sym(chrom_col_name)
  pos_sym <- rlang::sym(pos_col_name)

  df %>%
    dplyr::group_by(!!chrom_sym) %>%
    dplyr::arrange(!!pos_sym, .by_group = TRUE) %>%
    # Apply thinning logic to each chromosome group
    dplyr::group_modify(~ {
      # .x is the subset of data for the current chromosome, already sorted by position
      if (nrow(.x) == 0) {
        # Return an empty tibble with the same structure if the group is empty
        return(tibble::as_tibble(.x[0, , drop = FALSE]))
      }

      # Vector to store indices of rows to keep
      kept_indices <- integer(nrow(.x))
      num_kept <- 0

      # Always keep the first SNP in the sorted group
      num_kept <- num_kept + 1
      kept_indices[num_kept] <- 1
      last_selected_pos <- .x[[pos_col_name]][1] # Get position of the first SNP

      # Iterate through the rest of the SNPs in the current chromosome
      if (nrow(.x) > 1) {
        for (i in 2:nrow(.x)) {
          current_pos <- .x[[pos_col_name]][i]
          # If current SNP is far enough from the last selected SNP, keep it
          if (current_pos >= last_selected_pos + min_distance) {
            num_kept <- num_kept + 1
            kept_indices[num_kept] <- i
            last_selected_pos <- current_pos
          }
        }
      }
      # Subset the group to include only the kept SNPs
      .x[kept_indices[1:num_kept], , drop = FALSE]
    }) %>%
    dplyr::ungroup() # Ungroup to return a single dataframe
}
