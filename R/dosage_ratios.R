#' Calculate the Percentage of Each Dosage Value
#'
#' This function calculates the percentage of each dosage value within a genotype matrix.
#' It assumes that the samples are the columns, and the genomic markers are in rows. Missing data should
#' be set as NA, which will then be ignored for the calculations. All samples must have the same ploidy.
#'
#' @param data Genotype matrix or data.frame
#' @param ploidy The ploidy of the species being analyzed
#' @return A data.frame with percentages of dosage values in the genotype matrix
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
  percentages_df$Data <- "Dosages"
  melted_data <- percentages_df %>%
    pivot_longer(cols = -(Data),names_to = "Dosage", values_to = "Percentage")

  return(melted_data)
}

