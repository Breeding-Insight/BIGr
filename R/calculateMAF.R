#' Calculate Minor Allele Frequency from a Genotype Matrix
#'
#' This function calculates the allele frequency and minor allele frequency from a genotype matrix.
#' It assumes that the Samples are the columns, and the genomic markers are in rows. Missing data should
#' be set as NA, which will then be ignored for the calculations. All samples must have the same ploidy.
#'
#' @param df Genotype matrix or data.frame
#' @param ploidy The ploidy of the species being analyzed
#' @return A dataframe of AF and MAF values for each marker
#' @export
calculateMAF <- function(df, ploidy) {
  if (is.matrix(df)) {
    df <- as.data.frame(df)
  }
  
  allele_frequencies <- apply(df, 1, function(row) {
    non_na_count <- sum(!is.na(row))
    allele_sum <- sum(row, na.rm = TRUE)
    #print(paste("Non-NA count:", non_na_count, "Allele sum:", allele_sum))
    if (non_na_count > 0) {
      allele_sum / (ploidy * non_na_count)
    } else {
      NA
    }
  })
  
  maf <- ifelse(allele_frequencies <= 0.5, allele_frequencies, 1 - allele_frequencies)
  
  df$AF <- allele_frequencies
  df$MAF <- maf
  
  maf_df <- df[,c("AF", "MAF"), drop = FALSE]
  
  return(maf_df)
}