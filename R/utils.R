#Internal Functions

globalVariables(c(
  "ALT", "AlleleID", "CHROM", "Data", "ID", "MarkerName", "POS", 
  "QPseparate", "QPsolve_par", "REF", "Var1", "Variant", "geno", 
  "ind", "ref", "row_name", "size", "snp"
))

#' Convert GT format to numeric dosage
#' @param gt a genotype matrix with samples as columns and variants as rows
#' @noRd
convert_to_dosage <- function(gt) {
  # Split the genotype string
  alleles <- strsplit(gt, "[|/]")
  # Sum the alleles, treating NA values appropriately
  sapply(alleles, function(x) {
    if (any(is.na(x))) {
      return(NA)
    } else {
      return(sum(as.numeric(x), na.rm = TRUE))
    }
  })
}
