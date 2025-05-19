#Internal Functions

globalVariables(c(
  "ALT", "AlleleID", "CHROM", "Data", "ID", "MarkerName", "POS",
  "QPseparate", "QPsolve_par", "REF", "Var1", "Variant", "geno",
  "ind", "ref", "row_name", "size", "snp"
))

#' Convert GT format to numeric dosage
#' @param gt a genotype matrix with samples as columns and variants as rows
#' @return numeric genotype values
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

#' Check and Adjust Botloci and MADC Marker Compatibility
#'
#' This internal function checks the compatibility between botloci and MADC markers. It ensures that the marker IDs in the botloci file match those in the MADC file. If discrepancies are found, such as mismatched padding, the function attempts to adjust the IDs to ensure compatibility.
#'
#' @param botloci A data frame containing the botloci markers.
#' @param report A data frame containing the MADC markers.
#' @param verbose A logical value indicating whether to print detailed messages about the adjustments. Default is TRUE.
#'
#' @return A list containing the adjusted botloci and MADC data frames.
#'
#' @details
#' The function checks if the marker IDs in the botloci file are present in the MADC file. If no matches are found, it examines the padding (number of digits) in the marker IDs and adjusts them to match the longest padding. If the IDs still do not match after adjustment, an error is raised. This function is intended for internal use and helps ensure that the botloci and MADC files are compatible for downstream analysis.
#'
#' @keywords internal
#' @noRd
check_botloci <- function(botloci, report, verbose=TRUE){
  if(!any(botloci$V1 %in% report$CloneID)) {
    if(verbose) message("None of the botloci markers could be found in the MADC file. Checking padding match...\n")

    pad_madc <- unique(nchar(sub(".*_", "", report$CloneID)))
    pad_botloci <- unique(nchar(sub(".*_", "", botloci$V1)))

    if(length(pad_madc) > 1 | length(pad_botloci) > 1) stop("Check marker IDs in both MADC and botloci files. They should be the same.")

    if(pad_madc != pad_botloci) {
      if(verbose) message("Padding between MADC and botloci files do not match. Markers ID modified to match longest padding.\n")
      if (pad_madc < pad_botloci) {
        report$CloneID <- paste0(sub("_(.*)", "", report$CloneID), "_",
                                 sprintf(paste0("%0", pad_botloci, "d"), as.integer(sub(".*_", "", report$CloneID)))
        )
      } else {
        botloci$V1 <- paste0(sub("_(.*)", "", botloci$V1), "_",
                             sprintf(paste0("%0", pad_madc, "d"), as.integer(sub(".*_", "", botloci$V1)))
        )
        if(!any(botloci$V1 %in% report$CloneID)) stop("After matching padding, botloci markers still not found in MADC file. Check marker IDs.\n")
      }
    } else {
      stop("Check marker IDs in both MADC and botloci files. They should be the same.")
    }
  }
  return(list(botloci, report))
}
