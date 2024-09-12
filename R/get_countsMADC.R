#' Obtain Read Counts from MADC File
#'
#' This function takes the MADC file as input and retrieves the ref and alt counts for each sample,
#' and converts them to ref, alt, and size(total count) matrices for dosage calling tools. At the moment,
#' only the read counts for the Ref and Alt target loci are obtained while the additional loci are ignored.
#'
#'
#' @param madc_file Path to MADC file
#' @return A list of read count matrices for reference, alternate, and total read count values
#' @export
get_countsMADC <- function(madc_file) {
  # This function takes the MADC file as input and generates a Ref and Alt counts dataframe as output

  # Read the madc file
  madc_df <- read.csv(madc_file, sep = ',', skip = 7, check.names = FALSE)

  # Retain only the Ref and Alt haplotypes
  filtered_df <- madc_df[grep("\\|Ref$|\\|Alt$", madc_df$AlleleID), ]

  # Filter rows where 'AlleleID' ends with 'Ref'
  ref_df <- subset(filtered_df, grepl("Ref$", AlleleID))

  # Filter rows where 'AlleleID' ends with 'Alt'
  alt_df <- subset(filtered_df, grepl("Alt$", AlleleID))

  #Ensure that each has the same SNPs and that they are in the same order
  same <- identical(alt_df$CloneID,ref_df$CloneID)

  ###Convert the ref and alt counts into matrices with the CloneID as the index
  #Set SNP names as index
  row.names(ref_df) <- ref_df$CloneID
  row.names(alt_df) <- alt_df$CloneID

  #Retain only the rows in common if they are not identical and provide warning
  if (same == FALSE) {
    warning("Mismatch between Ref and Alt Markers. MADC likely altered. Markers without a Ref or Alt match removed.")
    # Find the common CloneIDs between the two dataframes
    common_ids <- intersect(rownames(ref_df), rownames(alt_df))
    # Subset both dataframes to retain only the common rows
    ref_df <- ref_df[common_ids, ]
    alt_df <- alt_df[common_ids, ]
  }

  #Remove unwanted columns and convert to matrix
  ref_matrix <- as.matrix(ref_df[, -c(1:16)])
  alt_matrix <- as.matrix(alt_df[, -c(1:16)])

  #Convert elements to numeric
  class(ref_matrix) <- "numeric"
  class(alt_matrix) <- "numeric"

  #Make the size matrix by combining the two matrices
  size_matrix <- (ref_matrix + alt_matrix)

  #Count the number of cells with 0 count to estimate missing data
  # Count the number of cells with the value 0
  count_zeros <- sum(size_matrix == 0)

  # Print the result
  ratio_missing_data <- count_zeros / length(size_matrix)
  cat("Ratio of missing data =", ratio_missing_data, "\n")

  # Return the ref and alt matrices as a list
  matrices_list <- list(ref_matrix = ref_matrix, alt_matrix = alt_matrix, size_matrix = size_matrix)
  return(matrices_list)
}
