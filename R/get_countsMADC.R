#' Obtain Read Counts from MADC File
#'
#' This function takes the MADC file as input and retrieves the ref and alt counts for each sample,
#' and converts them to ref, alt, and size(total count) matrices for dosage calling tools. At the moment,
#' only the read counts for the Ref and Alt target loci are obtained while the additional loci are ignored.
#'
#'
#' @param madc_file Path to MADC file
#' @import dplyr
#' @return A list of read count matrices for reference, alternate, and total read count values
#' @examples
#' # Get the path to the MADC file
#' madc_path <- system.file("iris_DArT_MADC.csv", package = "BIGr")
#'
#' # Extract the read count matrices
#' counts_matrices <- get_countsMADC(madc_path)
#'
#' # Access the reference, alternate, and size matrices
#'
#' # ref_matrix <- counts_matrices$ref_matrix
#' # alt_matrix <- counts_matrices$alt_matrix
#' # size_matrix <- counts_matrices$size_matrix
#'
#' rm(counts_matrices)
#' @export
get_countsMADC <- function(madc_file) {
  # This function takes the MADC file as input and generates a Ref and Alt counts dataframe as output
  update_df <- get_counts(madc_file)

  # Filter rows where 'AlleleID' ends with 'Ref'
  ref_df <- subset(update_df, grepl("Ref$", AlleleID))

  # Filter rows where 'AlleleID' ends with 'Alt'
  alt_df <- subset(update_df, grepl("Alt$", AlleleID))

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

  #Define columns to remove
  columns_to_remove <- c("AlleleID", "CloneID", "AlleleSequence","ClusterConsensusSequence",
                        "CallRate","OneRatioRef","OneRatioSnp","FreqHomRef","FreqHomSnp",
                        "FreqHets", "PICRef","PICSnp","AvgPIC","AvgCountRef","AvgCountSnp",
                        "RatioAvgCountRefAvgCountSnp") #Adjust as needed for different MADC versions

  #Identify columns that are in MADC
  col_remove <- intersect(columns_to_remove, colnames(ref_df))

  #Remove the specified columns and convert to matrix
  ref_matrix <- as.matrix(select(ref_df, -all_of(col_remove)))
  alt_matrix <- as.matrix(select(alt_df, -all_of(col_remove)))

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
  matrices_list <- list(ref_matrix = ref_matrix, size_matrix = size_matrix)

  return(matrices_list)

}

get_counts <- function(madc_file) {
  # Read the MADC file
  #Read only the first column for the first seven rows
  first_seven_rows <- read.csv(madc_file, header = FALSE, nrows = 7, colClasses = c(NA, "NULL"))

  #Check if all entries in the first column are either blank or "*"
  check_entries <- all(first_seven_rows[, 1] %in% c("", "*"))

  #Check if the MADC file has the filler rows or is processed from updated fixed allele ID pipeline
  if (check_entries) {
    #Note: This assumes that the first 7 rows are placeholder info from DArT processing

    #Read the madc file
    madc_df <- read.csv(madc_file, sep = ',', skip = 7, check.names = FALSE)

    #Retain only the Ref and Alt haplotypes
    filtered_df <- madc_df[!grepl("\\|AltMatch|\\|RefMatch", madc_df$AlleleID), ]

    #Remove extra text after Ref and Alt (_001 or _002)
    filtered_df$AlleleID <- sub("\\|Ref.*", "|Ref", filtered_df$AlleleID)
    filtered_df$AlleleID <- sub("\\|Alt.*", "|Alt", filtered_df$AlleleID)

  } else {

    #Read the madc file
    madc_df <- read.csv(madc_file, sep = ',', check.names = FALSE)

    # Retain only the Ref and Alt haplotypes
    filtered_df <- madc_df[!grepl("\\|AltMatch|\\|RefMatch", madc_df$AlleleID), ]

    #Remove extra text after Ref and Alt (_001 or _002)
    filtered_df$AlleleID <- sub("\\|Ref.*", "|Ref", filtered_df$AlleleID)
    filtered_df$AlleleID <- sub("\\|Alt.*", "|Alt", filtered_df$AlleleID)

  }

  return(filtered_df)
}
