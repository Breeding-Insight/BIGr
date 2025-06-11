#' Convert MADC Files to an Additive Genomic Relationship Matrix
#'
#' Scale and normalize MADC read count data and convert it to an additive genomic relationship matrix.
#'
#'@details
#' This function reads a MADC file, processes it to remove unnecessary columns, scales and normalizes the data, and
#' then converts it into an additive genomic relationship matrix using the `A.mat` function from the `rrBLUP` package.
#' The resulting matrix can be used for genomic selection or other genetic analyses.
#'
#'@importFrom utils read.csv write.csv
#'@importFrom rrBLUP A.mat
#'@import tibble stringr dplyr tidyr
#'
#'@param madc_file Path to the MADC file to be filtered
#'@param seed Optional seed for random number generation (default is NULL)
#'@param output.file Path to save the filtered data (if NULL, data will not be saved)
#'
#'@return data.frame or saved csv file
#'
#'@examples
#' #Input variables
#' madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
#'
#' #Calculations
#' temp <- tempfile()
#'
#' # Converting to additive relationship matrix
#' gmat <- madc2gmat(madc_file,
#'                  seed = 123,
#'                  output.file = NULL)
#'
#'@references
#'Endelman, J. B. (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. The Plant Genome, 4(3).
#'
#'@export
madc2gmat <- function(madc_file,
                      seed = NULL,
                      output.file = NULL) {
  #set seed if not null
  if (!is.null(seed)) {
    set.seed(seed)
  }

  #Need to first inspect the first 7 rows of the MADC to see if it has been preprocessed or not
  first_seven_rows <- read.csv(madc_file, header = FALSE, nrows = 7, colClasses = c(NA, "NULL"))

  #Check if all entries in the first column are either blank or "*"
  check_entries <- all(first_seven_rows[, 1] %in% c("", "*"))

  #Check if the MADC file has the filler rows or is processed from updated fixed allele ID pipeline
  if (check_entries) {
    #Note: This assumes that the first 7 rows are placeholder info from DArT processing

    warning("The MADC file has not been pre-processed for Fixed Allele IDs. The first 7 rows are placeholder info from DArT processing.")

    #Read the madc file
    filtered_df <- read.csv(madc_file, sep = ',', skip = 7, check.names = FALSE)

  } else {

    #Read the madc file
    filtered_df <- read.csv(madc_file, sep = ',', check.names = FALSE)
  }

  #Convert to data.table
  #filtered_df <- as.data.table(filtered_df)

  #Get allele ratios
  # --- Step 1: Calculate frequency of each specific allele using total locus count ---
  # The denominator MUST be the sum of all reads at the locus.
  allele_freq_df <- filtered_df %>%
    group_by(CloneID) %>%
    mutate(across(
      .cols = where(is.numeric),
      # The calculation is count / total_locus_count
      .fns = ~ . / sum(.),
      .names = "{.col}_freq"
    )) %>%
    ungroup()

  #Rm object
  rm(filtered_df)

  # --- Step 2: Prepare data for reshaping ---
  # Filter out the reference rows and create a unique marker name for each alternative allele
  markers_to_pivot <- allele_freq_df %>%
    filter(!str_ends(AlleleID, "\\|Ref_0001")) %>%
    # Create a new, unique column for each marker (e.g., "chr1.1_000194324_RefMatch_0001")
    # This makes each alternative allele its own column in the final matrix.
    mutate(MarkerID = str_replace(AlleleID, "\\|", "_")) %>%
    select(MarkerID, ends_with("_freq")) %>%
    tibble::column_to_rownames(var = "MarkerID")
  names(markers_to_pivot) <- str_remove(names(markers_to_pivot), "_freq$")
  markers_to_pivot <- t(markers_to_pivot)
  #Replace the NaN to NA
  markers_to_pivot[is.nan(markers_to_pivot)] <- NA

  #rm unneeded objects
  rm(allele_freq_df)

  #Step 3 is a robust transformation step, but probably not needed.
  # --- Step 3: Reshape data into the Marker Matrix format ---
  #marker_matrix_for_Amat <- markers_to_pivot %>%
  #  pivot_longer(
  #    cols = -MarkerID,
  #    names_to = "SampleID",
  #    values_to = "AlleleFreq"
  #  ) %>%
  #  mutate(SampleID = str_remove(SampleID, "_freq$")) %>%
  #  pivot_wider(
  #    names_from = MarkerID,
  #    values_from = AlleleFreq,
  # If a sample has zero reads for an allele, it won't appear after the filter.
  # values_fill = 0 ensures these are explicitly set to zero frequency.
  #    values_fill = 0
  #  ) %>%
  #  tibble::column_to_rownames("SampleID") %>%
  #  as.matrix()


  #Scale and normalized data
  message("Scaling and normalizing data to be -1,1")
  # Function to scale a matrix to be between -1 and 1 for rrBLUP
  scale_matrix <- function(mat) {
    #min_val <- min(mat)
    #max_val <- max(mat)

    # Normalize to [0, 1]
    #normalized <- (mat - min_val) / (max_val - min_val)

    # Scale to [-1, 1]
    #scaled <- 2 * normalized - 1
    scaled <- 2 * mat - 1

    return(scaled)
  }

  # Apply the scaling function
  markers_to_pivot <- scale_matrix(markers_to_pivot)

  #Making additive relationship matrix
  MADC.mat <- A.mat(markers_to_pivot)

  rm(markers_to_pivot)

  #Save the output to disk if file name provided
  if (!is.null(output.file)) {
    message("Saving filtered data to file")
    write.csv(MADC.mat, paste0(output.file,".csv"), row.names = TRUE)
  } else {
    return(MADC.mat)
  }

}
