#' Convert MADC Files to an Additive Genomic Relationship Matrix
#'
#' Scale and normalize MADC read count data and convert it to an additive genomic relationship matrix.
#'
#'@details
#' This function reads a MADC file, processes it to remove unnecessary columns, scales and normalizes the data, and
#' then converts it into an additive genomic relationship matrix using the `A.mat` function from the `rrBLUP` package.
#' The resulting matrix can be used for genomic selection or other genetic analyses.
#'
#'@import dplyr
#'@importFrom utils read.csv write.csv
#'@importFrom rrBLUP A.mat
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

  #Remove extra text after Ref and Alt (_001 or _002)
  filtered_df$AlleleID <- sub("\\|Ref_001*", "|Ref", filtered_df$AlleleID)
  filtered_df$AlleleID <- sub("\\|Alt_002", "|Alt", filtered_df$AlleleID)

  #Removing extra columns
  row.names(filtered_df) <- filtered_df$AlleleID
  filtered_df <- filtered_df %>%
    select(-c(AlleleID, CloneID, AlleleSequence))


  #Scale and normalized data
  message("Scaling and normalizing data to be -1,1")
  # Function to scale a matrix to be between -1 and 1 for rrBLUP
  scale_matrix <- function(mat) {
    min_val <- min(mat)
    max_val <- max(mat)

    # Normalize to [0, 1]
    normalized <- (mat - min_val) / (max_val - min_val)

    # Scale to [-1, 1]
    scaled <- 2 * normalized - 1

    return(scaled)
  }

  # Apply the scaling function
  filtered_df <- scale_matrix(filtered_df)

  #Making additive relationship matrix
  MADC.mat <- A.mat(t(filtered_df))

  rm(filtered_df)

  #Save the output to disk if file name provided
  if (!is.null(output.file)) {
    message("Saving filtered data to file")
    write.csv(MADC.mat, paste0(output.file,".csv"), row.names = TRUE)
  } else {
    return(MADC.mat)
  }

}
