#' Filter MADC Files
#'
#' Filter and process MADC files to remove low quality microhaplotypes
#'
#' @details
#' This function can filter raw MADC files or pre-processed MADC files with fixed allele IDs. Additionally,
#' it can filter based on mean read depth, number of mhaps per target loci, and other criteria. Optionally, users
#' can scale and normalize the data in preparation for conversion to relationship matrices,
#' plot summary statistics, and save the filtered data to a file.
#'
#'@import dplyr
#'@importFrom utils read.csv
#'
#'@param madc_file Path to the MADC file to be filtered
#'@param min.mean.reads Minimum mean read depth for filtering
#'@param max.mean.reads Maximum mean read depth for filtering
#'@param max.match.mhaps Maximum number of matching mhaps per target loci
#'@param min.reads.per.site Minimum number of reads per site for filtering
#'@param min.ind.with.reads Minimum number of individuals with reads for filtering
#'@param target_only Logical indicating whether to filter for target loci only
#'@param fixed_allele_ids Logical indicating whether the MADC file has been pre-processed for fixed allele IDs
#'@param plot.summary Logical indicating whether to plot summary statistics
#'@param output.file Path to save the filtered data (if NULL, data will not be saved)
#'@param verbose Logical indicating whether to print additional information during processing
#'
#'@return data.frame or saved csv file
#'
#'@examples
#' #Example...
#'
#' ##Plots
#' #Mean read depth
#' #Number of Altmatch and Refmatch mhaps per target loci
#'
#'
#'@export
filterMADC <- function(madc_file,
                       min.mean.reads = NULL,
                       max.mean.reads = NULL,
                       max.match.mhaps = 10,
                       min.reads.per.site = NULL,
                       min.ind.with.reads = NULL,
                       target_only = FALSE,
                       fixed_allele_ids = FALSE,
                       plot.summary = FALSE,
                       output.file = NULL) {


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

    #Remove extra text after Ref and Alt (_001 or _002)
    filtered_df$AlleleID <- sub("\\|Ref_.*", "|Ref", filtered_df$AlleleID)
    filtered_df$AlleleID <- sub("\\|Alt_.*", "|Alt", filtered_df$AlleleID)

  } else {

    #Read the madc file
    filtered_df <- read.csv(madc_file, sep = ',', check.names = FALSE)

    #Remove extra text after Ref and Alt (_001 or _002)
    filtered_df$AlleleID <- sub("\\|Ref_.*", "|Ref", filtered_df$AlleleID)
    filtered_df$AlleleID <- sub("\\|Alt_.*", "|Alt", filtered_df$AlleleID)

  }

  #Remove refmatch and altmatch if wanted
  if (target_only) {
    message("Retaining target markers only")
    #Retain only the Ref and Alt haplotypes
    filtered_df <- filtered_df[!grepl("\\|AltMatch|\\|RefMatch", filtered_df$AlleleID), ]
  }

  ## Filtering
  if (!is.null(min.mean.reads)) {
    message("Filtering for minimum mean reads across all samples")
    filtered_df <- filtered_df[filtered_df$MeanReads >= min.mean.reads, ]
  }

  #Save the output to disk if file name provided
  if (!is.null(output.file)) {
    message("Saving filtered data to file")
    write.csv(filtered_df, paste0(output.file,".csv"), row.names = FALSE)
  }

  return(filtered_df)
}
