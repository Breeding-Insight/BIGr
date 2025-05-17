#' Filter MADC Files
#'
#' Filter and process MADC files to remove low quality microhaplotypes
#'
#' @details
#' This function can filter raw MADC files or pre-processed MADC files with fixed allele IDs. Additionally,
#' it can filter based on mean read depth, number of mhaps per target loci, and other criteria. Optionally, users
#' can plot summary statistics and save the filtered data to a file.
#'
#'@import dplyr
#'@importFrom utils read.csv
#'
#'@param madc_file Path to the MADC file to be filtered
#'@param min.mean.reads Minimum mean read depth for filtering
#'@param max.mean.reads Maximum mean read depth for filtering
#'@param max.mhaps.per.loci Maximum number of matching mhaps per target loci. Retains only the target Ref and Alt loci at the sites that exceeds the \code{max.mhaps.per.loci} threshold.
#'@param min.reads.per.site Minimum number of reads per site for \code{min.ind.with.reads}. Otherwise, this parameter is ignored
#'@param min.ind.with.reads Minimum number of individuals with \code{min.reads.per.site} reads for filtering
#'@param target.only Logical indicating whether to filter for target loci only
#'@param n.summary.columns (optional) Number of summary columns to remove from MADC file not including the first three. Otherwise, the columns will be automatically detected and removed.
#'@param plot.summary Logical indicating whether to plot summary statistics
#'@param output.file Path to save the filtered data (if NULL, data will not be saved)
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
                       max.mhaps.per.loci = NULL,
                       min.reads.per.site = 1,
                       min.ind.with.reads = NULL,
                       target.only = FALSE,
                       n.summary.columns = NULL,
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
  #Check for extra columns
  #Save the three columns for later adding to the output
  saved_columns <- filtered_df[,1:3]

  if (!is.null(n.summary.columns)) {
    #Remove the first n.summary.columns columns
    filtered_df <- filtered_df[,-c(4:n.summary.columns)]
  }else{
    rm.col <- c("ClusterConsensusSequence",
                "CallRate", "OneRatioRef", "OneRatioSnp", "FreqHomRef", "FreqHomSnp",
                "FreqHets", "PICRef", "PICSnp", "AvgPIC", "AvgCountRef", "AvgCountSnp","RatioAvgCountRefAvgCountSnp")

    filtered_df <- filtered_df[, !(colnames(filtered_df) %in% rm.col)]
  }

  #Now add rownames
  rownames(filtered_df) <- saved_columns[,1]

  #Remove refmatch and altmatch if wanted
  if (target.only) {
    message("Retaining target markers only")
    #Retain only the Ref and Alt haplotypes
    filtered_df <- filtered_df[!grepl("\\|AltMatch|\\|RefMatch", filtered_df$AlleleID), ]
  }

  ## Filtering

  #Min mean reads
  if (!is.null(min.mean.reads)) {
    message("Filtering for minimum mean reads across all samples")
    #Get the mean value for each row, and remove the rows below that threshold
    filtered_df$MeanReads <- rowMeans(filtered_df[, -c(1:3)], na.rm = TRUE)
    filtered_df <- filtered_df[filtered_df$MeanReads >= min.mean.reads, ]
    #Remove the MeanReads column
    filtered_df <- filtered_df[, -which(colnames(filtered_df) == "MeanReads")]
  }

  #Max mean reads
  if (!is.null(max.mean.reads)) {
    message("Filtering for maximum mean reads across all samples")
    #Get the mean value for each row, and remove the rows above that threshold
    filtered_df$MeanReads <- rowMeans(filtered_df[, -c(1:3)], na.rm = TRUE)
    filtered_df <- filtered_df[filtered_df$MeanReads <= max.mean.reads, ]
    #Remove the MeanReads column
    filtered_df <- filtered_df[, -which(colnames(filtered_df) == "MeanReads")]
  }

  #Max mhaps per loci
  if (!is.null(max.mhaps.per.loci)) {
    message("Filtering for maximum number of matching mhaps per target loci")
    #Get the number of matching mhaps for loci, and remove the mhaps at those loci that exceed the max number
    mhap_counts <- filtered_df %>%
      group_by(CloneID) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      filter(Count > max.mhaps.per.loci)

    patterns_to_search <- "\\|AltMatch|\\|RefMatch"
    clone_ids_to_target <- mhap_counts$CloneID

    filtered_df <- filtered_df %>%
      filter(
        !( # "keep rows that DO NOT match both conditions"
          CloneID %in% clone_ids_to_target &  # Condition 1: CloneID is one of the targeted IDs
            grepl(patterns_to_search, AlleleID) # Condition 2: AlleleID contains one of the patterns
        )
      )
  }

  #Min individuals with reads
  if (!is.null(min.ind.with.reads)) {
    message("Filtering for minimum number of individuals with reads per site")
    message(past0("Minimum number of individuals with reads per site: ", min.ind.with.reads))
    message(past0("Minimum number of reads per site: ", min.reads.per.site))

    #Getting colnames
    cols_to_check <- colnames(filtered_df)[-(1:3)]

    filtered_df <- filtered_df %>%
      rowwise() %>%  # Process data row by row
      mutate(
        # For each row, count how many of the 'cols_to_check' meet the criterion
        qualifying_sites_count = sum(
          c_across(all_of(cols_to_check)) >= min.reads.per.site,
          na.rm = TRUE # Treats NAs in data as not meeting the criterion
        )
      ) %>%
      ungroup() %>% # Always ungroup after rowwise operations
      # Filter rows where this count meets the 'min.ind.with.reads' threshold
      filter(qualifying_sites_count >= min.ind.with.reads) %>%
      # Optionally, remove the temporary count column if it's no longer needed
      select(-qualifying_sites_count)
  }

  #Plots
  if (plot.summary) {
    message("Plotting summary statistics")
    #Plot mean read depth
    mean_reads <- rowMeans(filtered_df[, -c(1:3)], na.rm = TRUE)
    hist(mean_reads, main = "Mean Read Depth", xlab = "Mean Reads", ylab = "Frequency")

    #Plot number of Altmatch and Refmatch mhaps per target loci
    altmatch_counts <- filtered_df %>%
      filter(grepl("\\|AltMatch", AlleleID)) %>%
      group_by(CloneID) %>%
      summarise(Count = n(), .groups = 'drop')

    refmatch_counts <- filtered_df %>%
      filter(grepl("\\|RefMatch", AlleleID)) %>%
      group_by(CloneID) %>%
      summarise(Count = n(), .groups = 'drop')

    barplot(cbind(altmatch_counts$Count, refmatch_counts$Count), beside = TRUE,
            names.arg = altmatch_counts$CloneID, main = "Number of AltMatch and RefMatch Mhaps",
            xlab = "Clone ID", ylab = "Count")

    #Plot density of number of CloneID per site on a marker distribution plot

  }

  #Save the output to disk if file name provided
  if (!is.null(output.file)) {
    message("Saving filtered data to file")
    write.csv(filtered_df, paste0(output.file,".csv"), row.names = FALSE)
  }

  return(filtered_df)
}
