#' Fix MADC File Allele IDs
#'
#' Process raw MADC files to format and update the allele IDs with user supplied Chr and Pos information
#'
#' @details
#' This function can process raw MADC files to update the Allele IDs and Clone IDs to the Chr_Pos format with a user supplied file.
#' The output MADC will be the standard fixed allele ID format to support use in madc2vcf and BIGapp functions.
#'
#'@import dplyr
#'@import stringr
#'@importFrom tidyr replace_na
#'@importFrom utils read.csv write.csv
#'
#'@param madc.file Path to the MADC file to be filtered
#'@param marker.file Path to the three column marker ID file.
#' - The first column is the existing list of unique CloneIDs (obtained from raw MADC CloneID),
#'where each row is a unique CloneID.
#' - The second column is the chromosome that the marker is located (ie Chr01). No special characters (*#_-!.) are permitted in the
#' chromosome name.
#' - The third column is the numeric position of the marker within the chromosome (ie 1234). No special characters (*#_-!.) are permitted.
#'@param n.summary.columns (optional) Number of summary columns to remove from MADC file not including the first three. Otherwise, the columns will be automatically detected and removed.
#'@param output.file Path to save the fixed allele ID MADC file (if NULL, data will not be saved)
#'
#'@return data.frame or saved csv file
#'
#'@examples
#' #Example
#'
#' #Example MADC
#' madc_file <- system.file("iris_DArT_MADC.csv", package="BIGr")
#' marker_file <- system.file("iris_MADC_marker_file.csv", package="BIGr")
#'
#' #Fix the raw MADC file IDs to use the user provided Chr_Pos format
#' fixedMADC_df <- fixMADC(madc.file = madc_file,
#'                          marker.file = marker_file,
#'                          n.summary.columns = NULL,
#'                          output.file = NULL)
#'
#'
#'
#'@export
fixMADC <- function(madc.file,
                       marker.file,
                       n.summary.columns = NULL,
                       output.file = NULL) {


  #Need to first inspect the first 7 rows of the MADC to see if it has been preprocessed or not
  first_seven_rows <- read.csv(madc.file, header = FALSE, nrows = 7, colClasses = c(NA, "NULL"))

  #Check if all entries in the first column are either blank or "*"
  check_entries <- all(first_seven_rows[, 1] %in% c("", "*"))

  #Check if the MADC file has the filler rows or is processed from updated fixed allele ID pipeline
  if (check_entries) {
    #Note: This assumes that the first 7 rows are placeholder info from DArT processing

    #Read the madc file
    filtered_df <- read.csv(madc.file, sep = ',', skip = 7, check.names = FALSE)

    #Remove extra text after Ref and Alt (_001 or _002)
    #filtered_df$AlleleID <- sub("\\|Ref_.*", "|Ref", filtered_df$AlleleID)
    #filtered_df$AlleleID <- sub("\\|Alt_.*", "|Alt", filtered_df$AlleleID)

  } else {

    stop("This MADC file appears to already use fixed allele IDs and cannot be reprocessed with a raw-ID marker file.")

    #Read the madc file
    filtered_df <- read.csv(madc.file, sep = ',', check.names = FALSE)

    #Remove extra text after Ref and Alt (_001 or _002)
    #filtered_df$AlleleID <- sub("\\|Ref_.*", "|Ref", filtered_df$AlleleID)
    #filtered_df$AlleleID <- sub("\\|Alt_.*", "|Alt", filtered_df$AlleleID)

  }
  #Check for extra columns
  #Save the three columns for later adding to the output (currently unused)

  if (!is.null(n.summary.columns)) {
    #Remove the first n.summary.columns columns
    if (n.summary.columns > 0) {
      cols_to_remove <- 4:(3 + n.summary.columns)
      filtered_df <- filtered_df[, -cols_to_remove, drop = FALSE]
    }
  }else{
    rm.col <- c("ClusterConsensusSequence",
                "CallRate", "OneRatioRef", "OneRatioSnp", "FreqHomRef", "FreqHomSnp",
                "FreqHets", "PICRef", "PICSnp", "AvgPIC", "AvgCountRef", "AvgCountSnp","RatioAvgCountRefAvgCountSnp")

    filtered_df <- filtered_df[, !(colnames(filtered_df) %in% rm.col)]
  }

  #Trim whitespace if present for some reason
  filtered_df$CloneID <- trimws(as.character(filtered_df$CloneID))
  filtered_df$AlleleID <- trimws(as.character(filtered_df$AlleleID))

  #Read in the marker file
  marker_file <- read.csv(marker.file, sep = ',', check.names = FALSE)

  ### Verify marker file is formatted correctly ###

  marker_file[,1] <- trimws(as.character(marker_file[,1]))
  marker_file[,2] <- trimws(as.character(marker_file[,2]))
  marker_file[,3] <- trimws(as.character(marker_file[,3]))

  if (any(grepl("[*#_!.\\-]", marker_file[,3]))) {
    stop("Special characters (*#_-!.) detected in the position column (column 3). Please review the marker file.")
  }

  if (!all(grepl("^[0-9]+$", marker_file[,3]))) {
    stop("The position column (column 3) must be numeric. Please review the marker file.")
  }

  #Make marker IDs column and pad 0's for the position
  marker_file$new_ID <- paste0(
    marker_file[,2],
    "_",
    str_pad(marker_file[,3], width = 9, side = "left", pad = "0")
  )

  #Verify there are no duplicate IDs in the marker file
  if (length(unique(marker_file[,1])) != length(marker_file[,1])) {
    stop("There are duplicate marker IDs in the first column. Please review marker file.")
  }

  #Verify there are no duplicate position information
  if (length(unique(marker_file$new_ID)) != length(marker_file$new_ID)) {
    stop("There are duplicate Chr and Pos information where more than one marker has the same Chr_Pos. Please review the marker file.")
  }

  #Verify chromosome column (col 2) contains no special characters
  if (any(grepl("[*#_!.\\-]", marker_file[,2]))) {
    stop("Special characters (*#_-!.) detected in the chromosome column (column 2). Please review the marker file.")
  }


  ## Filtering

  #Identify MADC CloneIDs not found in the marker file
  missing_ids <- setdiff(filtered_df$CloneID, marker_file[,1])

  if (length(missing_ids) > 0) {
    warning(paste0(
      length(missing_ids), " CloneID(s) in the MADC file were not found in the marker file and will be removed from the output:\n",
      paste(missing_ids, collapse = "\n")
    ))

    # Remove unmatched IDs from MADC file
    filtered_df <- filtered_df[filtered_df$CloneID %in% marker_file[,1], ]
  }

  ###Replace old IDs with new IDs

  #Create a named lookup vector: old CloneID and AlleleID -> new_ID
  id_lookup <- setNames(marker_file$new_ID, marker_file[,1])

  filtered_df <- filtered_df %>%
    mutate(
      # Replace CloneID column directly via lookup
      CloneID = id_lookup[CloneID],

      # For AlleleID: extract the |suffix, replace the ID prefix, rejoin
      AlleleID = {
        old_id <- str_extract(AlleleID, "^[^|]+")
        suffix  <- str_extract(AlleleID, "\\|.*$")
        suffix  <- replace_na(suffix, "")  # Handle cases with no suffix
        paste0(id_lookup[old_id], suffix)
      }
    )

  ###Add the proper numbering suffix to allele IDs for unique IDs.
  # |Ref    -> |Ref_0001
  # |Alt    -> |Alt_0002
  # |RefMatch -> numbered _0001, _0002, ... within each CloneID
  # |AltMatch -> numbered _0001, _0002, ... within each CloneID

  filtered_df <- filtered_df %>%
    mutate(
      .suffix_type = str_extract(AlleleID, "(?<=\\|)[^_]+$")
    ) %>%
    group_by(CloneID, .suffix_type) %>%
    mutate(
      AlleleID = case_when(
        .suffix_type == "Ref" ~ paste0(str_remove(AlleleID, "\\|Ref$"), "|Ref_0001"),
        .suffix_type == "Alt" ~ paste0(str_remove(AlleleID, "\\|Alt$"), "|Alt_0002"),
        .suffix_type == "RefMatch" ~ paste0(
          str_remove(AlleleID, "\\|RefMatch$"),
          "|RefMatch_",
          sprintf("%04d", row_number())
        ),
        .suffix_type == "AltMatch" ~ paste0(
          str_remove(AlleleID, "\\|AltMatch$"),
          "|AltMatch_",
          sprintf("%04d", row_number())
        ),
        TRUE ~ AlleleID
      )
    ) %>%
    ungroup() %>%
    select(-.suffix_type)


  #Save the output to disk if file name provided
  if (!is.null(output.file)) {
    message("Saving fixed MADC data to file")
    write.csv(filtered_df, paste0(output.file,"_fixedID.csv"), row.names = FALSE)
  } else {
    message("No output file provided. Returning fixed MADC data.")
    return(filtered_df)
  }

}
