#' Obtain Read Counts from MADC File
#'
#' Reads a DArTag MADC report and returns reference and total read count matrices
#' per marker and sample. Only `Ref` and `Alt` target loci are retained;
#' `|AltMatch` / `|RefMatch` rows are either discarded or collapsed depending on
#' `collapse_matches_counts`.
#'
#' @details
#' Either `madc_file` or `madc_object` must be provided (not both `NULL`).
#' When `madc_object` is supplied it is passed directly to `get_counts()`,
#' skipping file I/O. The function constructs:
#' - `ref_matrix` — per-sample reference allele counts.
#' - `size_matrix` — per-sample total counts (ref + alt).
#'
#' Markers whose `CloneID` appears only in the `Ref` or only in the `Alt` rows
#' are removed with a warning. A summary of the proportion of zero-count
#' data points (missing data) is reported via `vmsg()`.
#'
#' @param madc_file character or `NULL`. Path to the input MADC CSV file.
#'   At least one of `madc_file` or `madc_object` must be provided.
#' @param madc_object data frame or `NULL`. A pre-read MADC data frame
#'   (e.g., as returned by `check_botloci()`). When supplied, file reading is
#'   skipped. At least one of `madc_file` or `madc_object` must be provided.
#' @param collapse_matches_counts logical. If `TRUE`, counts for `|AltMatch`
#'   and `|RefMatch` rows are summed into their corresponding `|Ref` and `|Alt`
#'   rows. If `FALSE` (default), `|AltMatch` and `|RefMatch` rows are discarded.
#' @param verbose logical. Whether to print progress messages. Default is `TRUE`.
#'
#' @return A named list with two numeric matrices, both with markers as rows
#'   and samples as columns:
#'   \describe{
#'     \item{`ref_matrix`}{Reference allele read counts.}
#'     \item{`size_matrix`}{Total read counts (reference + alternative).}
#'   }
#'
#' @examples
#' # Get the path to the MADC file
#' madc_path <- system.file("iris_DArT_MADC.csv", package = "BIGr")
#'
#' # Extract the read count matrices
#' counts_matrices <- get_countsMADC(madc_path)
#'
#' # Access the reference and size matrices
#' # ref_matrix  <- counts_matrices$ref_matrix
#' # size_matrix <- counts_matrices$size_matrix
#'
#' rm(counts_matrices)
#'
#' @seealso [get_counts()], [check_madc_sanity()]
#'
#' @import dplyr
#' @export
get_countsMADC <- function(madc_file = NULL, madc_object = NULL, collapse_matches_counts = FALSE, verbose = TRUE) {

  # Add check inputs
  if(is.null(madc_file) && is.null(madc_object)) stop("Please provide either madc_file or madc_object.")
  if(!is.null(madc_file) && !file.exists(madc_file)) stop("MADC file not found. Please provide a valid path.")
  if(!is.null(madc_object) && !is.data.frame(madc_object)) stop("madc_object must be a data frame.")

  vmsg(paste0("Extracting read counts from ", ifelse(!is.null(madc_file), paste0("MADC file: ", madc_file), "madc_object")), verbose = verbose, level = 0, type = ">>")
  vmsg(ifelse(collapse_matches_counts,
              "|AltMatch and |RefMatch counts will be collapsed into their respective |Ref and |Alt alleles.",
              "|AltMatch and |RefMatch rows will be discarded (collapse_matches_counts = FALSE)."),
       verbose = verbose, level = 1, type = ">>")

  # This function takes the MADC file as input and generates a Ref and Alt counts dataframe as output
  if (is.null(madc_object)) {
    update_df <- get_counts(madc_file = madc_file, collapse_matches_counts = collapse_matches_counts, verbose = verbose)
  } else {
    update_df <- get_counts(madc_object = madc_object, collapse_matches_counts = collapse_matches_counts, verbose = verbose)
  }
  # Ensure plain data.frame so row.names<- does not trigger tibble deprecation warning
  update_df <- as.data.frame(update_df)

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
    # Find the common CloneIDs between the two dataframes
    all_mks <- unique(c(rownames(ref_df), rownames(alt_df)))
    common_ids <- intersect(rownames(ref_df), rownames(alt_df))
    n_singles <- length(all_mks) - length(common_ids)

    vmsg(paste("There are", n_singles,"Ref tags without corresponding Alt tags, or vice versa"), verbose = verbose, level = 2, type = ">>")
    vmsg("Only the markers with both Ref and Alt tags will be retained for the conversion", verbose = verbose, level = 1, type = ">>")

    warning(paste("There are", n_singles,"Ref tags without corresponding Alt tags, or vice versa. Only the markers with both Ref and Alt tags will be retained for the conversion"))

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
  vmsg(paste0("Percentage of missing data (datapoints with 0 total count): ", round(ratio_missing_data * 100, 2), "%"), verbose = verbose, level = 2, type = ">>")

  # Return the ref and alt matrices as a list
  matrices_list <- list(ref_matrix = ref_matrix, size_matrix = size_matrix)

  return(matrices_list)

}

#' Read and Pre-process a MADC File
#'
#' Reads a DArTag MADC CSV file (or accepts a pre-read data frame), detects the
#' file format, and returns a filtered data frame containing only `Ref` and `Alt`
#' haplotype rows ready for count-matrix construction.
#'
#' @details
#' **Input**: either `madc_file` (path to CSV) or `madc_object` (pre-read data
#' frame) must be supplied; at least one is required.
#'
#' **Format detection** (applied to file or object alike): the first seven rows
#' of the first column are inspected:
#' - **Standard format**: all entries are blank or `"*"` — the first 7 rows are
#'   treated as DArT placeholder rows and skipped.
#' - **Fixed-allele-ID format**: no filler rows — data are used as-is.
#'
#' **`|AltMatch` / `|RefMatch` handling** (controlled by `collapse_matches_counts`):
#' - `FALSE` (default): these rows are simply discarded.
#' - `TRUE`: their counts are summed into the corresponding `|Ref` or `|Alt`
#'   row for the same `CloneID`.
#'
#' In all cases, trailing suffixes on `AlleleID` (e.g., `|Ref_001`, `|Alt_002`)
#' are stripped to the canonical `|Ref` / `|Alt` form.
#'
#' @param madc_file character or `NULL`. Path to the input MADC CSV file.
#'   At least one of `madc_file` or `madc_object` must be provided.
#' @param madc_object data frame or `NULL`. A pre-read MADC data frame
#'   (e.g., from `check_botloci()`). When supplied, file reading is skipped.
#'   At least one of `madc_file` or `madc_object` must be provided.
#' @param collapse_matches_counts logical. If `TRUE`, counts for `|AltMatch`
#'   and `|RefMatch` rows are summed into their corresponding `|Ref` and `|Alt`
#'   rows. If `FALSE` (default), those rows are discarded.
#' @param verbose logical. Whether to print progress messages. Default is `TRUE`.
#'
#' @return A data frame with one row per `Ref` or `Alt` allele entry, retaining
#'   all original columns (`AlleleID`, `CloneID`, `AlleleSequence`, sample
#'   count columns, etc.).
#'
#' @importFrom dplyr mutate group_by summarise across where select
#' @importFrom dplyr %>%
#'
#' @keywords internal
get_counts <- function(madc_file = NULL, madc_object = NULL, collapse_matches_counts = FALSE, verbose = TRUE) {

  # Add check inputs
  if(is.null(madc_file) && is.null(madc_object)) stop("Please provide either madc_file or madc_object.")
  if(!is.null(madc_file) && !file.exists(madc_file)) stop("MADC file not found. Please provide a valid path.")
  if(!is.null(madc_object) && !is.data.frame(madc_object)) stop("madc_object must be a data frame.")

  # Read the MADC file

  if(!is.null(madc_file)){
  #Read only the first column for the first seven rows
  first_seven_rows <- read.csv(madc_file, header = FALSE, nrows = 7, colClasses = c(NA, "NULL"))

  #Check if all entries in the first column are either blank or "*"
  check_entries <- all(first_seven_rows[, 1] %in% c("", "*"))

  } else {
    check_entries <- all(madc_object[1:min(7L, nrow(madc_object)), 1] %in% c("", "*"))
  }

  #Check if the MADC file has the filler rows or is processed from updated fixed allele ID pipeline
  if (check_entries) {
    #Note: This assumes that the first 7 rows are placeholder info from DArT processing
    #Read the madc file
    vmsg("Detected raw MADC format with 7-row header. Reading file while skipping the first 7 rows.", verbose = verbose, level = 1, type = ">>")
    if(!is.null(madc_file)){
      madc_df <- read.csv(madc_file, sep = ',', skip = 7, check.names = FALSE)
    } else {
      madc_df <- madc_object[-(1:7), ]
    }
  } else {

    #Read the madc file
    vmsg("Detected fixed allele IDs MADC format", verbose = verbose, level = 1, type = ">>")
    if(!is.null(madc_file)){
      madc_df <- read.csv(madc_file, sep = ',', check.names = FALSE)
    } else {
      madc_df <- madc_object
    }
  }

  if(collapse_matches_counts){
    filtered_df <- madc_df[order(madc_df$AlleleID),] %>%
      mutate(Type = ifelse(grepl("Alt", AlleleID), "Alt", "Ref")) %>%
      group_by(CloneID, Type) %>%
      summarise(
        AlleleID = paste0(unique(CloneID), "|", unique(Type)),
        AlleleSequence = first(AlleleSequence),
        across(where(is.numeric), sum),
        .groups = "drop"
      ) %>%
      select(AlleleID, CloneID, AlleleSequence, everything(), -Type)

  } else {
    #Retain only the Ref and Alt haplotypes
    filtered_df <- madc_df[!grepl("\\|AltMatch|\\|RefMatch", madc_df$AlleleID), ]
  }

  #Remove extra text after Ref and Alt (_001 or _002)
  filtered_df$AlleleID <- sub("\\|Ref.*", "|Ref", filtered_df$AlleleID)
  filtered_df$AlleleID <- sub("\\|Alt.*", "|Alt", filtered_df$AlleleID)

  return(filtered_df)
}
