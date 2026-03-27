#' Run basic sanity checks on a MADC-style allele report
#'
#' @description
#' Performs nine quick validations on an allele report:
#' 1) **Columns** - required columns are present (`CloneID`, `AlleleID`, `AlleleSequence`);
#' 2) **FixAlleleIDs** - first column's first up-to-6 rows are not all blank or `"*"`
#'    *and* both `_0001` and `_0002` appear in `AlleleID`;
#' 3) **IUPACcodes** - presence of non-ATCG characters in `AlleleSequence`;
#' 4) **LowerCase** - presence of lowercase a/t/c/g in `AlleleSequence`;
#' 5) **Indels** - reference/alternate allele lengths differ for the same `CloneID`,
#'    or a `"-"` character is present in `AlleleSequence`;
#' 6) **ChromPos** - all `CloneID` values follow the `Chr_Position` format
#'    (prefix matches `"chr"` case-insensitively, suffix is a positive integer);
#' 7) **allNAcol** - at least one column contains only `NA` or empty values;
#' 8) **allNArow** - at least one row contains only `NA` or empty values;
#' 9) **RefAltSeqs** - every `CloneID` has at least one `Ref` and one `Alt` allele row.
#'
#' @param report A `data.frame` with at least the columns
#'   `CloneID`, `AlleleID`, and `AlleleSequence`. The first column is also
#'   used in the `FixAlleleIDs` check to inspect its first up to six entries.
#'   If `FixAlleleIDs` is `FALSE` (raw DArT format), the first 7 rows are
#'   treated as header filler and skipped before further checks are run.
#'
#' @details
#' - **FixAlleleIDs:** When the first six rows of the first column are all blank
#'   or `"*"` (raw DArT format), row 7 is promoted to column headers and the
#'   first 7 rows are dropped before subsequent checks are run. The check is
#'   `TRUE` when the file has already been processed by HapApp (fixed IDs with
#'   `_0001`/`_0002` suffixes).
#' - **IUPAC check:** Flags any character outside `A`, `T`, `C`, `G` and `"-"`
#'   (case-insensitive), which includes ambiguity codes (`N`, `R`, `Y`, etc.).
#' - **Indels:** Rows are split by `AlleleID` containing `"Ref_0001"` vs
#'   `"Alt_0002"`, merged by `CloneID`, and flagged as indels if either (a) the
#'   lengths of `AlleleSequence` differ, (b) the sequences have the same length
#'   but more than one character differs between them (complex indel / local
#'   rearrangement), or (c) a `"-"` character is present anywhere in
#'   `AlleleSequence`.
#' - **ChromPos:** Each `CloneID` is split on `"_"` into exactly two parts; the
#'   first part must match `"Chr"` (case-insensitive) and the second must be a
#'   positive integer. Returns `FALSE` when any `CloneID` is `NA`.
#' - **allNAcol / allNArow:** Detected via `apply()` over columns/rows
#'   respectively; a cell is treated as missing when it is `NA` or an empty
#'   string (`""`). Useful for flagging empty or corrupt entries.
#' - **RefAltSeqs:** For each unique `CloneID`, checks whether at least one row
#'   with a `Ref` (`|Ref_` when `FixAlleleIDs = TRUE`, `|Ref$` otherwise) and
#'   one row with an `Alt` (`|Alt_` / `|Alt$`) allele exist. `CloneID`s that
#'   lack a `Ref` row are stored in `missRef`; those lacking an `Alt` row in
#'   `missAlt`. The check is `TRUE` when both sets are empty.
#' - If required columns are missing (`Columns = FALSE`), only `Columns` and
#'   `FixAlleleIDs` are evaluated; all other checks remain `NA` and
#'   `indel_clone_ids`, `missRef`, and `missAlt` are returned as `NULL`.
#'
#' @return A named list with five elements:
#' \describe{
#'   \item{checks}{Named logical vector with nine entries:
#'     `Columns`, `FixAlleleIDs`, `IUPACcodes`, `LowerCase`, `Indels`,
#'     `ChromPos`, `allNAcol`, `allNArow`, `RefAltSeqs`.
#'     `TRUE` means the condition was detected (or passed for `Columns`,
#'     `FixAlleleIDs`, `ChromPos`, and `RefAltSeqs`); `NA` means the check
#'     was skipped.}
#'   \item{messages}{Named list of length-2 character vectors, one per check.
#'     Element `[[1]]` is the message when the check is `TRUE`, element `[[2]]`
#'     when it is `FALSE`. Indexed by the same names as `checks`.}
#'   \item{indel_clone_ids}{Character vector of `CloneID`s where ref/alt
#'     lengths differ. Returns `character(0)` if none are found, or `NULL`
#'     when required columns are missing.}
#'   \item{missRef}{Character vector of `CloneID`s that have no `Ref` allele
#'     row. Returns `character(0)` if all `CloneID`s have a `Ref` row, or
#'     `NULL` when required columns are missing.}
#'   \item{missAlt}{Character vector of `CloneID`s that have no `Alt` allele
#'     row. Returns `character(0)` if all `CloneID`s have an `Alt` row, or
#'     `NULL` when required columns are missing.}
#' }
#'
#' @export
check_madc_sanity <- function(report) {

  # Initialize
  checks <- c(Columns = NA, FixAlleleIDs = NA, IUPACcodes = NA, LowerCase = NA, Indels = NA,
              ChromPos = NA, allNAcol = NA, allNArow = NA, RefAltSeqs = NA, OtherAlleles = NA)
  messages <-  list(Columns = NA, FixAlleleIDs = NA, IUPACcodes = NA, LowerCase = NA, Indels = NA,
                    ChromPos = NA, allNAcol = NA, allNArow = NA, RefAltSeqs = NA, OtherAlleles = NA)

  # ---- FixAlleleIDs ----
  # Check if first up-to-6 entries in the *first column* are all "" or "*"
  n <- nrow(report)
  idx <- seq_len(min(6L, n))
  first_col_vals <- report[[1]][idx]
  all_blank_or_star <- all(first_col_vals %in% c("", "*"), na.rm = TRUE)
  # Also require that both _0001 and _0002 appear in AlleleID
  if(all_blank_or_star) {
    colnames(report) <- report[7,]
    report <- report[-c(1:7),]
  }
  has_0001 <- any(grepl("_0001", report$AlleleID, fixed = TRUE), na.rm = TRUE)
  has_0002 <- any(grepl("_0002", report$AlleleID, fixed = TRUE), na.rm = TRUE)
  checks["FixAlleleIDs"] <- (!all_blank_or_star) & has_0001 & has_0002

  # Validate required columns
  required <- c("CloneID", "AlleleID", "AlleleSequence")
  missing_cols <- setdiff(required, names(report))
  checks["Columns"] <- length(missing_cols) == 0

  if(checks[["Columns"]]){
    # ---- IUPACcodes ----
    iu <- grepl("[^ATCG-]", report$AlleleSequence, ignore.case = TRUE)
    checks["IUPACcodes"] <- any(iu, na.rm = TRUE)

    # ---- LowerCase ----
    lc <- grepl("[atcg]", report$AlleleSequence)
    checks["LowerCase"] <- any(lc, na.rm = TRUE)

    # ---- Indels ----
    refs <- subset(report, grepl("Ref_0001", AlleleID, fixed = TRUE),
                   select = c(CloneID, AlleleID, AlleleSequence))
    alts <- subset(report, grepl("Alt_0002", AlleleID, fixed = TRUE),
                   select = c(CloneID, AlleleID, AlleleSequence))

    merged <- merge(refs, alts, by = "CloneID", suffixes = c("_ref", "_alt"), all = FALSE)

    if (nrow(merged) > 0) {
      ref_len <- nchar(merged$AlleleSequence_ref, keepNA = TRUE)
      alt_len <- nchar(merged$AlleleSequence_alt, keepNA = TRUE)
      cmp_ok <- !is.na(ref_len) & !is.na(alt_len)

      # Classic indel: different lengths
      indel_mask <- cmp_ok & (ref_len != alt_len)

      # Complex indel: same length but >1 character difference between sequences
      same_len <- cmp_ok & (ref_len == alt_len)
      if (any(same_len)) {
        n_diffs <- mapply(function(r, a) {
          r_chars <- strsplit(r, "")[[1]]
          a_chars <- strsplit(a, "")[[1]]
          sum(r_chars != a_chars)
        }, merged$AlleleSequence_ref[same_len], merged$AlleleSequence_alt[same_len])
        indel_mask[same_len] <- n_diffs > 1
      }

      checks["Indels"] <- any(indel_mask) | any(grepl("-", report$AlleleSequence), na.rm = TRUE)
      indels <- if (any(indel_mask)) merged$CloneID[indel_mask] else character(0)

    } else {
      checks["Indels"] <- FALSE
      indels <- character(0)
    }

    # --- All NA ----
    checks["allNArow"] <- any(apply(report, 1, function(x) all(is.na(x) | x == "")))
    checks["allNAcol"] <- any(apply(report, 2, function(x) all(is.na(x) | x == "")))

    # ---- Chrom Pos ----
    if(!any(is.na(report$CloneID))) {
      pos <- strsplit(report$CloneID, "_")
      format <- all(sapply(pos, length) == 2)
      first <- all(grepl("^[A-Za-z]", sapply(pos, "[", 1)))
      second <- suppressWarnings(all(sapply(pos, function(x) as.numeric(x[2])) > 0))
      checks["ChromPos"] <- all(format, first, second)
    } else checks["ChromPos"] <- FALSE

    # ---- RefAltSeqs ----
    all_clones <- unique(report$CloneID)
    if (isTRUE(checks["FixAlleleIDs"])) {
      has_ref <- unique(report$CloneID[grepl("\\|Ref_",  report$AlleleID)])
      has_alt <- unique(report$CloneID[grepl("\\|Alt_",  report$AlleleID)])
    } else {
      has_ref <- unique(report$CloneID[grepl("\\|Ref$",  report$AlleleID)])
      has_alt <- unique(report$CloneID[grepl("\\|Alt$",  report$AlleleID)])
    }
    missRef <- setdiff(all_clones, has_ref)
    missAlt <- setdiff(all_clones, has_alt)
    checks["RefAltSeqs"] <- length(missRef) == 0 & length(missAlt) == 0

    # ---- OtherAlleles ----
    checks["OtherAlleles"] <- any(grepl("[|]Other", report$AlleleID))

  } else {
    indels  <- NULL
    missRef <- NULL
    missAlt <- NULL
  }

  messages[["Columns"]] <- c("Required columns are present",
                             "One or more required columns missing. Verify if your file has columns: CloneID, AlleleID, AlleleSequence")
  messages[["FixAlleleIDs"]] <- c("Fixed Allele IDs look good",
                                  "MADC not processed by HapApp.")
  messages[["IUPACcodes"]] <- c("IUPAC (non-ATCG) codes found in AlleleSequence. This codes are not currently supported",
                                "No IUPAC (non-ATCG) codes found in AlleleSequence")
  messages[["LowerCase"]] <- c("Lowercase bases found in AlleleSequence",
                               "No lowercase bases found in AlleleSequence")
  messages[["Indels"]] <- c(paste("Indels found (ref/alt lengths differ or >1 mismatch between same-length sequences) for the CloneIDs:",paste(indels, collapse = " ")),
                            "No indels found (ref/alt lengths match and at most 1 mismatch) for all CloneIDs")
  messages[["ChromPos"]] <- c("Chromosome and Position format in CloneID look good",
                              "CloneID does not have the expected Chromosome_Position format. Please check your CloneIDs or provide a file with this information")
  messages[["allNArow"]] <- c("One or more rows contain all NA values.",
                              "No rows with all NA values")
  messages[["allNAcol"]] <- c("One or more columns contain all NA values.",
                              "No columns with all NA values")
  messages[["RefAltSeqs"]] <- c("All CloneIDs have both Ref and Alt alleles",
                                paste0("Some CloneIDs are missing Ref and/or Alt alleles. ",
                                       "Missing Ref: ", paste(missRef, collapse = " "), ". ",
                                       "Missing Alt: ", paste(missAlt, collapse = " "), "."))
  messages[["OtherAlleles"]] <- c("Alleles other than Ref and Alt were found in AlleleID.",
                                  "No alleles other than Ref and Alt found in AlleleID.")

  list(checks = checks, messages = messages, indel_clone_ids = indels,
       missRef = missRef, missAlt = missAlt)
}

#' Check and Adjust Botloci and MADC Marker Compatibility
#'
#' This internal function checks the compatibility between botloci and MADC markers. It ensures that the marker IDs in the botloci file match those in the MADC file. If discrepancies are found, such as mismatched padding, the function attempts to adjust the IDs to ensure compatibility.
#'
#' @param botloci A data frame containing the botloci markers.
#' @param report A data frame containing the MADC markers.
#' @param ChromPos logical value indicating whether the CloneID in the MADC file contains chromosome and position information in the format "Chr_Pos". Default is TRUE
#' @param mi_df A data frame containing marker information with columns CloneID, Chr, and Pos. Required if `ChromPos` is FALSE.
#' @param verbose A logical value indicating whether to print detailed messages about the adjustments. Default is TRUE. Required if `ChromPos` is FALSE.
#'
#' @return A list containing the adjusted botloci and MADC data frames.
#'
#' @details
#' The function checks if the marker IDs in the botloci file are present in the MADC file. If no matches are found, it examines the padding (number of digits) in the marker IDs and adjusts them to match the longest padding. If the IDs still do not match after adjustment, an error is raised. This function is intended for internal use and helps ensure that the botloci and MADC files are compatible for downstream analysis.
#'
#' @keywords internal
#' @noRd
check_botloci <- function(botloci, report, ChromPos=TRUE, mi_df = NULL, verbose=TRUE){

  # Check inputs
  if(!ChromPos) {
    if(is.null(mi_df)) stop("When MADC CloneID don't follow the format Chr_Pos, a marker_info file with CloneID, Chr and Pos columns must be provided.")
    # if exists, it must contain CloneID or BI_markerID that matches the report$CloneID, and Chr and Pos columns
    if(!any(mi_df$CloneID %in% report$CloneID) & !any(mi_df$BI_markerID %in% report$CloneID)) {
      stop("None of the MADC CloneID could be found in the markers_info CloneID or BI_markerID. Please make sure they match.")
    } else {
      use_col <- if(any(mi_df$CloneID %in% report$CloneID)) "CloneID" else "BI_markerID"
      vmsg(paste("Using", use_col, "column in marker_info to match MADC CloneID"), verbose = verbose, level = 1, type = ">>")
    }
    if(is.null(mi_df$Chr) | is.null(mi_df$Pos)) stop("When MADC CloneID don't follow the format Chr_Pos, Chr and Pos columns must be provided in the markers_info file.")
  }

  if(!any(botloci$V1 %in% report$CloneID)) { # First check if any botloci markers are found in MADC file. If not, check for padding mismatch.
    vmsg("No botloci markers found in MADC file. Checking for padding mismatch...", verbose = verbose, level = 1, type = ">>")

    pad_madc <- unique(nchar(sub(".*_", "", report$CloneID)))
    pad_botloci <- unique(nchar(sub(".*_", "", botloci$V1)))

    if(length(pad_madc) > 1 | length(pad_botloci) > 1) stop("Check marker IDs in both MADC and botloci files. They should be the same.")

    if(pad_madc != pad_botloci) {
      vmsg("Padding between MADC and botloci files do not match. Markers ID modified to match longest padding.", verbose = verbose, level = 1, type = ">>")
      if (pad_madc < pad_botloci) {
        report$CloneID <- paste0(sub("_(.*)", "", report$CloneID), "_",
                                 sprintf(paste0("%0", pad_botloci, "d"), as.integer(sub(".*_", "", report$CloneID)))
        )
        report$AlleleID <- paste0(report$CloneID, "|", sapply(strsplit(report$AlleleID, "[|]"), "[[",2))
      } else {
        botloci$V1 <- paste0(sub("_(.*)", "", botloci$V1), "_",
                             sprintf(paste0("%0", pad_madc, "d"), as.integer(sub(".*_", "", botloci$V1)))
        )
        if(!any(botloci$V1 %in% report$CloneID)) stop("After matching padding, botloci markers still not found in MADC file. Check marker IDs.\n")
      }
    } else if (!(is.null(mi_df$Chr) | is.null(mi_df$Pos))){
      vmsg("It is not a padding mismatch issue.", verbose = verbose, level = 1, type = ">>")
      vmsg("Checking if jointing provided Chromosome and Position information in marker_file solve the issue", verbose = verbose, level = 1, type = ">>")
      if(!any(mi_df$CloneID %in% report$CloneID) & !any(mi_df$BI_markerID %in% report$CloneID)) {
        stop("None of the MADC CloneID could be found in the markers_info CloneID or BI_markerID. Please make sure they match.")
      } else {
        use_col <- if(any(mi_df$CloneID %in% report$CloneID)) "CloneID" else "BI_markerID"
        vmsg(paste("Using", use_col, "column in marker_info to match MADC CloneID"), verbose = verbose, level = 1, type = ">>")
      }
      mk_info_CloneID <- paste0(mi_df$Chr, "_", sprintf(paste0("%0",pad_botloci, "d"), as.integer(mi_df$Pos)))

      if(!any(botloci$V1 %in% mk_info_CloneID)){
        vmsg("It is not a padding mismatch issue.", verbose = verbose, level = 1, type = ">>")
        vmsg("Chromosome and Position information in marker_file don't solve the issue.", verbose = verbose, level = 1, type = ">>")
        stop("Check marker IDs in both MADC and botloci files. They should be the same.")
      } else {
        vmsg("Chromosome and Position information in marker_file solve the issue.", verbose = verbose, level = 1, type = ">>")
        vmsg("Using this information to modify MADC CloneIDs to match botloci markers.", verbose = verbose, level = 1, type = ">>")
        report$CloneID <- mk_info_CloneID[match(report$CloneID, mi_df[[use_col]])]
      }
    } else {
      vmsg("It is not a padding mismatch issue.", verbose = verbose, level = 1, type = ">>")
      vmsg("Chromosome and Position information in marker_file not provided.", verbose = verbose, level = 1, type = ">>")
      stop("Check marker IDs in both MADC and botloci files. They should be the same.")
    }
  }
  return(list(botloci, report))
}
