#' Run basic sanity checks on a MADC-style allele report
#'
#' @description
#' Performs five quick validations on an allele report:
#' 1) **Columns** – required columns are present (`CloneID`, `AlleleID`, `AlleleSequence`);
#' 2) **FixAlleleIDs** – first column’s first up-to-6 rows are not all blank or "*"
#'    *and* both `_0001` and `_0002` appear in `AlleleID`;
#' 3) **IUPACcodes** – presence of non-ATCG characters in `AlleleSequence`;
#' 4) **LowerCase** – presence of lowercase a/t/c/g in `AlleleSequence`;
#' 5) **Indels** – reference/alternate allele lengths differ for the same `CloneID`.
#'
#' @param report A `data.frame` with at least the columns
#'   `CloneID`, `AlleleID`, and `AlleleSequence`. The first column is also
#'   used in the “FixAlleleIDs” check to inspect its first up to six entries.
#'
#' @details
#' - **IUPAC check:** Flags any character outside `ATCG` (case-insensitive),
#'   which will include ambiguity codes (`N`, `R`, `Y`, etc.) and symbols like `-`.
#' - **Indels:** Rows are split by `AlleleID` containing `"Ref_0001"` vs `"Alt_0002"`,
#'   merged by `CloneID`, and the lengths of `AlleleSequence` are compared.
#' - If required columns are missing, only **Columns** is evaluated (`FALSE`) and
#'   `indel_clone_ids` is returned as `NULL`.
#'
#' @return A list with:
#' \describe{
#'   \item{checks}{Named logical vector with entries
#'     `Columns`, `FixAlleleIDs`, `IUPACcodes`, `LowerCase`, `Indels`.}
#'   \item{indel_clone_ids}{Character vector of `CloneID`s where ref/alt lengths differ.
#'     Returns `character(0)` if none, or `NULL` when required columns are missing.}
#' }
#'
#'
#' @export
check_madc_sanity <- function(report) {

  # Initialize
  checks <- c(Columns = NA, FixAlleleIDs = NA, IUPACcodes = NA, LowerCase = NA, Indels = NA, ChromPos = NA)
  messages <-  list(Columns = NA, FixAlleleIDs = NA, IUPACcodes = NA, LowerCase = NA, Indels = NA, ChromPos = NA)

  # Validate required columns
  required <- c("CloneID", "AlleleID", "AlleleSequence")
  missing_cols <- setdiff(required, names(report))
  checks["Columns"] <- length(missing_cols) == 0

  if(checks[["Columns"]]){
    # ---- FixAlleleIDs ----
    # Check if first up-to-6 entries in the *first column* are all "" or "*"
    n <- nrow(report)
    idx <- seq_len(min(6L, n))
    first_col_vals <- report[[1]][idx]
    all_blank_or_star <- all(first_col_vals %in% c("", "*"), na.rm = TRUE)
    # Also require that both _0001 and _0002 appear in AlleleID
    has_0001 <- any(grepl("_0001", report$AlleleID, fixed = TRUE), na.rm = TRUE)
    has_0002 <- any(grepl("_0002", report$AlleleID, fixed = TRUE), na.rm = TRUE)
    checks["FixAlleleIDs"] <- (!all_blank_or_star) & has_0001 & has_0002

    # ---- IUPACcodes ----
    iu <- grepl("[^ATCG]", report$AlleleSequence, ignore.case = TRUE)
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
      indel_mask <- cmp_ok & (ref_len != alt_len)
      checks["Indels"] <- any(indel_mask)
      indels <- if (any(indel_mask)) merged$CloneID[indel_mask] else character(0)
    } else {
      checks["Indels"] <- FALSE
      indels <- character(0)
    }

    # ---- Chrom Pos ----
    pos <- strsplit(report[,2], "_")
    checks["ChromPos"] <- all(sapply(pos, length) == 2)

  } else indels <- NULL

  messages[["Columns"]] <- c("Required columns are present\n",
                           "One or more required columns missing. Verify if your file has columns: CloneID, AlleleID, AlleleSequence\n")
  messages[["FixAlleleIDs"]] <- c("Fixed Allele IDs look good\n",
                               "MADC not processed by BI. Please contact us to assign allele IDs to your MADC according to the specie haplotype dabatase. This guarantee reproducibility between diferent datasets\n")
  messages[["IUPACcodes"]] <- c("IUPAC (non-ATCG) codes found in AlleleSequence. This codes are not currently supported\n",
                             "No IUPAC (non-ATCG) codes found in AlleleSequence\n")
  messages[["LowerCase"]] <- c("Lowercase bases found in AlleleSequence\n",
                            "No lowercase bases found in AlleleSequence\n")
  messages[["Indels"]] <- c(paste("Indels found (ref/alt lengths differ) for the CloneIDs:",paste(indels, collapse = " ")),
                         "No indels found (ref/alt lengths match) for all CloneIDs\n")
  messages[["ChromPos"]] <- c("Chromosome and Position format in CloneID look good\n",
                             "CloneID does not have the expected Chromosome_Position format. Please check your CloneIDs or provide a file with this information\n")

  list(checks = checks, messages = messages, indel_clone_ids = indels)
}
