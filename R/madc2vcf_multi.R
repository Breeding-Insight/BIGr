#' Convert MADC file to VCF using polyRAD for multiallelic genotyping
#'
#' This function converts a DArTag MADC file to a VCF using the polyRAD package's
#' `readDArTag` and `RADdata2VCF` pipeline. It runs `check_madc_sanity` before
#' loading the data, applies corrections for lowercase sequences and all-NA
#' rows/columns, and sets `n.header.rows` automatically based on whether the
#' MADC file follows the raw DArT format (6 header rows) or the fixed allele ID
#' format (no header rows).
#'
#' @param madc_file character. Path or URL to the input MADC CSV file.
#' @param botloci_file character. Path or URL to the botloci file listing target
#'   IDs designed on the bottom strand.
#' @param outfile character. Path for the output VCF file.
#' @param markers_info character or NULL. Optional path or URL to a CSV file
#'   with marker metadata. Required when CloneIDs do not follow the
#'   \code{Chr_Pos} format; must contain \code{CloneID} (or
#'   \code{BI_markerID}), \code{Chr}, and \code{Pos} columns.
#' @param ploidy integer. Ploidy level of the samples passed to \code{taxaPloidy}.
#'   Default is 2.
#' @param verbose logical. Whether to print progress messages. Default is TRUE.
#'
#' @return Invisible NULL. Writes a VCF file to \code{outfile}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reads the MADC file and runs \code{check_madc_sanity}.
#'   \item Validates the botloci file against MADC CloneIDs using
#'         \code{check_botloci}, fixing any padding mismatches automatically.
#'   \item Converts lowercase bases in AlleleSequence to uppercase if detected.
#'   \item Removes all-NA rows and columns if detected.
#'   \item Writes the corrected data to a temporary file and passes it to
#'         \code{polyRAD::readDArTag}.
#'   \item Estimates overdispersion with \code{polyRAD::TestOverdispersion} and
#'         calls \code{polyRAD::IterateHWE}, then exports the result with
#'         \code{polyRAD::RADdata2VCF}.
#' }
#'
#' @importFrom utils read.csv write.csv read.table
#'
#' @export
madc2vcf_multi <- function(madc_file,
                           botloci_file,
                           outfile,
                           markers_info = NULL,
                           ploidy = 2L,
                           verbose = TRUE) {

  vmsg("Running BIGr madc2vcf_multi", verbose = verbose, level = 0, type = ">>")
  vmsg("madc_file    : %s", verbose = verbose, level = 1, madc_file)
  vmsg("botloci_file : %s", verbose = verbose, level = 1, botloci_file)
  vmsg("markers_info : %s", verbose = verbose, level = 1, if (is.null(markers_info)) "NULL" else markers_info)
  vmsg("outfile      : %s", verbose = verbose, level = 1, outfile)
  vmsg("ploidy       : %s", verbose = verbose, level = 1, ploidy)

  vmsg("Checking inputs", verbose = verbose, level = 0, type = ">>")

  if (!(file.exists(madc_file) | url_exists(madc_file)))       stop("MADC file not found. Please provide a valid path or URL.")
  if (!(file.exists(botloci_file) | url_exists(botloci_file))) stop("Botloci file not found. Please provide a valid path or URL.")
  if (!is.null(markers_info) && !(file.exists(markers_info) | url_exists(markers_info))) stop("markers_info file not found. Please provide a valid path or URL.")
  if (!is.numeric(ploidy) || ploidy < 1) stop("ploidy must be a positive integer.")

  # ---- Load markers_info if provided ----
  mi_df <- if (!is.null(markers_info)) read.csv(markers_info) else NULL

  # ---- Read and sanity-check ----
  report <- read.csv(madc_file, check.names = FALSE)
  checks <- check_madc_sanity(report)

  messages_results <- mapply(function(check, message) {
    if (check) message[1] else message[2]
  }, checks$checks, checks$messages)

  for (i in seq_along(messages_results))
    vmsg(messages_results[i], verbose = verbose, level = 1, type = ">>")

  if (!checks$checks["Columns"])
    stop("The MADC file is missing required columns (CloneID, AlleleID, AlleleSequence)")

  if (checks$checks["IUPACcodes"])
    stop("MADC Allele Sequences contain IUPAC (non-ATCG) codes. Please run HapApp to clean MADC file before using this function.")

  if (!isTRUE(checks$checks["RefAltSeqs"]))
    stop("Not all Ref sequences have a corresponding Alt or vice versa. Please provide a complete MADC file before using this function.")

  if (!isTRUE(checks$checks["FixAlleleIDs"]))
    stop("The MADC file does not have fixed AlleleIDs. Please process the MADC file through HapApp before using this function.")

  if (!isTRUE(checks$checks["ChromPos"])) {
    if (is.null(markers_info))
      stop("CloneID column does not follow the 'Chr_Pos' format. ",
           "Please provide a markers_info file with at least 'CloneID'/'BI_markerID', ",
           "'Chr', and 'Pos' columns.")
    if (!all(c("Chr", "Pos") %in% colnames(mi_df)))
      stop("CloneID column does not follow the 'Chr_Pos' format. ",
           "markers_info must contain at least 'Chr' and 'Pos' columns to remap marker IDs.")
  }

  # ---- Check botloci vs MADC CloneIDs ----
  vmsg("Checking botloci file", verbose = verbose, level = 0, type = ">>")
  cloneids_before <- report$CloneID
  botloci_df      <- read.table(botloci_file, header = FALSE)
  botloci_before  <- botloci_df$V1
  checked_botloci <- check_botloci(botloci_df, report, ChromPos = checks$checks["ChromPos"], mi_df = mi_df, verbose = verbose)
  botloci_df      <- checked_botloci[[1]]
  report          <- checked_botloci[[2]]
  mi_df           <- checked_botloci[[3]]
  cloneid_changed  <- !identical(report$CloneID, cloneids_before)
  botloci_changed  <- !identical(botloci_df$V1, botloci_before)

  # ---- Botloci temp file (if IDs were remapped) ----
  if (botloci_changed) {
    tmp_botloci <- tempfile()
    on.exit(unlink(tmp_botloci), add = TRUE)
    write.table(botloci_df, tmp_botloci, row.names = FALSE, col.names = FALSE, quote = FALSE)
    botloci_input <- tmp_botloci
  } else {
    botloci_input <- botloci_file
  }

  # ---- Corrections: only create a temp file if needed ----
  need_temp <- isTRUE(checks$checks["allNArow"]) || isTRUE(checks$checks["allNAcol"]) || cloneid_changed

  if (need_temp) {
    if (checks$checks["LowerCase"]) {
      vmsg("MADC Allele Sequences contain lowercase characters. Converting to uppercase",
           verbose = verbose, level = 1, type = ">>")
      report$AlleleSequence <- toupper(report$AlleleSequence)
    }

    if (checks$checks["allNArow"]) {
      idx <- apply(report, 1, function(x) all(is.na(x) | x == ""))
      vmsg("Removing %s all-NA row(s)", verbose = verbose, level = 1, type = ">>", sum(idx))
      report <- report[!idx, ]
    }

    if (checks$checks["allNAcol"]) {
      idx <- apply(report, 2, function(x) all(is.na(x) | x == ""))
      vmsg("Removing %s all-NA column(s)", verbose = verbose, level = 1, type = ">>", sum(idx))
      report <- report[, !idx]
    }

    tmp_madc <- tempfile(fileext = ".csv")
    on.exit(unlink(tmp_madc), add = TRUE)
    write.csv(report, tmp_madc, row.names = FALSE, quote = TRUE)
    input_file <- tmp_madc
  } else {
    if (checks$checks["LowerCase"])
      vmsg("MADC Allele Sequences contain lowercase characters. polyRAD will handle them",
           verbose = verbose, level = 1, type = ">>")
    input_file <- madc_file
  }

  vmsg("Loading MADC into polyRAD", verbose = verbose, level = 0, type = ">>")

  raddat <- polyRAD::readDArTag(
    file              = input_file,
    botloci           = botloci_input,
    n.header.rows     = 0L,
    sample.name.row   = 1,
    trim.sample.names = "",
    taxaPloidy        = as.integer(ploidy)
  )

  overdispersionP <- polyRAD::TestOverdispersion(raddat)
  my_ovdisp <- overdispersionP$optimal

  vmsg("Running HWE iteration (overdispersion = %s)", verbose = verbose, level = 0, type = ">>", my_ovdisp)

  raddat_hwe <- polyRAD::IterateHWE(raddat, overdispersion = my_ovdisp)

  vmsg("Writing VCF to %s", verbose = verbose, level = 0, type = ">>", outfile)

  polyRAD::RADdata2VCF(raddat_hwe, file = outfile, asSNPs = FALSE, hindhe = FALSE)

  vmsg("Done!", verbose = verbose, level = 0, type = ">>")

  invisible(NULL)
}

