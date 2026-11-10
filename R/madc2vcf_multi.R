#' Convert MADC file to VCF using polyRAD for multiallelic genotyping
#'
#' This function converts a DArTag fixed allele ID MADC file to a VCF
#' containing multiallelic markers based on the microhaplotypes using
#' the polyRAD package's `readDArTag`, `IterateHWE` population model
#' and `RADdata2VCF` pipeline.
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
#' @param vcf_header character or NULL. Optional path or URL to a VCF header file
#'   with the \code{##reference} and \code{##contig} lines of the reference
#'   genome, such as the VCF header files provided with the microhaplotype
#'   databases. Its \code{##contig} lines replace the ones written by polyRAD,
#'   which have no contig lengths. The other lines of the file are not used. A
#'   warning is given if a chromosome in the output is not listed in the
#'   \code{##contig} lines.
#' @param reference character or NULL. Optional reference genome ID written to
#'   the \code{##reference} line (e.g. \code{"GCF_020283565.1"}). Use it when no
#'   \code{vcf_header} file is available, or to replace the \code{##reference}
#'   line of \code{vcf_header}. If neither \code{vcf_header} nor
#'   \code{reference} is given, a warning is given and the header has no
#'   \code{##reference} line and only the \code{##contig} lines written by polyRAD.
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
#'   \item Adds the reference genome lines from \code{vcf_header} and/or
#'         \code{reference} and a \code{##BIGrCommandLine} line to the VCF header.
#' }
#'
#' **Sanity check behaviour and requirements**
#'
#' The function always stops if IUPAC codes, unpaired Ref/Alt sequences, or
#' unfixed AlleleIDs are detected (see \code{check_madc_sanity}). For the
#' remaining checks the required inputs are:
#'
#' | Check | Status | Required |
#' |---|---|---|
#' | **Indels** | detected | `botloci_file` |
#' | | not detected | `botloci_file` |
#' | **ChromPos** | valid | `botloci_file` |
#' | | invalid | `markers_info` with `Chr`/`Pos` + `botloci_file` |
#'
#' @importFrom utils read.csv write.csv read.table
#'
#' @export
madc2vcf_multi <- function(madc_file,
                           botloci_file,
                           outfile,
                           markers_info = NULL,
                           ploidy = 2L,
                           vcf_header = NULL,
                           reference = NULL,
                           verbose = TRUE) {

  vmsg("Running BIGr madc2vcf_multi", verbose = verbose, level = 0, type = ">>")
  vmsg("madc_file    : %s", verbose = verbose, level = 1, type = ">>", madc_file)
  vmsg("botloci_file : %s", verbose = verbose, level = 1, type = ">>", botloci_file)
  vmsg("markers_info : %s", verbose = verbose, level = 1, type = ">>", if (is.null(markers_info)) "NULL" else markers_info)
  vmsg("outfile      : %s", verbose = verbose, level = 1, type = ">>", outfile)
  vmsg("ploidy       : %s", verbose = verbose, level = 1, type = ">>", ploidy)
  vmsg("vcf_header   : %s", verbose = verbose, level = 1, type = ">>", if (is.null(vcf_header)) "NULL" else vcf_header)
  vmsg("reference    : %s", verbose = verbose, level = 1, type = ">>", if (is.null(reference)) "NULL" else reference)

  vmsg("Checking inputs", verbose = verbose, level = 0, type = ">>")

  if (!(file.exists(madc_file) | url_exists(madc_file)))       stop("MADC file not found. Please provide a valid path or URL.")
  if (!(file.exists(botloci_file) | url_exists(botloci_file))) stop("Botloci file not found. Please provide a valid path or URL.")
  if (!is.null(markers_info) && !(file.exists(markers_info) | url_exists(markers_info))) stop("markers_info file not found. Please provide a valid path or URL.")
  if (!is.numeric(ploidy) || ploidy < 1) stop("ploidy must be a positive integer.")

  # Create a VCF header line with metadata about the command and its parameters
  bigr_meta <- paste0('##BIGrCommandLine.madc2vcf_multi=<ID=madc2vcf_multi,Version="',
                      packageVersion("BIGr"), '",Data="',
                      Sys.time(),'",CommandLine="', vcf_quote_escape(paste0('> madc2vcf_multi(',deparse1(substitute(madc_file)),', ',
                      "botloci_file= ", botloci_file, ', ',
                      "outfile= ", outfile, ', ',
                      "markers_info= ", markers_info, ', ',
                      "ploidy= ", ploidy, ', ',
                      "vcf_header= ", vcf_header, ', ',
                      "reference= ", reference, ', ',
                      "verbose= ", verbose,')')), '">')

  # Reference genome lines for the VCF header
  ref_lines <- vcf_reference_lines(vcf_header, reference)

  # ---- Load markers_info if provided ----
  mi_df <- if (!is.null(markers_info)) read.csv(markers_info) else NULL

  # ---- Read and sanity-check ----
  report <- read.csv(madc_file, check.names = FALSE)
  checks <- check_madc_sanity(report)

  messages_results <- mapply(function(check, message) {
    if (isTRUE(check)) message[1] else message[2]
  }, checks$checks, checks$messages)

  for (i in seq_along(messages_results))
    vmsg(messages_results[i], verbose = verbose, level = 1, type = ">>")

  if (!checks$checks["Columns"])
    stop("The MADC file is missing required columns (CloneID, AlleleID, AlleleSequence)")

  if (checks$checks["IUPACcodes"])
    stop("MADC Allele Sequences contain different IUPAC (non-ATCG) codes in REF and ALT sequences. Please run HapApp to clean MADC file before using this function.")

  if (!isTRUE(checks$checks["RefAltSeqs"]))
    stop("Not all Ref sequences have a corresponding Alt or vice versa. Please provide a complete MADC file before using this function.")

  if (!isTRUE(checks$checks["FixAlleleIDs"]))
    stop("The MADC file does not have fixed AlleleIDs. Please process the MADC file through HapApp before using this function.")

  if (!isTRUE(checks$checks["ChromPos"])) {
    if (is.null(markers_info))
      stop("CloneID column does not follow the 'Chr_Pos' format. ",
           "Please provide a markers_info file with at least 'CloneID'/'Marker_ID'/'BI_markerID', ",
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

  if (!requireNamespace("polyRAD", quietly = TRUE)) {
    stop("Package 'polyRAD' is required for madc2vcf_multi(). Please install it with install.packages('polyRAD').", call. = FALSE)
  }
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

  # Add the reference genome lines and the BIGr command line to the header
  # written by polyRAD
  meta <- update_vcf_reference_lines(read_vcf_meta(outfile), ref_lines)
  replace_vcf_meta(outfile, c(meta, bigr_meta))
  check_vcf_contigs(raddat_hwe$locTable$Chr, ref_lines$contigs)

  vmsg("Done!", verbose = verbose, level = 0, type = ">>")

  invisible(NULL)
}

#' Read the meta lines (starting with ##) of a VCF file
#'
#' @param file Path to an uncompressed VCF file.
#' @return Character vector with the meta lines.
#' @noRd
read_vcf_meta <- function(file) {
  con <- file(file, "r")
  on.exit(close(con))
  meta <- character(0)
  repeat {
    lines <- readLines(con, n = 1000L)
    if (length(lines) == 0) break
    is_meta <- startsWith(lines, "##")
    meta <- c(meta, lines[is_meta])
    if (!all(is_meta)) break
  }
  meta
}

#' Replace the meta lines (starting with ##) of a VCF file
#'
#' The #CHROM line and the records are copied unchanged, in chunks, so the file
#' is not read into memory at once.
#'
#' @param file Path to an uncompressed VCF file.
#' @param meta Character vector with the new meta lines.
#' @return No return value, called for side effects.
#' @noRd
replace_vcf_meta <- function(file, meta) {
  tmp <- tempfile(fileext = ".vcf")
  con_in <- file(file, "r")
  con_out <- file(tmp, "w")
  writeLines(meta, con_out)
  repeat {
    lines <- readLines(con_in, n = 10000L)
    if (length(lines) == 0) break
    lines <- lines[!startsWith(lines, "##")]
    if (length(lines) > 0) writeLines(lines, con_out)
  }
  close(con_in)
  close(con_out)
  file.copy(tmp, file, overwrite = TRUE)
  unlink(tmp)
  invisible()
}
