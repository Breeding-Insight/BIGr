#Internal Functions

utils::globalVariables(c(
  "ALT", "AlleleID", "AlleleSequence", "CHROM", "Concordance", "Data", "ID",
  "MarkerName", "POS",
  "REF", "Type", "Var1", "Variant", "geno",
  "ind", "ref", "row_name", "size", "snp",
  "CloneID", "Count", "qualifying_sites_count",
  "MarkerID", "SampleID", "Dosage",
  "pos", "alt", "match_key",
  "chr",
  # madc_plot() aes symbols
  "n_mhaps", "depth", "len", "marker", "sample", "value",
  "missing_rate", "group", "PCx", "PCy", "Var2",
  "gx", "x0", "x1", "shape", "alt_ratio", "count",
  "ratio_lo", "ratio_hi", "grp"
))

#' Convert GT format to numeric dosage
#' @param gt a genotype matrix with samples as columns and variants as rows
#' @return numeric genotype values
#' @noRd
convert_to_dosage <- function(gt) {
  # Split the genotype string
  alleles <- strsplit(gt, "[|/]")
  # Sum the alleles, treating NA values appropriately
  sapply(alleles, function(x) {
    if (any(is.na(x))) {
      return(NA)
    } else {
      return(sum(as.numeric(x), na.rm = TRUE))
    }
  })
}

#' Verbose Message Utility
#'
#' Prints a formatted verbose message with timestamp, indentation, and type label, if verbose is TRUE.
#'
#' @param text Character string, the message to print (supports sprintf formatting).
#' @param verbose Logical. If TRUE, prints the message; if FALSE, suppresses output.
#' @param level Integer, indentation level (0=header, 1=main step, 2=detail, 3=sub-detail).
#' @param type Character string, message type (e.g., "INFO", "WARN", "ERROR"). Only shown for level 0.
#' @param ... Additional arguments passed to sprintf for formatting.
#'
#' @details Use the verbose argument to control message output. Typically, pass the function's verbose parameter to vmsg.
#'
#' @return No return value, called for side effects.
#' @noRd
vmsg <- function(text, verbose = FALSE, level = 1, type = ">>", ...) {
  if (!verbose) return(invisible())
  # Format timestamp
  timestamp <- format(Sys.time(), "[%H:%M:%S]")

  # Create indentation based on level
  indent <- switch(as.character(level),
                   "0" = "",           # Section headers
                   "1" = "  \u2219 ",       # Main steps (medium bullet)
                   "2" = "    - ",     # Details
                   "3" = "      > ",   # Sub-details
                   paste0(paste(rep("  ", level), collapse = ""), "\u2022 ")  # Fallback for level > 3
  )

  # Format type label (only show for level 0)
  type_label <- if (level == 0) sprintf("%-1s ", type) else ""

  # Format message text
  dots <- list(...)
  if(length(dots) == 0) {
    msg_text <- text
  } else {
    msg_text <- sprintf(text, ...)
  }
  # Combine everything
  formatted_msg <- sprintf("%s %s%s%s", timestamp, type_label, indent, msg_text)
  message(formatted_msg)
}


#' Check Whether a URL Is Accessible
#'
#' Attempts to open a connection to the given URL and returns `TRUE` if
#' successful, `FALSE` otherwise. Errors and warnings are both treated as
#' inaccessible.
#'
#' @param u character. The URL to test.
#'
#' @return A single logical: `TRUE` if the URL can be opened, `FALSE` if not.
#'
#' @keywords internal
#' @noRd
#'
url_exists <- function(u) {
  tryCatch({
    con <- url(u, open = "rb")
    close(con)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
}

#' Build the reference genome lines of a VCF header
#'
#' Returns the ##reference and ##contig lines to write in a VCF header, using
#' the most informative input available:
#' 1. `vcf_header`: a VCF header file with the ##reference and ##contig lines
#'    (reference genome and scaffold lengths). Only these lines are taken from
#'    the file; the calling function writes its own ##fileformat, ##INFO and
#'    ##FORMAT lines.
#' 2. `reference`: the reference genome ID only, written as the ##reference line.
#' 3. Nothing: no ##reference or ##contig lines are written, since placeholder
#'    lines such as "##reference=NA" break tools that parse the header (e.g.,
#'    GIGWA). A warning is given when `warn` is TRUE.
#'
#' @param vcf_header Path or URL to a VCF header file, or NULL.
#' @param reference Reference genome ID, or NULL. If given together with
#'   `vcf_header`, it is used for the ##reference line.
#' @param warn Logical. Warn when neither `vcf_header` nor `reference` is given.
#'
#' @return A list with `lines` (character vector of header lines, possibly
#'   empty) and `contigs` (contig IDs from the ##contig lines, or NULL when no
#'   ##contig lines were provided).
#' @noRd
vcf_reference_lines <- function(vcf_header = NULL, reference = NULL, warn = TRUE) {

  if (!is.null(reference) && (!is.character(reference) || length(reference) != 1L ||
                              is.na(reference) || !nzchar(trimws(reference)))) {
    stop("reference must be a single character string with the reference genome ID.", call. = FALSE)
  }

  ref_line <- if (!is.null(reference)) paste0("##reference=", trimws(reference)) else character(0)
  contig_lines <- character(0)

  if (!is.null(vcf_header)) {
    if (!is.character(vcf_header) || length(vcf_header) != 1L) {
      stop("vcf_header must be the path or URL to a VCF header file.", call. = FALSE)
    }
    if (!(file.exists(vcf_header) | url_exists(vcf_header))) {
      stop("vcf_header file not found. Please provide a valid path or URL.", call. = FALSE)
    }
    header <- trimws(readLines(vcf_header, warn = FALSE), which = "right")

    file_ref <- grep("^##reference=", header, value = TRUE)
    contig_lines <- grep("^##contig=", header, value = TRUE)

    if (length(file_ref) == 0 && length(contig_lines) == 0) {
      stop("vcf_header has no ##reference or ##contig lines: ", vcf_header, call. = FALSE)
    }
    if (length(file_ref) > 1) {
      stop("vcf_header has more than one ##reference line: ", vcf_header, call. = FALSE)
    }
    bad_contig <- !grepl("^##contig=<ID=[^,>]+[,>]", contig_lines)
    if (any(bad_contig)) {
      stop("vcf_header has ##contig lines without an ID, e.g. ", contig_lines[bad_contig][1],
           ". Expected ##contig=<ID=...,length=...>.", call. = FALSE)
    }
    if (length(ref_line) == 0) ref_line <- file_ref
  }

  if (warn && is.null(vcf_header) && is.null(reference)) {
    warning("No reference genome information was provided (vcf_header or reference). ",
            "The ##reference and ##contig lines are left out of the VCF header.", call. = FALSE)
  }

  contigs <- if (length(contig_lines) > 0) sub("^##contig=<ID=([^,>]+).*", "\\1", contig_lines) else NULL

  list(lines = c(ref_line, contig_lines), contigs = contigs)
}

#' Escape a value for a quoted field of a VCF header line
#'
#' Double quotes and backslashes are escaped with a backslash, as required by the
#' VCF specification for quoted values (e.g., the CommandLine field).
#'
#' @param x Character string.
#' @return The escaped string.
#' @noRd
vcf_quote_escape <- function(x) {
  gsub('(["\\\\])', '\\\\\\1', x)
}

#' Replace the reference genome lines of an existing VCF header
#'
#' The ##reference and ##contig lines from `vcf_reference_lines()` replace the
#' ones in `meta` and are placed after the ##fileformat line. Lines of `meta` that
#' are not replaced are kept.
#'
#' @param meta Meta lines (starting with ##) of a VCF header.
#' @param ref_lines List returned by `vcf_reference_lines()`.
#'
#' @return Character vector with the updated meta lines.
#' @noRd
update_vcf_reference_lines <- function(meta, ref_lines) {
  new_ref <- grep("^##reference=", ref_lines$lines, value = TRUE)
  new_contig <- grep("^##contig=", ref_lines$lines, value = TRUE)
  if (length(new_ref) > 0) meta <- meta[!grepl("^##reference=", meta)]
  if (length(new_contig) > 0) meta <- meta[!grepl("^##contig=", meta)]
  fileformat_pos <- grep("^##fileformat=", meta)[1]
  if (is.na(fileformat_pos)) fileformat_pos <- 0
  append(meta, c(new_ref, new_contig), after = fileformat_pos)
}

#' Warn when VCF chromosomes are missing from the ##contig lines
#'
#' @param chrom CHROM values of the VCF records.
#' @param contigs Contig IDs from `vcf_reference_lines()`. Nothing is checked when NULL.
#'
#' @return No return value, called for side effects.
#' @noRd
check_vcf_contigs <- function(chrom, contigs) {
  if (is.null(contigs)) return(invisible())
  missing_chrom <- setdiff(unique(as.character(chrom)), contigs)
  missing_chrom <- missing_chrom[!is.na(missing_chrom)]
  if (length(missing_chrom) > 0) {
    warning(length(missing_chrom), " chromosome(s) in the VCF are not listed in the ##contig lines of vcf_header: ",
            paste(utils::head(missing_chrom, 5), collapse = ", "),
            if (length(missing_chrom) > 5) paste0(" and ", length(missing_chrom) - 5, " more") else "",
            ". Check that the header is for the same reference genome.", call. = FALSE)
  }
  invisible()
}
