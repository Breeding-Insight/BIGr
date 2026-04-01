#Internal Functions

utils::globalVariables(c(
  "ALT", "AlleleID", "AlleleSequence", "CHROM", "Concordance", "Data", "ID",
  "MarkerName", "POS",
  "QPseparate", "QPsolve_par", "REF", "Type", "Var1", "Variant", "geno",
  "ind", "ref", "row_name", "size", "snp",
  "CloneID", "Count", "qualifying_sites_count",
  "MarkerID", "SampleID", "Dosage",
  "pos", "alt", "match_key"
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

##' Verbose Message Utility
##'
##' Prints a formatted verbose message with timestamp, indentation, and type label, if verbose is TRUE.
##'
##' @param text Character string, the message to print (supports sprintf formatting).
##' @param verbose Logical. If TRUE, prints the message; if FALSE, suppresses output.
##' @param level Integer, indentation level (0=header, 1=main step, 2=detail, 3=sub-detail).
##' @param type Character string, message type (e.g., "INFO", "WARN", "ERROR"). Only shown for level 0.
##' @param ... Additional arguments passed to sprintf for formatting.
##'
##' @details Use the verbose argument to control message output. Typically, pass the function's verbose parameter to vmsg.
##'
##' @return No return value, called for side effects.
##' @export
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
