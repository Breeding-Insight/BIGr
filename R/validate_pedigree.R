#' Validate Pedigree Trios Using Mendelian Error Analysis
#'
#' Validates parent-offspring trios by calculating Mendelian error rates from
#' SNP genotype data. Identifies incorrect parentage assignments and optionally
#' suggests or fills in best-matching replacements.
#'
#' @param pedigree_file Character. Path to the pedigree file (TSV/CSV/TXT) with
#'   columns: \code{ID}, \code{Male_Parent}, \code{Female_Parent}.
#' @param genotypes_file Character. Path to the genotypes file (TSV/CSV/TXT)
#'   with an \code{ID} column followed by marker columns coded as 0, 1, 2
#'   (additive allele dosage).
#' @param trio_error_threshold Numeric. Maximum Mendelian error percentage to
#'   classify a trio as \code{PASS} (default: \code{5.0}). Must be between
#'   0 and 100.
#' @param min_markers Integer. Minimum number of non-missing markers required
#'   to evaluate a trio; below this the trio is flagged \code{LOW_MARKERS}
#'   (default: \code{10}).
#' @param single_parent_error_threshold Numeric. Maximum homozygous-marker
#'   mismatch percentage for a parent to be considered acceptable in a failed
#'   trio (default: \code{2.0}). Must be between 0 and 100.
#' @param fill_pedigree Logical. If \code{TRUE}, writes an additional file
#'   (\code{filled_pedigree.txt}) with failed parents replaced by the
#'   best-matching candidate IDs (default: \code{FALSE}).
#' @param verbose Logical. If \code{TRUE}, prints progress messages, a summary
#'   table, and the results to the console (default: \code{TRUE}).
#' @param write_txt Logical. If \code{TRUE}, writes the validation results
#'   to \code{output_filename} (default: \code{TRUE}).
#' @param output_filename Character. Path/name of the output file for
#'   validation results (default: \code{"pedigree_validation_results.txt"}).
#'
#' @return A \code{data.table} (returned invisibly) with one row per trio and
#'   the following columns:
#'   \describe{
#'     \item{ID}{Individual ID (first column of the pedigree input).}
#'     \item{Male_Parent}{Declared male parent ID.}
#'     \item{Female_Parent}{Declared female parent ID.}
#'     \item{Mendelian_Error_Pct}{Trio-level Mendelian error percentage.}
#'     \item{Markers_Tested}{Number of markers with non-missing genotypes in
#'           all three individuals.}
#'     \item{Status}{One of \code{PASS}, \code{FAIL}, \code{LOW_MARKERS}, or
#'           \code{NO_DATA}.}
#'     \item{Correction_Decision}{One of \code{NONE}, \code{KEEP_BOTH},
#'           \code{REMOVE_MALE_PARENT}, \code{REMOVE_FEMALE_PARENT}, or
#'           \code{REMOVE_BOTH}.}
#'     \item{Male_Parent_Hom_Error_Pct}{Male parent homozygous-marker mismatch
#'           percentage (\code{NA} unless \code{Status == "FAIL"}).}
#'     \item{Female_Parent_Hom_Error_Pct}{Female parent homozygous-marker
#'           mismatch percentage (\code{NA} unless \code{Status == "FAIL"}).}
#'     \item{Best_Male_Parent}{Best-matching male parent candidate ID
#'           (\code{NA} unless male parent is removed).}
#'     \item{Best_Male_Parent_Error_Pct}{Homozygous mismatch percentage for
#'           \code{Best_Male_Parent}.}
#'     \item{Best_Female_Parent}{Best-matching female parent candidate ID
#'           (\code{NA} unless female parent is removed).}
#'     \item{Best_Female_Parent_Error_Pct}{Homozygous mismatch percentage for
#'           \code{Best_Female_Parent}.}
#'   }
#'
#' @details
#' Trios are filtered to individuals present in the genotype file before
#' analysis. Mendelian errors are counted as genotype combinations impossible
#' under Mendelian inheritance (e.g. both parents homozygous reference but the
#' offspring carries the alternate allele). Failed trios are further dissected
#' using homozygous-only markers to identify which parent is likely incorrect.
#'
#' A corrected pedigree with failed parents replaced by \code{0} is always
#' written to \code{corrected_pedigree.txt} in the working directory. If
#' \code{fill_pedigree = TRUE}, a second file (\code{filled_pedigree.txt})
#' replaces those zeros with the best genomic match found across all
#' genotyped individuals.
#'
#' @examples
#' \dontrun{
#' # Basic run with defaults
#' results <- validate_pedigree("pedigree.txt", "genotypes.txt")
#'
#' # Stricter thresholds, fill replacements, suppress console output
#' results <- validate_pedigree(
#'   pedigree_file                 = "pedigree.txt",
#'   genotypes_file                = "genotypes.txt",
#'   trio_error_threshold          = 2.0,
#'   single_parent_error_threshold = 1.0,
#'   fill_pedigree                 = TRUE,
#'   verbose                       = FALSE,
#'   output_filename               = "my_validation.txt"
#' )
#' }
#'
#' @importFrom data.table fread fwrite rbindlist copy data.table := set
#' @export
validate_pedigree <- function(pedigree_file, genotypes_file,
                              trio_error_threshold = 5.0,
                              min_markers = 10,
                              single_parent_error_threshold = 2.0,
                              fill_pedigree = FALSE,
                              verbose = TRUE,
                              write_txt = TRUE,
                              output_filename = "pedigree_validation_results.txt") {

  #### Input Validation ####
  if (trio_error_threshold < 0 || trio_error_threshold > 100) {
    stop("trio_error_threshold must be between 0 and 100")
  }
  if (single_parent_error_threshold < 0 || single_parent_error_threshold > 100) {
    stop("single_parent_error_threshold must be between 0 and 100")
  }

  tryCatch({
    pedigree <- data.table::fread(pedigree_file,  sep = "auto")
    genos    <- data.table::fread(genotypes_file, sep = "auto")
  }, error = function(e) {
    stop("Error reading input files. Ensure paths are correct and files are TXT/TSV/CSV.")
  })

  #### Check Required Columns ####
  required_ped_cols <- c("ID", "Male_Parent", "Female_Parent")   # <-- changed
  missing_cols <- base::setdiff(required_ped_cols, base::names(pedigree))
  if (base::length(missing_cols) > 0) {
    stop("Pedigree file missing required columns: ",
         base::paste(missing_cols, collapse = ", "))
  }
  if (!"ID" %in% base::names(genos)) {
    stop("Genotypes file must have an 'ID' column")
  }

  # Keep original pedigree for correction output
  original_pedigree <- data.table::copy(pedigree)

  #### Filter to Individuals with Genotype Data ####
  valid_ids     <- genos$ID
  initial_trios <- base::nrow(pedigree)
  pedigree <- pedigree[ID %in% valid_ids &                        # <-- changed
                         Male_Parent %in% valid_ids &
                         Female_Parent %in% valid_ids]
  removed_trios <- initial_trios - base::nrow(pedigree)
  if (removed_trios > 0 && verbose) {
    base::cat("Removed", removed_trios, "trios due to missing genotype data.\n")
  }
  if (base::nrow(pedigree) == 0) {
    stop("No valid trios remain after filtering for genotype availability.")
  }

  #### Mendelian Error Calculation ####
  genos_mat <- base::as.matrix(genos, rownames = "ID")

  # Homozygous-only matrix
  genos_hom   <- data.table::copy(genos)
  marker_cols <- base::setdiff(base::names(genos_hom), "ID")
  for (col in marker_cols) {
    genos_hom[base::get(col) == 1, (col) := NA_integer_]
  }
  genos_hom_mat <- base::as.matrix(genos_hom, rownames = "ID")

  #### Helper: Find Best Matching Parent ####
  find_best_parent <- function(prog_id, exclude_ids = base::character(0)) {
    candidates <- base::setdiff(base::rownames(genos_hom_mat),
                                c(prog_id, exclude_ids))
    if (base::length(candidates) == 0) {
      return(base::list(id = NA_character_, error_pct = NA_real_))
    }
    prog_hom <- genos_hom_mat[prog_id, ]
    errors <- base::sapply(candidates, function(cand_id) {
      cand_hom    <- genos_hom_mat[cand_id, ]
      comparisons <- base::sum(!base::is.na(cand_hom) & !base::is.na(prog_hom))
      if (comparisons == 0) return(NA_real_)
      (base::sum(cand_hom != prog_hom, na.rm = TRUE) / comparisons) * 100
    })
    best_idx <- base::which.min(errors)
    base::list(id        = candidates[best_idx],
               error_pct = base::round(errors[best_idx], 2))
  }

  results_list <- base::lapply(base::seq_len(base::nrow(pedigree)), function(i) {
    prog_id          <- pedigree$ID[i]                            # <-- changed
    male_parent_id   <- pedigree$Male_Parent[i]
    female_parent_id <- pedigree$Female_Parent[i]

    progeny_vec       <- genos_mat[prog_id,          ]
    male_parent_vec   <- genos_mat[male_parent_id,   ]
    female_parent_vec <- genos_mat[female_parent_id, ]

    mismatches <- base::sum(
      (male_parent_vec == 0 & female_parent_vec == 0 & progeny_vec > 0) |
        (male_parent_vec == 2 & female_parent_vec == 2 & progeny_vec < 2) |
        ((male_parent_vec == 0 & female_parent_vec == 1) |
           (male_parent_vec == 1 & female_parent_vec == 0)) & (progeny_vec == 2) |
        ((male_parent_vec == 2 & female_parent_vec == 1) |
           (male_parent_vec == 1 & female_parent_vec == 2)) & (progeny_vec == 0) |
        ((male_parent_vec == 0 & female_parent_vec == 2) |
           (male_parent_vec == 2 & female_parent_vec == 0)) & (progeny_vec != 1),
      na.rm = TRUE
    )

    markers_tested <- base::sum(!base::is.na(male_parent_vec) &
                                  !base::is.na(female_parent_vec) &
                                  !base::is.na(progeny_vec))

    male_parent_error_pct   <- NA_real_
    female_parent_error_pct <- NA_real_
    best_male_parent        <- NA_character_
    best_male_parent_pct    <- NA_real_
    best_female_parent      <- NA_character_
    best_female_parent_pct  <- NA_real_

    if (markers_tested == 0) {
      error_pct           <- NA_real_
      status              <- "NO_DATA"
      correction_decision <- "NONE"
    } else if (markers_tested < min_markers) {
      error_pct           <- (mismatches / markers_tested) * 100
      status              <- "LOW_MARKERS"
      correction_decision <- "NONE"
    } else {
      error_pct <- (mismatches / markers_tested) * 100
      if (error_pct <= trio_error_threshold) {
        status              <- "PASS"
        correction_decision <- "NONE"
      } else {
        status <- "FAIL"

        progeny_hom       <- genos_hom_mat[prog_id,          ]
        male_parent_hom   <- genos_hom_mat[male_parent_id,   ]
        female_parent_hom <- genos_hom_mat[female_parent_id, ]

        male_comparisons      <- base::sum(!base::is.na(male_parent_hom) &
                                             !base::is.na(progeny_hom))
        male_parent_error_pct <- if (male_comparisons == 0) NA_real_ else
          base::round((base::sum(male_parent_hom != progeny_hom,
                                 na.rm = TRUE) / male_comparisons) * 100, 2)

        female_comparisons      <- base::sum(!base::is.na(female_parent_hom) &
                                               !base::is.na(progeny_hom))
        female_parent_error_pct <- if (female_comparisons == 0) NA_real_ else
          base::round((base::sum(female_parent_hom != progeny_hom,
                                 na.rm = TRUE) / female_comparisons) * 100, 2)

        male_acceptable   <- !base::is.na(male_parent_error_pct) &&
          male_parent_error_pct   <= single_parent_error_threshold
        female_acceptable <- !base::is.na(female_parent_error_pct) &&
          female_parent_error_pct <= single_parent_error_threshold

        if (male_acceptable && female_acceptable) {
          correction_decision <- "KEEP_BOTH"
        } else if (male_acceptable && !female_acceptable) {
          correction_decision    <- "REMOVE_FEMALE_PARENT"
          best                   <- find_best_parent(prog_id,
                                                     exclude_ids = c(male_parent_id))
          best_female_parent     <- best$id
          best_female_parent_pct <- best$error_pct
        } else if (!male_acceptable && female_acceptable) {
          correction_decision  <- "REMOVE_MALE_PARENT"
          best                 <- find_best_parent(prog_id,
                                                   exclude_ids = c(female_parent_id))
          best_male_parent     <- best$id
          best_male_parent_pct <- best$error_pct
        } else {
          correction_decision    <- "REMOVE_BOTH"
          best_m                 <- find_best_parent(prog_id,
                                                     exclude_ids = base::character(0))
          best_male_parent       <- best_m$id
          best_male_parent_pct   <- best_m$error_pct
          best_f                 <- find_best_parent(prog_id,
                                                     exclude_ids = c(best_m$id))
          best_female_parent     <- best_f$id
          best_female_parent_pct <- best_f$error_pct
        }
      }
    }

    data.table::data.table(
      ID                           = prog_id,           # <-- changed
      Male_Parent                  = male_parent_id,
      Female_Parent                = female_parent_id,
      Mendelian_Error_Pct          = base::round(error_pct, 2),
      Markers_Tested               = markers_tested,
      Status                       = status,
      Correction_Decision          = correction_decision,
      Male_Parent_Hom_Error_Pct    = male_parent_error_pct,
      Female_Parent_Hom_Error_Pct  = female_parent_error_pct,
      Best_Male_Parent             = best_male_parent,
      Best_Male_Parent_Error_Pct   = best_male_parent_pct,
      Best_Female_Parent           = best_female_parent,
      Best_Female_Parent_Error_Pct = best_female_parent_pct
    )
  })

  final_df <- data.table::rbindlist(results_list)

  #### Always Write Corrected Pedigree ####
  corrected_pedigree <- data.table::copy(original_pedigree)
  for (i in base::seq_len(base::nrow(final_df))) {
    prog_id  <- final_df$ID[i]                                   # <-- changed
    decision <- final_df$Correction_Decision[i]
    row_idx  <- base::which(corrected_pedigree$ID == prog_id)    # <-- changed
    if (decision == "REMOVE_MALE_PARENT") {
      data.table::set(corrected_pedigree, row_idx, "Male_Parent", 0L)
    } else if (decision == "REMOVE_FEMALE_PARENT") {
      data.table::set(corrected_pedigree, row_idx, "Female_Parent", 0L)
    } else if (decision == "REMOVE_BOTH") {
      data.table::set(corrected_pedigree, row_idx, "Male_Parent",   0L)
      data.table::set(corrected_pedigree, row_idx, "Female_Parent", 0L)
    }
  }
  tryCatch({
    data.table::fwrite(corrected_pedigree, file = "corrected_pedigree.txt",
                       sep = "\t", quote = FALSE)
    if (verbose) base::cat("Corrected pedigree written to: corrected_pedigree.txt\n")
  }, error = function(e) {
    warning("Could not write corrected pedigree. Error: ", e$message, call. = FALSE)
  })

  #### Optionally Write Filled Pedigree ####
  if (fill_pedigree) {
    filled_pedigree <- data.table::copy(original_pedigree)
    for (i in base::seq_len(base::nrow(final_df))) {
      prog_id  <- final_df$ID[i]                                 # <-- changed
      decision <- final_df$Correction_Decision[i]
      row_idx  <- base::which(filled_pedigree$ID == prog_id)     # <-- changed
      if (decision == "REMOVE_MALE_PARENT") {
        data.table::set(filled_pedigree, row_idx,
                        "Male_Parent", final_df$Best_Male_Parent[i])
      } else if (decision == "REMOVE_FEMALE_PARENT") {
        data.table::set(filled_pedigree, row_idx,
                        "Female_Parent", final_df$Best_Female_Parent[i])
      } else if (decision == "REMOVE_BOTH") {
        data.table::set(filled_pedigree, row_idx,
                        "Male_Parent",   final_df$Best_Male_Parent[i])
        data.table::set(filled_pedigree, row_idx,
                        "Female_Parent", final_df$Best_Female_Parent[i])
      }
    }
    tryCatch({
      data.table::fwrite(filled_pedigree, file = "filled_pedigree.txt",
                         sep = "\t", quote = FALSE)
      if (verbose) base::cat("Filled pedigree written to: filled_pedigree.txt\n")
    }, error = function(e) {
      warning("Could not write filled pedigree. Error: ", e$message, call. = FALSE)
    })
  }

  #### Summary ####
  if (verbose) {
    total_trios   <- base::nrow(final_df)
    status_counts <- base::table(final_df$Status)
    base::cat("\n--- Trio Validation Summary ---\n")
    base::cat("Total trios tested:", total_trios, "\n")
    for (s in base::names(status_counts)) {
      base::cat(base::sprintf("%-12s: %d (%.1f%%)\n", s,
                              status_counts[s],
                              (status_counts[s] / total_trios) * 100))
    }
    base::cat("Error threshold:", trio_error_threshold, "%\n")
    base::cat("Homozygous threshold:", single_parent_error_threshold, "%\n")
    base::cat("Minimum markers required:", min_markers, "\n\n")
    corrections <- base::table(final_df$Correction_Decision)
    base::cat("Correction summary:\n")
    for (decision in base::names(corrections)) {
      if (decision != "NONE")
        base::cat("  ", decision, ":", corrections[decision], "\n")
    }
    base::cat("\n")
    base::print(final_df)
  }

  #### Output ####
  if (write_txt) {
    tryCatch({
      data.table::fwrite(final_df, file = output_filename,
                         sep = "\t", quote = FALSE)
      if (verbose) base::cat("Results written to:", output_filename, "\n")
    }, error = function(e) {
      warning("Could not write results. Error: ", e$message, call. = FALSE)
    })
  }

  return(base::invisible(final_df))
}
