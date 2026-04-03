#' Validate Pedigree Using Mendelian Error Analysis
#'
#' Validates parent-offspring trios by calculating Mendelian error rates from
#' SNP genotype data. Trios exceeding the error threshold are flagged and
#' analysed further using homozygous-marker comparisons to identify which
#' parent(s) are likely incorrect. Optionally writes corrected pedigree files
#' with erroneous parents replaced by zeros or by the best-matching candidate IDs.
#'
#' @param pedigree_file Character. Path to the pedigree file (TSV or CSV).
#'   Must contain columns \code{Progeny}, \code{Sire}, and \code{Dam}.
#' @param genotypes_file Character. Path to the genotypes file (TSV or CSV).
#'   Must contain an \code{ID} column followed by one column per marker,
#'   with genotypes coded as 0, 1, or 2 (allele dosage).
#' @param error_threshold Numeric. Maximum acceptable Mendelian error percentage
#'   (0-100) for a trio to be considered a \code{PASS}. Default is \code{5.0}.
#' @param min_markers Integer. Minimum number of non-missing markers required
#'   across all three individuals for a trio to be fully evaluated. Trios below
#'   this threshold receive status \code{LOW_MARKERS}. Default is \code{10}.
#' @param homozygous_threshold Numeric. Maximum acceptable homozygous-marker
#'   mismatch percentage (0-100) for a parent to be considered acceptable in a
#'   failed trio. Default is \code{2.0}.
#' @param fill_pedigree Logical. If \code{TRUE}, writes an additional corrected
#'   pedigree file (\code{id_corrected_pedigree.txt}) in which removed parents
#'   are replaced by the best-matching candidate ID rather than zero.
#'   Default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, prints a summary of results and
#'   correction decisions to the console. Default is \code{TRUE}.
#' @param write_txt Logical. If \code{TRUE}, writes the full validation results
#'   table to \code{pedigree_validation_results.txt}. Default is \code{TRUE}.
#'
#' @return A \code{data.table} with one row per validated trio containing the
#'   following columns: \code{Progeny} (progeny ID); \code{Sire} (declared sire
#'   ID); \code{Dam} (declared dam ID); \code{Mendelian_Error_Pct} (overall
#'   Mendelian error rate \% across all three individuals);
#'   \code{Markers_Tested} (number of markers with non-missing genotypes in all
#'   three individuals); \code{Status} (validation outcome: \code{PASS},
#'   \code{FAIL}, \code{LOW_MARKERS}, or \code{NO_DATA});
#'   \code{Correction_Decision} (action taken for failed trios: \code{NONE},
#'   \code{KEEP_BOTH}, \code{REMOVE_SIRE}, \code{REMOVE_DAM}, or
#'   \code{REMOVE_BOTH}); \code{Sire_Hom_Error_Pct} (homozygous-marker mismatch
#'   \% between sire and progeny, \code{NA} unless status is \code{FAIL});
#'   \code{Dam_Hom_Error_Pct} (homozygous-marker mismatch \% between dam and
#'   progeny, \code{NA} unless status is \code{FAIL}); \code{Best_Sire} (ID of
#'   the best-matching sire candidate, populated only when
#'   \code{Correction_Decision} is \code{REMOVE_SIRE} or \code{REMOVE_BOTH});
#'   \code{Best_Sire_Error_Pct} (homozygous mismatch \% for \code{Best_Sire});
#'   \code{Best_Dam} (ID of the best-matching dam candidate, populated only when
#'   \code{Correction_Decision} is \code{REMOVE_DAM} or \code{REMOVE_BOTH});
#'   \code{Best_Dam_Error_Pct} (homozygous mismatch \% for \code{Best_Dam}).
#'   The function also writes \code{zero_corrected_pedigree.txt} (always) and,
#'   if \code{fill_pedigree = TRUE}, \code{id_corrected_pedigree.txt}. If
#'   \code{write_txt = TRUE}, results are written to
#'   \code{pedigree_validation_results.txt}. The return value is invisible when
#'   \code{verbose = TRUE}.
#'
#' @details
#' Mendelian errors are identified using standard allele-dosage rules, e.g. a
#' progeny cannot carry an allele absent in both parents. Only homozygous
#' parental markers (coded 0 or 2) are used in the per-parent mismatch
#' analysis, as heterozygous markers are uninformative for tracing allele
#' origin. Trios in the pedigree that lack genotype data for any of the three
#' individuals are removed prior to analysis.
#'
#' @importFrom data.table fread fwrite copy as.data.table rbindlist data.table
#'
#' @examples
#' results <- validate_pedigree(
#'   pedigree_file        = "pedigree.txt",
#'   genotypes_file       = "genotypes.txt",
#'   error_threshold      = 5.0,
#'   min_markers          = 10,
#'   homozygous_threshold = 2.0,
#'   fill_pedigree        = TRUE,
#'   verbose              = TRUE,
#'   write_txt            = TRUE
#' )
validate_pedigree <- function(pedigree_file, genotypes_file,
                                    error_threshold = 5.0,
                                    min_markers = 10,
                                    homozygous_threshold = 2.0,
                                    fill_pedigree = FALSE,
                                    verbose = TRUE,
                                    write_txt = TRUE) {

  #### Input Validation and Data Loading ####
  if (error_threshold < 0 || error_threshold > 100) {
    stop("error_threshold must be between 0 and 100")
  }
  if (homozygous_threshold < 0 || homozygous_threshold > 100) {
    stop("homozygous_threshold must be between 0 and 100")
  }

  tryCatch({
    pedigree <- fread(pedigree_file)
    genos    <- fread(genotypes_file)
  }, error = function(e) {
    stop("Error reading input files. Ensure paths are correct and files are TSV/CSV.")
  })

  # Keep original pedigree for correction if needed
  original_pedigree <- copy(pedigree)

  # Check required columns
  required_ped_cols <- c("Progeny", "Sire", "Dam")
  missing_cols <- setdiff(required_ped_cols, names(pedigree))
  if (length(missing_cols) > 0) {
    stop("Pedigree file missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!"ID" %in% names(genos)) {
    stop("Genotypes file must have an 'ID' column")
  }

  # Filter pedigree to only include individuals with genotype data
  valid_ids     <- genos$ID
  initial_trios <- nrow(pedigree)

  pedigree <- pedigree[Progeny %in% valid_ids & Sire %in% valid_ids & Dam %in% valid_ids]

  removed_trios <- initial_trios - nrow(pedigree)
  if (removed_trios > 0 && verbose) {
    cat("Removed", removed_trios, "trios due to missing genotype data.\n")
  }

  if (nrow(pedigree) == 0) {
    stop("No valid trios remain after filtering for genotype availability.")
  }

  #### Mendelian Error Calculation ####
  genos_mat <- as.matrix(genos, rownames = "ID")

  # Create homozygous-only matrix for parent analysis
  genos_hom <- copy(genos)
  marker_cols <- setdiff(names(genos_hom), "ID")
  for (col in marker_cols) {
    genos_hom[get(col) == 1, (col) := NA_integer_]
  }
  genos_hom_mat <- as.matrix(genos_hom, rownames = "ID")

  #### Helper: Find Best Matching Parent ####
  # Returns list(id, error_pct) for the candidate with lowest homozygous mismatch vs progeny
  find_best_parent <- function(prog_id, exclude_ids = character(0)) {
    candidates <- setdiff(rownames(genos_hom_mat), c(prog_id, exclude_ids))

    if (length(candidates) == 0) return(list(id = NA_character_, error_pct = NA_real_))

    prog_hom <- genos_hom_mat[prog_id, ]

    errors <- sapply(candidates, function(cand_id) {
      cand_hom    <- genos_hom_mat[cand_id, ]
      comparisons <- sum(!is.na(cand_hom) & !is.na(prog_hom))
      if (comparisons == 0) return(NA_real_)
      (sum(cand_hom != prog_hom, na.rm = TRUE) / comparisons) * 100
    })

    best_idx <- which.min(errors)
    list(id = candidates[best_idx], error_pct = round(errors[best_idx], 2))
  }

  results_list <- lapply(seq_len(nrow(pedigree)), function(i) {
    prog_id <- pedigree$Progeny[i]
    sire_id <- pedigree$Sire[i]
    dam_id  <- pedigree$Dam[i]

    # Extract genotype vectors
    progeny_vec <- genos_mat[prog_id, ]
    sire_vec    <- genos_mat[sire_id, ]
    dam_vec     <- genos_mat[dam_id,  ]

    # Calculate Mendelian errors using same logic as original function
    mismatches <- sum(
      (sire_vec == 0 & dam_vec == 0 & progeny_vec > 0) |
        (sire_vec == 2 & dam_vec == 2 & progeny_vec < 2) |
        ((sire_vec == 0 & dam_vec == 1) | (sire_vec == 1 & dam_vec == 0)) & (progeny_vec == 2) |
        ((sire_vec == 2 & dam_vec == 1) | (sire_vec == 1 & dam_vec == 2)) & (progeny_vec == 0) |
        ((sire_vec == 0 & dam_vec == 2) | (sire_vec == 2 & dam_vec == 0)) & (progeny_vec != 1),
      na.rm = TRUE
    )

    # Count comparable markers (non-NA in all three individuals)
    markers_tested <- sum(!is.na(sire_vec) & !is.na(dam_vec) & !is.na(progeny_vec))

    # Initialise per-parent and best-parent fields (populated only for FAILs)
    sire_error_pct <- NA_real_
    dam_error_pct  <- NA_real_
    best_sire      <- NA_character_
    best_sire_pct  <- NA_real_
    best_dam       <- NA_character_
    best_dam_pct   <- NA_real_

    # Calculate error percentage and determine status
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

      if (error_pct <= error_threshold) {
        status              <- "PASS"
        correction_decision <- "NONE"
      } else {
        status <- "FAIL"

        # Per-parent homozygous analysis for failed trios
        progeny_hom <- genos_hom_mat[prog_id, ]
        sire_hom    <- genos_hom_mat[sire_id, ]
        dam_hom     <- genos_hom_mat[dam_id,  ]

        # Sire homozygous error
        sire_comparisons <- sum(!is.na(sire_hom) & !is.na(progeny_hom))
        sire_error_pct   <- if (sire_comparisons == 0) NA_real_ else
          round((sum(sire_hom != progeny_hom, na.rm = TRUE) / sire_comparisons) * 100, 2)

        # Dam homozygous error
        dam_comparisons <- sum(!is.na(dam_hom) & !is.na(progeny_hom))
        dam_error_pct   <- if (dam_comparisons == 0) NA_real_ else
          round((sum(dam_hom != progeny_hom, na.rm = TRUE) / dam_comparisons) * 100, 2)

        sire_acceptable <- !is.na(sire_error_pct) && sire_error_pct <= homozygous_threshold
        dam_acceptable  <- !is.na(dam_error_pct)  && dam_error_pct  <= homozygous_threshold

        if (sire_acceptable && dam_acceptable) {
          correction_decision <- "KEEP_BOTH"
        } else if (sire_acceptable && !dam_acceptable) {
          correction_decision <- "REMOVE_DAM"
          best         <- find_best_parent(prog_id, exclude_ids = c(sire_id))
          best_dam     <- best$id
          best_dam_pct <- best$error_pct
        } else if (!sire_acceptable && dam_acceptable) {
          correction_decision <- "REMOVE_SIRE"
          best          <- find_best_parent(prog_id, exclude_ids = c(dam_id))
          best_sire     <- best$id
          best_sire_pct <- best$error_pct
        } else {
          correction_decision <- "REMOVE_BOTH"
          best_s        <- find_best_parent(prog_id, exclude_ids = character(0))
          best_sire     <- best_s$id
          best_sire_pct <- best_s$error_pct
          best_d        <- find_best_parent(prog_id, exclude_ids = c(best_s$id))
          best_dam      <- best_d$id
          best_dam_pct  <- best_d$error_pct
        }
      }
    }

    data.table(
      Progeny             = prog_id,
      Sire                = sire_id,
      Dam                 = dam_id,
      Mendelian_Error_Pct = round(error_pct, 2),
      Markers_Tested      = markers_tested,
      Status              = status,
      Correction_Decision = correction_decision,
      Sire_Hom_Error_Pct  = sire_error_pct,      # NA unless FAIL
      Dam_Hom_Error_Pct   = dam_error_pct,        # NA unless FAIL
      Best_Sire           = best_sire,
      Best_Sire_Error_Pct = best_sire_pct,
      Best_Dam            = best_dam,
      Best_Dam_Error_Pct  = best_dam_pct
    )
  })

  final_df <- rbindlist(results_list)

  #### Always Write Corrected Pedigree (zeros for failed parents) ####
  zero_corrected_pedigree <- copy(original_pedigree)

  for (i in seq_len(nrow(final_df))) {
    prog_id  <- final_df$Progeny[i]
    decision <- final_df$Correction_Decision[i]

    if (decision == "REMOVE_SIRE") {
      zero_corrected_pedigree[Progeny == prog_id, Sire := 0]
    } else if (decision == "REMOVE_DAM") {
      zero_corrected_pedigree[Progeny == prog_id, Dam := 0]
    } else if (decision == "REMOVE_BOTH") {
      zero_corrected_pedigree[Progeny == prog_id, `:=`(Sire = 0, Dam = 0)]
    }
    # KEEP_BOTH and NONE require no changes
  }

  tryCatch({
    fwrite(zero_corrected_pedigree, file = "zero_corrected_pedigree.txt", sep = "\t", quote = FALSE)
    if (verbose) cat("Corrected pedigree (zeros) written to: zero_corrected_pedigree.txt\n")
  }, error = function(e) {
    warning("Could not write corrected pedigree to file. Error: ", e$message, call. = FALSE)
  })

  #### Optionally Write Filled Pedigree (best matching IDs for failed parents) ####
  if (fill_pedigree) {
    id_corrected_pedigree <- copy(original_pedigree)

    for (i in seq_len(nrow(final_df))) {
      prog_id  <- final_df$Progeny[i]
      decision <- final_df$Correction_Decision[i]

      if (decision == "REMOVE_SIRE") {
        id_corrected_pedigree[Progeny == prog_id, Sire := final_df$Best_Sire[i]]
      } else if (decision == "REMOVE_DAM") {
        id_corrected_pedigree[Progeny == prog_id, Dam := final_df$Best_Dam[i]]
      } else if (decision == "REMOVE_BOTH") {
        id_corrected_pedigree[Progeny == prog_id, `:=`(Sire = final_df$Best_Sire[i],
                                                  Dam  = final_df$Best_Dam[i])]
      }
      # KEEP_BOTH and NONE require no changes
    }

    tryCatch({
      fwrite(id_corrected_pedigree, file = "id_corrected_pedigree.txt", sep = "\t", quote = FALSE)
      if (verbose) cat("Filled pedigree (best IDs) written to: id_corrected_pedigree.txt\n")
    }, error = function(e) {
      warning("Could not write filled pedigree to file. Error: ", e$message, call. = FALSE)
    })
  }

  #### Summary Statistics ####
  if (verbose) {
    total_trios   <- nrow(final_df)
    status_counts <- table(final_df$Status)

    cat("\n--- Trio Validation Summary ---\n")
    cat("Total trios tested:", total_trios, "\n")
    for (status in names(status_counts)) {
      cat(sprintf("%-12s: %d (%.1f%%)\n", status, status_counts[status],
                  (status_counts[status] / total_trios) * 100))
    }
    cat("Error threshold:", error_threshold, "%\n")
    cat("Homozygous threshold:", homozygous_threshold, "%\n")
    cat("Minimum markers required:", min_markers, "\n\n")

    corrections <- table(final_df$Correction_Decision)
    cat("Correction summary:\n")
    for (decision in names(corrections)) {
      if (decision != "NONE") {
        cat("  ", decision, ":", corrections[decision], "\n")
      }
    }
    cat("\n")
  }

  #### Output ####
  if (write_txt) {
    output_filename <- "pedigree_validation_results.txt"
    tryCatch({
      fwrite(final_df, file = output_filename, sep = "\t", quote = FALSE)
      if (verbose) cat("Results written to:", output_filename, "\n")
    }, error = function(e) {
      warning("Could not write results to file. Error: ", e$message, call. = FALSE)
    })
  }

  if (verbose) {
    print(final_df)
    return(invisible(final_df))
  } else {
    return(final_df)
  }
}
