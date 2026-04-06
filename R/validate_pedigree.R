#' Validate Pedigree Trios Using Mendelian Error Analysis
#'
#' Validates parent-offspring trios by calculating Mendelian error rates from
#' SNP genotype data. Identifies incorrect parentage assignments and optionally
#' suggests or fills in best-matching replacements.
#'
#' @param pedigree_file  Character. Path to the pedigree file (TSV/CSV) with
#'                       columns: \code{Progeny}, \code{Sire}, \code{Dam}.
#' @param genotypes_file Character. Path to the genotypes file (TSV/CSV) with
#'                       an \code{ID} column followed by marker columns coded
#'                       as 0, 1, 2 (additive allele dosage).
#' @param trio_error_threshold     Numeric. Maximum Mendelian error percentage to
#'                            classify a trio as \code{PASS} (default: \code{5.0}).
#' @param min_markers         Integer. Minimum number of non-missing markers
#'                            required to evaluate a trio; below this the trio
#'                            is flagged \code{LOW_MARKERS} (default: \code{10}).
#' @param single_parent_error_threshold Numeric. Maximum homozygous-marker mismatch
#'                             percentage for a parent to be considered acceptable
#'                             in a failed trio (default: \code{2.0}).
#' @param fill_pedigree Logical. If \code{TRUE}, writes an additional file with
#'                      failed parents replaced by the best-matching candidate
#'                      IDs (default: \code{FALSE}).
#' @param verbose    Logical. If \code{TRUE}, prints progress messages, summary
#'                   statistics, and the results table to the console
#'                   (default: \code{TRUE}).
#' @param write_txt  Logical. If \code{TRUE}, writes the validation results to
#'                   \code{output_filename} (default: \code{TRUE}).
#' @param output_filename Character. Name of the output file for validation
#'                        results (default: \code{"trio_validation_results.txt"}).
#'
#' @return A \code{data.table} (returned invisibly) with one row per trio and
#'   the following columns:
#'   \describe{
#'     \item{Progeny}{Progeny ID.}
#'     \item{Sire}{Declared sire ID.}
#'     \item{Dam}{Declared dam ID.}
#'     \item{Mendelian_Error_Pct}{Trio-level Mendelian error percentage.}
#'     \item{Markers_Tested}{Number of markers compared across all three individuals.}
#'     \item{Status}{One of \code{PASS}, \code{FAIL}, \code{LOW_MARKERS}, or \code{NO_DATA}.}
#'     \item{Correction_Decision}{One of \code{NONE}, \code{KEEP_BOTH}, \code{REMOVE_SIRE},
#'           \code{REMOVE_DAM}, or \code{REMOVE_BOTH}.}
#'     \item{Sire_Hom_Error_Pct}{Sire homozygous-marker mismatch percentage (\code{NA} unless \code{FAIL}).}
#'     \item{Dam_Hom_Error_Pct}{Dam homozygous-marker mismatch percentage (\code{NA} unless \code{FAIL}).}
#'     \item{Best_Sire}{Best-matching sire candidate ID (\code{NA} unless sire removed).}
#'     \item{Best_Sire_Error_Pct}{Homozygous mismatch percentage for \code{Best_Sire}.}
#'     \item{Best_Dam}{Best-matching dam candidate ID (\code{NA} unless dam removed).}
#'     \item{Best_Dam_Error_Pct}{Homozygous mismatch percentage for \code{Best_Dam}.}
#'   }
#'
#' @details
#' Trios are filtered to individuals present in the genotype file before
#' analysis. Mendelian errors are counted as genotype combinations impossible
#' under Mendelian inheritance (e.g. both parents homozygous reference but
#' progeny carries the alternate allele). Failed trios are further dissected
#' using homozygous-only markers to identify which parent is likely incorrect.
#' A corrected pedigree with failed parents set to \code{0} is always written
#' to \code{corrected_pedigree.txt}. If \code{fill_pedigree = TRUE}, a second
#' file (\code{filled_pedigree.txt}) replaces those zeros with the best
#' genomic match.
#'
#' @examples
#' \dontrun{
#' # Basic run with defaults
#' results <- validate_pedigree("pedigree.txt", "genotypes.txt")
#'
#' # Stricter thresholds, custom output name, no console output
#' results <- validate_pedigree(
#'   pedigree_file      = "pedigree.txt",
#'   genotypes_file     = "genotypes.txt",
#'   trio_error_threshold    = 2.0,
#'   single_parent_error_threshold = 1.0,
#'   fill_pedigree      = TRUE,
#'   verbose            = FALSE,
#'   output_filename    = "my_validation.txt"
#' )
#' }
#'
#' @import data.table
#' @export
 validate_pedigree <- function(pedigree_file, genotypes_file,
trio_error_threshold = 5.0,
min_markers = 10,
single_parent_error_threshold = 2.0,
fill_pedigree = FALSE,
verbose = TRUE,
write_txt = TRUE,
output_filename = "trio_validation_results.txt") {  
  
  #### Input Validation and Data Loading ####
  if (trio_error_threshold < 0 || trio_error_threshold > 100) {
    stop("trio_error_threshold must be between 0 and 100")
  }
  if (single_parent_error_threshold < 0 || single_parent_error_threshold > 100) {
    stop("single_parent_error_threshold must be between 0 and 100")
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
      
      if (error_pct <= trio_error_threshold) {
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
        
        sire_acceptable <- !is.na(sire_error_pct) && sire_error_pct <= single_parent_error_threshold
        dam_acceptable  <- !is.na(dam_error_pct)  && dam_error_pct  <= single_parent_error_threshold
        
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
  corrected_pedigree <- copy(original_pedigree)
  
  for (i in seq_len(nrow(final_df))) {
    prog_id  <- final_df$Progeny[i]
    decision <- final_df$Correction_Decision[i]
    
    if (decision == "REMOVE_SIRE") {
      corrected_pedigree[Progeny == prog_id, Sire := 0]
    } else if (decision == "REMOVE_DAM") {
      corrected_pedigree[Progeny == prog_id, Dam := 0]
    } else if (decision == "REMOVE_BOTH") {
      corrected_pedigree[Progeny == prog_id, `:=`(Sire = 0, Dam = 0)]
    }
    # KEEP_BOTH and NONE require no changes
  }
  
  tryCatch({
    fwrite(corrected_pedigree, file = "corrected_pedigree.txt", sep = "\t", quote = FALSE)
    if (verbose) cat("Corrected pedigree (zeros) written to: corrected_pedigree.txt\n")
  }, error = function(e) {
    warning("Could not write corrected pedigree to file. Error: ", e$message, call. = FALSE)
  })
  
  #### Optionally Write Filled Pedigree (best matching IDs for failed parents) ####
  if (fill_pedigree) {
    filled_pedigree <- copy(original_pedigree)
    
    for (i in seq_len(nrow(final_df))) {
      prog_id  <- final_df$Progeny[i]
      decision <- final_df$Correction_Decision[i]
      
      if (decision == "REMOVE_SIRE") {
        filled_pedigree[Progeny == prog_id, Sire := final_df$Best_Sire[i]]
      } else if (decision == "REMOVE_DAM") {
        filled_pedigree[Progeny == prog_id, Dam := final_df$Best_Dam[i]]
      } else if (decision == "REMOVE_BOTH") {
        filled_pedigree[Progeny == prog_id, `:=`(Sire = final_df$Best_Sire[i],
                                                 Dam  = final_df$Best_Dam[i])]
      }
      # KEEP_BOTH and NONE require no changes
    }
    
    tryCatch({
      fwrite(filled_pedigree, file = "filled_pedigree.txt", sep = "\t", quote = FALSE)
      if (verbose) cat("Filled pedigree (best IDs) written to: filled_pedigree.txt\n")
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
    cat("Error threshold:", trio_error_threshold, "%\n")
    cat("Homozygous threshold:", single_parent_error_threshold, "%\n")
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
    tryCatch({
      fwrite(final_df, file = output_filename, sep = "\t", quote = FALSE)  # <-- uses new arg
      if (verbose) cat("Results written to:", output_filename, "\n")
    }, error = function(e) {
      warning("Could not write results to file. Error: ", e$message, call. = FALSE)
    })
  }
  
  if (verbose) print(final_df)
  
  return(invisible(final_df))  
 }
 