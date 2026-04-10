#' Find Parentage Assignments for Progeny
#'
#' Assigns the most likely parent(s) to each progeny individual based on
#' genotypic data using Mendelian error rates or homozygous mismatch rates.
#'
#' @param genotypes_file Path to a TSV/CSV/TXT file containing genotype data.
#'   Must include an 'ID' column followed by marker columns coded as 0, 1, 2
#'   (allele dosage).
#' @param parents_file Path to a TSV/CSV/TXT file listing candidate parent IDs.
#'   Must include an 'ID' column. An optional 'Sex' column with values
#'   'M' (male parent), 'F' (female parent), or 'A' (ambiguous) determines
#'   which parents are tested for each role. If absent, all parents are treated
#'   as ambiguous.
#' @param progeny_file Path to a TSV/CSV/TXT file listing progeny IDs to assign.
#'   Must include an 'ID' column.
#' @param method Character. Parentage assignment method. One of:
#'   \itemize{
#'     \item \code{"best_male_parent"}  — finds the best male parent for each
#'       progeny using homozygous mismatch rate.
#'     \item \code{"best_female_parent"} — finds the best female parent for each
#'       progeny using homozygous mismatch rate.
#'     \item \code{"best_match"} — finds the single best parent (either sex)
#'       using homozygous mismatch rate.
#'     \item \code{"best_pair"}  — finds the best male-female parent pair for
#'       each progeny using full Mendelian error rate (default).
#'   }
#' @param min_markers Integer. Minimum number of non-missing markers required
#'   to report a parentage assignment. Progeny-parent comparisons with fewer
#'   markers are flagged as \code{LOW_MARKERS} and no assignment is made
#'   (default: \code{10}).
#' @param error_threshold Numeric. Maximum mismatch percentage to report a
#'   parentage assignment as confident. Assignments above this threshold are
#'   flagged as \code{HIGH_ERROR} in the \code{Assignment_Status} column
#'   (default: \code{5.0}). Must be between 0 and 100.
#' @param show_ties Logical. If \code{TRUE}, all tied best pairs are reported
#'   as additional columns (\code{Male_Parent_1}, \code{Male_Parent_2}, etc.)
#'   when \code{method = "best_pair"}. If \code{FALSE}, only one tied pair is
#'   reported with a warning. Default is \code{TRUE}.
#' @param allow_selfing Logical. If \code{FALSE}, male-female parent pairs where
#'   both IDs are identical are excluded when \code{method = "best_pair"}.
#'   Default is \code{TRUE}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages, summary
#'   statistics, and the results table to the console. Default is \code{TRUE}.
#' @param write_txt Logical. If \code{TRUE}, writes results to
#'   \code{parentage_results_dt.txt} in the working directory. Default is
#'   \code{TRUE}.
#'
#' @return A \code{data.table} with one row per progeny (or more if ties are
#'   reported). Columns depend on the method used:
#'   \itemize{
#'     \item \code{best_male_parent} / \code{best_female_parent} / \code{best_match}:
#'       \code{Progeny}, \code{Best_Match}, \code{Mendelian_Error_Pct},
#'       \code{Markers_Tested}, \code{Assignment_Status}.
#'     \item \code{best_pair} (no ties): \code{Progeny}, \code{Male_Parent},
#'       \code{Female_Parent}, \code{Mendelian_Error_Pct}, \code{Markers_Tested},
#'       \code{Assignment_Status}.
#'     \item \code{best_pair} (with ties): columns are suffixed \code{_1},
#'       \code{_2}, etc. for each tied pair, plus \code{Assignment_Status}.
#'   }
#'   \code{Assignment_Status} is one of \code{PASS}, \code{HIGH_ERROR}, or
#'   \code{LOW_MARKERS}. Returned invisibly when \code{verbose = TRUE}.
#'
#' @details
#' For \code{"best_male_parent"}, \code{"best_female_parent"}, and
#' \code{"best_match"}, only homozygous markers (coded 0 or 2) are used for
#' comparison; heterozygous markers (coded 1) are set to \code{NA}. This
#' reduces false mismatches caused by phase ambiguity.
#'
#' For \code{"best_pair"}, all markers are used and full Mendelian inheritance
#' rules are applied across all possible male-female parent combinations via
#' \code{data.table::CJ()}.
#'
#' Individuals in \code{parents_file} or \code{progeny_file} that are absent
#' from \code{genotypes_file} are removed with a warning.
#'
#' Progeny with fewer non-missing markers than \code{min_markers} are flagged
#' \code{LOW_MARKERS} and no parent assignment is reported. Progeny whose best
#' match exceeds \code{error_threshold} are flagged \code{HIGH_ERROR}.
#'
#' @examples
#' \dontrun{
#' # Assign best male-female parent pair to each progeny
#' results <- find_parentage(
#'   genotypes_file = "genotypes.txt",
#'   parents_file   = "parents.txt",
#'   progeny_file   = "progeny.txt",
#'   method         = "best_pair",
#'   min_markers    = 50,
#'   error_threshold = 5.0,
#'   show_ties      = TRUE,
#'   allow_selfing  = FALSE
#' )
#'
#' # Find best individual parent match (ignoring sex)
#' results <- find_parentage(
#'   genotypes_file  = "genotypes.txt",
#'   parents_file    = "parents.txt",
#'   progeny_file    = "progeny.txt",
#'   method          = "best_match",
#'   min_markers     = 30,
#'   error_threshold = 3.0
#' )
#' }
#'
#' @importFrom data.table fread fwrite copy CJ rbindlist set data.table as.data.table
#' @export
find_parentage <- function(genotypes_file, parents_file, progeny_file,
                           method = "best_pair",
                           min_markers = 10,
                           error_threshold = 5.0,
                           show_ties = TRUE,
                           allow_selfing = TRUE,
                           verbose = TRUE,
                           write_txt = TRUE) {

  #### Input Validation and Data Loading ####
  allowed_methods <- c("best_male_parent", "best_female_parent", "best_match", "best_pair")
  if (!method %in% allowed_methods)
    stop("Method must be one of: ", paste(allowed_methods, collapse = ", "))
  if (min_markers < 1)
    stop("min_markers must be a positive integer.")
  if (error_threshold < 0 || error_threshold > 100)
    stop("error_threshold must be between 0 and 100.")

  tryCatch({
    genos              <- data.table::fread(genotypes_file, sep = "auto")
    all_parents        <- data.table::fread(parents_file,   sep = "auto")
    progeny_candidates <- data.table::fread(progeny_file,   sep = "auto")
  }, error = function(e) {
    stop("Error reading input files. Ensure paths are correct and files are TXT/TSV/CSV.")
  })

  valid_ids       <- genos$ID
  removed_parents <- base::setdiff(all_parents$ID, valid_ids)
  if (length(removed_parents) > 0) {
    warning("The following parent IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_parents, collapse = ", "), call. = FALSE)
    all_parents <- all_parents[ID %in% valid_ids]
  }

  removed_progeny <- base::setdiff(progeny_candidates$ID, valid_ids)
  if (length(removed_progeny) > 0) {
    warning("The following progeny IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_progeny, collapse = ", "), call. = FALSE)
    progeny_candidates <- progeny_candidates[ID %in% valid_ids]
  }

  if (!"Sex" %in% base::colnames(all_parents)) {
    warning("No 'Sex' column in parents file. All parents treated as ambiguous ('A').")
    all_parents[, Sex := "A"]
  }

  all_parents[, Sex := base::toupper(Sex)]
  male_parent_candidates   <- all_parents[Sex %in% c("M", "A", "NA"), .SD]
  female_parent_candidates <- all_parents[Sex %in% c("F", "A", "NA")]

  if (base::nrow(male_parent_candidates) == 0 && method %in% c("best_male_parent", "best_pair"))
    warning("No valid male parent candidates remain after filtering.", call. = FALSE)
  if (base::nrow(female_parent_candidates) == 0 && method %in% c("best_female_parent", "best_pair"))
    warning("No valid female parent candidates remain after filtering.", call. = FALSE)
  if (base::nrow(progeny_candidates) == 0)
    stop("No valid progeny candidates remain after filtering.")

  #### Helper: assign Assignment_Status from markers and error rate ####
  ## Returns LOW_MARKERS, HIGH_ERROR, or PASS
  assign_status <- function(markers, error_pct) {
    base::ifelse(markers < min_markers, "LOW_MARKERS",
                 base::ifelse(error_pct > error_threshold, "HIGH_ERROR", "PASS"))
  }

  #### Logic for Homozygous Matching Methods ####
  if (method %in% c("best_male_parent", "best_female_parent", "best_match")) {

    genos_hom   <- data.table::copy(genos)
    marker_cols <- base::setdiff(base::names(genos_hom), "ID")
    for (col in marker_cols)
      genos_hom[base::get(col) == 1, (col) := NA_integer_]

    parent_ids <- base::switch(method,
                               "best_male_parent"   = male_parent_candidates$ID,
                               "best_female_parent" = female_parent_candidates$ID,
                               "best_match"         = base::union(male_parent_candidates$ID,
                                                                  female_parent_candidates$ID))

    parent_genos  <- base::as.matrix(genos_hom[ID %in% parent_ids],              rownames = "ID")
    progeny_genos <- base::as.matrix(genos_hom[ID %in% progeny_candidates$ID],   rownames = "ID")

    results_list <- base::lapply(base::rownames(progeny_genos), function(progeny_id) {
      progeny_vec      <- progeny_genos[progeny_id, ]
      mismatches       <- base::rowSums(parent_genos != progeny_vec, na.rm = TRUE)
      comparisons      <- base::rowSums(!base::is.na(parent_genos) & !base::is.na(progeny_vec))
      percent_mismatch <- (mismatches / comparisons) * 100
      percent_mismatch[base::is.nan(percent_mismatch)] <- NA

      best_idx <- base::which.min(percent_mismatch)

      # No candidate found — return NA row flagged LOW_MARKERS
      if (base::length(best_idx) == 0) {
        return(data.table::data.table(
          Progeny             = progeny_id,
          Best_Match          = NA_character_,
          Mendelian_Error_Pct = NA_real_,
          Markers_Tested      = 0L,
          Assignment_Status   = "LOW_MARKERS"
        ))
      }

      best_markers <- comparisons[best_idx]
      best_error   <- base::round(percent_mismatch[best_idx], 2)

      data.table::data.table(
        Progeny             = progeny_id,
        Best_Match          = base::rownames(parent_genos)[best_idx],
        Mendelian_Error_Pct = best_error,
        Markers_Tested      = best_markers,
        Assignment_Status   = assign_status(best_markers, best_error)
      )
    })

    final_df <- data.table::rbindlist(results_list)
  }

  #### Logic for Best Pair Method ####
  if (method == "best_pair") {

    genos_mat    <- base::as.matrix(genos, rownames = "ID")
    parent_pairs <- data.table::CJ(Male_Parent   = male_parent_candidates$ID,
                                   Female_Parent = female_parent_candidates$ID)

    if (!allow_selfing) {
      parent_pairs <- parent_pairs[Male_Parent != Female_Parent]
      if (verbose) base::cat("Selfing is disallowed. Pairs with identical parents are removed.\n")
    }

    if (base::nrow(parent_pairs) == 0) stop("No valid parent pairs to test.")

    male_parent_genos_mat   <- genos_mat[parent_pairs$Male_Parent,   , drop = FALSE]
    female_parent_genos_mat <- genos_mat[parent_pairs$Female_Parent, , drop = FALSE]

    results_list <- base::lapply(progeny_candidates$ID, function(prog_id) {

      progeny_vec <- genos_mat[prog_id, ]

      mismatches <- base::rowSums(
        (male_parent_genos_mat == 0 & female_parent_genos_mat == 0 & progeny_vec > 0) |
          (male_parent_genos_mat == 2 & female_parent_genos_mat == 2 & progeny_vec < 2) |
          ((male_parent_genos_mat == 0 & female_parent_genos_mat == 1) |
             (male_parent_genos_mat == 1 & female_parent_genos_mat == 0)) & (progeny_vec == 2) |
          ((male_parent_genos_mat == 2 & female_parent_genos_mat == 1) |
             (male_parent_genos_mat == 1 & female_parent_genos_mat == 2)) & (progeny_vec == 0) |
          ((male_parent_genos_mat == 0 & female_parent_genos_mat == 2) |
             (male_parent_genos_mat == 2 & female_parent_genos_mat == 0)) & (progeny_vec != 1),
        na.rm = TRUE
      )

      comparisons      <- base::rowSums(!base::is.na(male_parent_genos_mat) &
                                          !base::is.na(female_parent_genos_mat) &
                                          !base::is.na(progeny_vec))
      percent_mismatch <- (mismatches / comparisons) * 100
      percent_mismatch[base::is.nan(percent_mismatch)] <- NA

      min_mismatch_val <- base::min(percent_mismatch, na.rm = TRUE)

      # No markers overlap at all — flag LOW_MARKERS, return early
      if (base::is.infinite(min_mismatch_val)) {
        return(data.table::data.table(
          Progeny           = prog_id,
          Markers_Tested    = 0L,
          Assignment_Status = "LOW_MARKERS"
        ))
      }

      best_indices <- base::which(percent_mismatch == min_mismatch_val)
      best_pairs   <- parent_pairs[best_indices]
      best_markers <- comparisons[best_indices[1]]
      best_error   <- base::round(min_mismatch_val, 2)
      a_status     <- assign_status(best_markers, best_error)

      if (!show_ties && base::nrow(best_pairs) > 1) {
        warning("Progeny '", prog_id, "' has ", base::nrow(best_pairs),
                " tied best pairs. Only one is reported as show_ties=FALSE.", call. = FALSE)
      }

      num_to_report <- base::min(base::nrow(best_pairs),
                                 if (show_ties) base::nrow(best_pairs) else 1)
      result_row    <- base::list(Progeny = prog_id)

      if (num_to_report == 1) {
        result_row[["Male_Parent"]]         <- best_pairs$Male_Parent[1]
        result_row[["Female_Parent"]]       <- best_pairs$Female_Parent[1]
        result_row[["Mendelian_Error_Pct"]] <- base::sprintf("%.2f", min_mismatch_val)
        result_row[["Markers_Tested"]]      <- best_markers
        result_row[["Assignment_Status"]]   <- a_status

      } else {
        for (k in base::seq_len(num_to_report)) {
          result_row[[base::paste0("Male_Parent_",         k)]] <- best_pairs$Male_Parent[k]
          result_row[[base::paste0("Female_Parent_",       k)]] <- best_pairs$Female_Parent[k]
          result_row[[base::paste0("Mendelian_Error_Pct_", k)]] <- min_mismatch_val
          result_row[[base::paste0("Markers_Tested_",      k)]] <- comparisons[best_indices[k]]
        }
        result_row[["Assignment_Status"]] <- a_status
      }

      data.table::as.data.table(result_row)
    })

    final_df <- data.table::rbindlist(results_list, fill = TRUE)
  }

  #### Summary ####
  if (verbose) {
    total     <- base::nrow(final_df)
    a_counts  <- base::table(final_df$Assignment_Status)
    base::cat("\n--- Parentage Assignment Summary ---\n")
    base::cat("Total progeny evaluated:", total, "\n")
    for (s in base::names(a_counts))
      base::cat(base::sprintf("  %-14s: %d (%.1f%%)\n", s,
                              a_counts[s], (a_counts[s] / total) * 100))
    base::cat("Min markers threshold :", min_markers,      "\n")
    base::cat("Error threshold       :", error_threshold,  "%\n\n")
  }

  #### Output ####
  if (write_txt) {
    output_filename <- "parentage_results_dt.txt"
    tryCatch({
      data.table::fwrite(final_df, file = output_filename, sep = "\t", quote = FALSE)
      if (verbose) base::cat("Results successfully written to:", output_filename, "\n")
    }, error = function(e) {
      warning("Could not write results to file. Error: ", e$message, call. = FALSE)
    })
  }

  if (verbose) {
    base::cat("\n--- Parentage Assignment Results ---\n")
    base::print(final_df)
    return(base::invisible(final_df))
  } else {
    return(final_df)
  }
}
