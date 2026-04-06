#' Find Parentage Assignments for Progeny
#'
#' Assigns the most likely parent(s) to each progeny individual based on
#' genotypic data using Mendelian error rates or homozygous mismatch rates.
#'
#' @param genotypes_file Path to a TSV/CSV file containing genotype data.
#'   Must include an 'ID' column followed by marker columns coded as 0, 1, 2
#'   (allele dosage).
#' @param parents_file Path to a TSV/CSV file listing candidate parent IDs.
#'   Must include an 'ID' column. An optional 'Sex' column with values
#'   'M' (sire), 'F' (dam), or 'A' (ambiguous) determines which parents are
#'   tested for each role. If absent, all parents are treated as ambiguous.
#' @param progeny_file Path to a TSV/CSV file listing progeny IDs to assign.
#'   Must include an 'ID' column.
#' @param method Character. Parentage assignment method. One of:
#'   \itemize{
#'     \item \code{"best.sire"}  — finds the best sire for each progeny using
#'       homozygous mismatch rate.
#'     \item \code{"best.dam"}   — finds the best dam for each progeny using
#'       homozygous mismatch rate.
#'     \item \code{"best.match"} — finds the single best parent (either sex)
#'       using homozygous mismatch rate.
#'     \item \code{"best.pair"}  — finds the best sire-dam pair for each
#'       progeny using full Mendelian error rate (default).
#'   }
#' @param show_ties Logical. If \code{TRUE}, all tied best pairs are reported
#'   as additional columns (\code{Sire_1}, \code{Sire_2}, etc.) when
#'   \code{method = "best.pair"}. If \code{FALSE}, only one tied pair is
#'   reported with a warning. Default is \code{TRUE}.
#' @param allow_selfing Logical. If \code{FALSE}, sire-dam pairs where both
#'   IDs are identical are excluded when \code{method = "best.pair"}.
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
#'     \item \code{best.sire} / \code{best.dam} / \code{best.match}: \code{Progeny},
#'       \code{Best_Match}, \code{Mendelian_Error_Pct}, \code{Markers_Tested}.
#'     \item \code{best.pair} (no ties): \code{Progeny}, \code{Sire}, \code{Dam},
#'       \code{Mendelian_Error_Pct}, \code{Markers_Tested}.
#'     \item \code{best.pair} (with ties): columns are suffixed \code{_1}, \code{_2},
#'       etc. for each tied pair.
#'   }
#'   Returned invisibly when \code{verbose = TRUE}.
#'
#' @details
#' For \code{"best.sire"}, \code{"best.dam"}, and \code{"best.match"}, only
#' homozygous markers (coded 0 or 2) are used for comparison; heterozygous
#' markers (coded 1) are set to \code{NA}. This reduces false mismatches caused
#' by phase ambiguity.
#'
#' For \code{"best.pair"}, all markers are used and full Mendelian inheritance
#' rules are applied across all possible sire-dam combinations via
#' \code{data.table::CJ()}.
#'
#' Individuals in \code{parents_file} or \code{progeny_file} that are absent
#' from \code{genotypes_file} are removed with a warning.
#'
#' @examples
#' \dontrun{
#' # Assign best sire-dam pair to each progeny
#' results <- find_parentage(
#'   genotypes_file = "genotypes.txt",
#'   parents_file   = "parents.txt",
#'   progeny_file   = "progeny.txt",
#'   method         = "best.pair",
#'   show_ties      = TRUE,
#'   allow_selfing  = FALSE
#' )
#'
#' # Find best individual parent match (ignoring sex)
#' results <- find_parentage(
#'   genotypes_file = "genotypes.txt",
#'   parents_file   = "parents.txt",
#'   progeny_file   = "progeny.txt",
#'   method         = "best.match"
#' )
#' }
#'
#' @importFrom data.table fread fwrite copy CJ rbindlist set data.table
#' @export
find_parentage <- function(genotypes_file, parents_file, progeny_file,
                           method = "best.pair",
                           show_ties = TRUE,
                           allow_selfing = TRUE,
                           verbose = TRUE,
                           write_txt = TRUE) {

  #### Input Validation and Data Loading ####
  allowed_methods <- c("best.sire", "best.dam", "best.match", "best.pair")
  if (!method %in% allowed_methods) {
    stop("Method must be one of: ", paste(allowed_methods, collapse = ", "))
  }

  tryCatch({
    genos <- fread(genotypes_file)
    all_parents <- fread(parents_file)
    progeny_candidates <- fread(progeny_file)
  }, error = function(e) {
    stop("Error reading input files. Ensure paths are correct and files are TSV/CSV.")
  })

  valid_ids <- genos$ID
  removed_parents <- setdiff(all_parents$ID, valid_ids)
  if (length(removed_parents) > 0) {
    warning("The following parent IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_parents, collapse = ", "), call. = FALSE)
    all_parents <- all_parents[ID %in% valid_ids]
  }

  removed_progeny <- setdiff(progeny_candidates$ID, valid_ids)
  if (length(removed_progeny) > 0) {
    warning("The following progeny IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_progeny, collapse = ", "), call. = FALSE)
    progeny_candidates <- progeny_candidates[ID %in% valid_ids]
  }

  if (!"Sex" %in% colnames(all_parents)) {
    warning("No 'Sex' column in parents file. All parents treated as ambiguous ('A').")
    all_parents[, Sex := "A"]
  }

  all_parents[, Sex := toupper(Sex)]
  sire_candidates <- all_parents[Sex %in% c("M", "A", "NA")]
  dam_candidates <- all_parents[Sex %in% c("F", "A", "NA")]

  if (nrow(sire_candidates) == 0 && method %in% c("best.sire", "best.pair")) {
    warning("No valid sire candidates remain after filtering.", call. = FALSE)
  }
  if (nrow(dam_candidates) == 0 && method %in% c("best.dam", "best.pair")) {
    warning("No valid dam candidates remain after filtering.", call. = FALSE)
  }
  if (nrow(progeny_candidates) == 0) {
    stop("No valid progeny candidates remain after filtering.")
  }

  #### Logic for Homozygous Matching Methods ####
  if (method %in% c("best.sire", "best.dam", "best.match")) {
    genos_hom <- copy(genos)
    marker_cols <- setdiff(names(genos_hom), "ID")
    for (col in marker_cols) {
      genos_hom[get(col) == 1, (col) := NA_integer_]
    }

    parent_ids <- switch(method,
                         "best.sire" = sire_candidates$ID,
                         "best.dam" = dam_candidates$ID,
                         "best.match" = union(sire_candidates$ID, dam_candidates$ID))

    parent_genos <- as.matrix(genos_hom[ID %in% parent_ids], rownames = "ID")
    progeny_genos <- as.matrix(genos_hom[ID %in% progeny_candidates$ID], rownames = "ID")

    results_list <- lapply(rownames(progeny_genos), function(progeny_id) {
      progeny_vec <- progeny_genos[progeny_id, ]
      mismatches <- rowSums(parent_genos != progeny_vec, na.rm = TRUE)
      comparisons <- rowSums(!is.na(parent_genos) & !is.na(progeny_vec))
      percent_mismatch <- (mismatches / comparisons) * 100
      percent_mismatch[is.nan(percent_mismatch)] <- NA

      best_idx <- which.min(percent_mismatch)
      if (length(best_idx) == 0) {
        data.table(Progeny = progeny_id, Best_Match = NA, Mendelian_Error_Pct = NA, Markers_Tested = NA)
      } else {
        data.table(Progeny = progeny_id,
                   Best_Match = rownames(parent_genos)[best_idx],
                   Mendelian_Error_Pct = round(percent_mismatch[best_idx], 2),
                   Markers_Tested = comparisons[best_idx])
      }
    })
    final_df <- rbindlist(results_list)
  }

  #### Logic for Best Pair Method ####
  if (method == "best.pair") {
    genos_mat <- as.matrix(genos, rownames = "ID")

    parent_pairs <- CJ(Sire = sire_candidates$ID, Dam = dam_candidates$ID)

    if (!allow_selfing) {
      parent_pairs <- parent_pairs[Sire != Dam]
      if (verbose) cat("Selfing is disallowed. Pairs with identical parents are removed.\n")
    }
    if (nrow(parent_pairs) == 0) stop("No valid parent pairs to test.")

    sire_genos_mat <- genos_mat[parent_pairs$Sire, , drop = FALSE]
    dam_genos_mat <- genos_mat[parent_pairs$Dam, , drop = FALSE]

    results_list <- lapply(progeny_candidates$ID, function(prog_id) {
      progeny_vec <- genos_mat[prog_id, ]

      mismatches <- rowSums(
        (sire_genos_mat == 0 & dam_genos_mat == 0 & progeny_vec > 0) |
          (sire_genos_mat == 2 & dam_genos_mat == 2 & progeny_vec < 2) |
          ((sire_genos_mat == 0 & dam_genos_mat == 1) | (sire_genos_mat == 1 & dam_genos_mat == 0)) & (progeny_vec == 2) |
          ((sire_genos_mat == 2 & dam_genos_mat == 1) | (sire_genos_mat == 1 & dam_genos_mat == 2)) & (progeny_vec == 0) |
          ((sire_genos_mat == 0 & dam_genos_mat == 2) | (sire_genos_mat == 2 & dam_genos_mat == 0)) & (progeny_vec != 1),
        na.rm = TRUE
      )

      comparisons <- rowSums(!is.na(sire_genos_mat) & !is.na(dam_genos_mat) & !is.na(progeny_vec))
      percent_mismatch <- (mismatches / comparisons) * 100
      percent_mismatch[is.nan(percent_mismatch)] <- NA

      min_mismatch_val <- min(percent_mismatch, na.rm = TRUE)

      if (is.infinite(min_mismatch_val)) {
        return(data.table(Progeny = prog_id, Markers_Tested = 0))
      }

      best_indices <- which(percent_mismatch == min_mismatch_val)
      best_pairs <- parent_pairs[best_indices]

      if (!show_ties && nrow(best_pairs) > 1) {
        warning("Progeny '", prog_id, "' has ", nrow(best_pairs), " tied best pairs. Only one is reported as show_ties=FALSE.", call. = FALSE)
      }

      num_to_report <- if (show_ties) nrow(best_pairs) else 1
      num_to_report <- min(nrow(best_pairs), num_to_report)

      result_row <- list(Progeny = prog_id)
      if (num_to_report == 1) {
        result_row[['Sire']] <- best_pairs$Sire[1]
        result_row[['Dam']] <- best_pairs$Dam[1]
        result_row[['Mendelian_Error_Pct']] <- sprintf("%.2f", min_mismatch_val)
        result_row[['Markers_Tested']] <- comparisons[best_indices[1]]
      } else if (num_to_report > 1) {
        for (k in 1:num_to_report) {
          result_row[[paste0("Sire_", k)]] <- best_pairs$Sire[k]
          result_row[[paste0("Dam_", k)]] <- best_pairs$Dam[k]
          result_row[[paste0("Mendelian_Error_Pct_", k)]] <- min_mismatch_val
          result_row[[paste0("Markers_Tested_", k)]] <- comparisons[best_indices[k]]
        }
      }
      as.data.table(result_row)
    })
    final_df <- rbindlist(results_list, fill = TRUE)
  }

  #### Output ####
  if (write_txt) {
    output_filename <- "parentage_results_dt.txt"
    tryCatch({
      fwrite(final_df, file = output_filename, sep = "\t", quote = FALSE)
      if (verbose) cat("\nResults successfully written to:", output_filename, "\n")
    }, error = function(e) {
      warning("Could not write results to file. Error: ", e$message, call. = FALSE)
    })
  }

  if (verbose) {
    cat("\n--- Parentage Assignment Results ---\n")
    print(final_df)
    return(invisible(final_df))
  } else {
    return(final_df)
  }
}
