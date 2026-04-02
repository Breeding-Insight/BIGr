#' Perform Parentage Assignment from Genotypic Data
#'
#' @description
#' Assigns parents to progeny based on genetic marker data coded as 0, 1, or 2
#' representing allele counts. The function supports several methods, including
#' finding the best single parent (sire, dam, or best match) based on
#' homozygous loci mismatches, or identifying the best parent pair by
#' minimizing Mendelian inheritance errors.
#'
#' @details
#' The function operates in one of two main modes depending on the `method` argument:
#' \itemize{
#'   \item \strong{Homozygous Mismatch Methods} (`"best.sire"`, `"best.dam"`, `"best.match"`):
#'     These methods work by considering only homozygous loci (coded as 0 or 2).
#'     They calculate the percentage of mismatching homozygous loci between each
#'     progeny and the potential parents. Heterozygous loci (coded as 1) are
#'     ignored in this comparison.
#'   \item \strong{Best Pair Method} (`"best.pair"`): This method evaluates all
#'     possible sire-dam pairs for each progeny. It counts the number of loci
#'     that show a Mendelian error given the parental pair's genotypes (e.g.,
#'     two parents with genotype 0 cannot produce a progeny with genotype 1 or 2).
#'     The pair(s) with the minimum error percentage is/are reported as the best match.
#' }
#'
#' @param genotypes_file A character string. Path to the tab-separated file
#'   containing genotypic data for all individuals. Must have an 'ID' column
#'   followed by marker columns.
#' @param parents_file A character string. Path to the tab-separated file listing
#'   potential parent individuals. Must have an 'ID' column and optionally a
#'   'Sex' column ('M' for male, 'F' for female).
#' @param progeny_file A character string. Path to the tab-separated file listing
#'   the progeny individuals to be analyzed. Must have an 'ID' column.
#' @param method A character string specifying the assignment method. Must be one
#'   of `"best.pair"` (default), `"best.sire"`, `"best.dam"`, or `"best.match"`.
#' @param show.ties A logical value. If `TRUE` (default), all tied best matches
#'   for a progeny are reported in wide format. If `FALSE`, only the first
#'   tied match is reported and a warning is issued.
#' @param allow.selfing A logical value. If `TRUE` (default), an individual can
#'   be assigned as both sire and dam. Only applicable when `method = "best.pair"`.
#' @param verbose A logical value. If `TRUE` (default), progress messages and the
#'   final results table are printed to the console.
#' @param write.txt A logical value. If `TRUE` (default), the results table is
#'   written to a file named "parentage_results.txt" in the current directory.
#'
#' @return A `tibble` (data frame) containing the parentage assignment results.
#'   If `verbose = TRUE`, the function prints the results to the console and
#'   returns the `tibble` invisibly.
#'
#' @importFrom dplyr filter mutate across bind_rows
#' @importFrom tibble column_to_rownames tibble as_tibble
#'
#' @export
#'
find_parentage <- function(genotypes_file, parents_file, progeny_file,
                           method = "best.pair",
                           show.ties = TRUE,
                           allow.selfing = TRUE,
                           verbose = TRUE,
                           write.txt = TRUE) {
  #### Input Validation and Data Loading ####
  allowed_methods <- c("best.sire", "best.dam", "best.match", "best.pair")
  if (!method %in% allowed_methods) {
    stop("Method must be one of: ", paste(allowed_methods, collapse = ", "))
  }

  tryCatch({
    genos <- read.table(genotypes_file, header = TRUE, stringsAsFactors = FALSE)
    all_parents <- read.table(parents_file, header = TRUE, stringsAsFactors = FALSE)
    progeny_candidates <- read.table(progeny_file, header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop("Error reading input files. Please ensure paths are correct and files are properly formatted.")
  })

  valid_ids <- genos$ID
  removed_parents <- setdiff(all_parents$ID, valid_ids)
  if (length(removed_parents) > 0) {
    warning("The following parent IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_parents, collapse = ", "), call. = FALSE)
    all_parents <- all_parents %>% filter(ID %in% valid_ids)
  }

  removed_progeny <- setdiff(progeny_candidates$ID, valid_ids)
  if (length(removed_progeny) > 0) {
    warning("The following progeny IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_progeny, collapse = ", "), call. = FALSE)
    progeny_candidates <- progeny_candidates %>% filter(ID %in% valid_ids)
  }

  if (!"Sex" %in% colnames(all_parents)) {
    warning("No 'Sex' column in parents file. All parents treated as ambiguous ('A').")
    all_parents$Sex <- "A"
  }

  all_parents <- all_parents %>% mutate(Sex = toupper(Sex))
  sire_candidates <- all_parents %>% filter(Sex %in% c("M", "A"))
  dam_candidates <- all_parents %>% filter(Sex %in% c("F", "A"))

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
    genos_hom <- genos %>% mutate(across(-ID, ~ ifelse(.x == 1, NA, .x)))
    parent_ids <- switch(method,
                         "best.sire" = sire_candidates$ID,
                         "best.dam" = dam_candidates$ID,
                         "best.match" = union(sire_candidates$ID, dam_candidates$ID)
    )

    parent_genos <- genos_hom %>% filter(ID %in% parent_ids) %>% column_to_rownames("ID")
    progeny_genos <- genos_hom %>% filter(ID %in% progeny_candidates$ID) %>% column_to_rownames("ID")
    results_list <- list()

    for (i in seq_len(nrow(progeny_genos))) {
      progeny_id <- rownames(progeny_genos)[i]
      progeny_vec <- as.numeric(progeny_genos[i, ])

      mismatches <- rowSums(parent_genos != progeny_vec, na.rm = TRUE)
      comparisons <- rowSums(!is.na(parent_genos) & !is.na(progeny_vec))
      percent_mismatch <- (mismatches / comparisons) * 100
      percent_mismatch[is.nan(percent_mismatch)] <- NA

      best_idx <- which.min(percent_mismatch)
      if (length(best_idx) == 0) {
        best_parent_id <- NA
        min_percent <- NA
        markers_tested <- NA
      } else {
        best_parent_id <- rownames(parent_genos)[best_idx]
        min_percent <- percent_mismatch[best_idx]
        markers_tested <- comparisons[best_idx]
      }

      results_list[[progeny_id]] <- tibble(
        Progeny = progeny_id,
        Best_Match = best_parent_id,
        Mendelian_Error_Pct = round(min_percent,2),
        Markers_Tested = markers_tested
      )
    }
    final_df <- bind_rows(results_list)
  }

  #### Logic for Best Pair Method ####
  if (method == "best.pair") {
    genos_mat <- genos %>% column_to_rownames("ID") %>% as.matrix()
    parent_pairs <- expand.grid(Sire = sire_candidates$ID, Dam = dam_candidates$ID, stringsAsFactors = FALSE)

    if (!allow.selfing) {
      parent_pairs <- parent_pairs %>% filter(Sire != Dam)
      if(verbose) cat("Selfing is disallowed. Pairs with identical parents are removed.\n")
    }
    if (nrow(parent_pairs) == 0) stop("No valid parent pairs to test.")

    sire_genos_mat <- genos_mat[parent_pairs$Sire, , drop = FALSE]
    dam_genos_mat <- genos_mat[parent_pairs$Dam, , drop = FALSE]
    results_list <- list()

    for (prog_id in progeny_candidates$ID) {
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
        results_list[[prog_id]] <- tibble(Progeny = prog_id, Markers_Tested = 0)
        next
      }

      best_indices <- which(percent_mismatch == min_mismatch_val)
      best_pairs <- parent_pairs[best_indices, ]

      if (!show.ties && nrow(best_pairs) > 1) {
        warning("Progeny '", prog_id, "' has ", nrow(best_pairs), " tied best pairs. Only one is reported as show.ties=FALSE.", call. = FALSE)
      }

      num_to_report <- if (show.ties) nrow(best_pairs) else 1
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
      results_list[[prog_id]] <- as_tibble(result_row)
    }
    final_df <- bind_rows(results_list)
  }

  #### Output ####
  # write .txt output
  if (write.txt) {
    output_filename <- "parentage_results.txt"
    tryCatch({
      write.table(final_df,
                  file = output_filename,
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
      if (verbose) {
        cat("\nResults successfully written to:", output_filename, "\n")
      }
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

