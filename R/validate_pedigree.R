#' Validate Pedigree Trios Using Mendelian Error Analysis
#'
#' Validates parent-offspring trios against SNP genotype data using Mendelian
#' error rates. Identifies incorrect parentage assignments, suggests
#' best-matching replacements, and outputs a corrected pedigree. Founder trios
#' (both parents coded as 0) are preserved unchanged if a founders file is
#' supplied. Trios absent from the genotype file are retained as
#' no_genotype_data.
#'
#' @param pedigree_file Character. Path to the pedigree file (TSV/CSV/TXT)
#'   with columns: id, male_parent, female_parent.
#' @param genotypes_file Character. Path to the genotypes file (TSV/CSV/TXT)
#'   with an id column followed by marker columns coded as 0, 1, 2.
#' @param founders_file Character, optional. Path to a one-column file listing
#'   founder IDs. Founders with both parents coded as 0 are left unchanged.
#'   Defaults to NULL.
#' @param trio_error_threshold Numeric. Maximum Mendelian error percentage to
#'   classify a trio as pass (default: 5.0). Must be between 0 and 100.
#' @param min_markers Integer. Minimum non-missing markers required to evaluate
#'   a trio (default: 10).
#' @param single_parent_error_threshold Numeric. Maximum homozygous-marker
#'   mismatch percentage for a parent to be considered acceptable
#'   (default: 2.0). Must be between 0 and 100.
#' @param verbose Logical. If TRUE, prints progress, summary, and results to
#'   the console (default: TRUE).
#' @param plot_results Logical. If TRUE, prints a histogram of trio Mendelian
#'   error percentages with a threshold line (default: TRUE).
#'
#' @return An invisible named list with the following elements:
#' \describe{
#'   \item{pass}{Trios that passed the Mendelian error threshold.}
#'   \item{fail}{Trios that failed the Mendelian error threshold.}
#'   \item{low_markers}{Trios with insufficient markers for evaluation.}
#'   \item{no_genotype_data}{Trios absent from the genotype file.}
#'   \item{founders}{Trios identified as founders.}
#'   \item{missing_parents}{Trios with one or both parents coded as 0 (non-founders).}
#'   \item{full_results}{Complete data.table with all trios and all output columns.}
#'   \item{corrected_pedigree}{Pedigree table after applying recommended corrections.}
#'   \item{plot}{ggplot object. Use ggsave() to save if desired.}
#' }
#'
#' @author Josue Chinchilla-Vargas
#'
#' @importFrom data.table fread copy data.table set rbindlist
#' @export
validate_pedigree <- function(pedigree_file, genotypes_file,
                              founders_file                 = NULL,
                              trio_error_threshold          = 5.0,
                              min_markers                   = 10,
                              single_parent_error_threshold = 2.0,
                              verbose                       = TRUE,
                              plot_results                  = TRUE) {

  #### Input validation ####
  if (trio_error_threshold < 0 || trio_error_threshold > 100)
    stop("trio_error_threshold must be between 0 and 100")
  if (single_parent_error_threshold < 0 || single_parent_error_threshold > 100)
    stop("single_parent_error_threshold must be between 0 and 100")

  tryCatch({
    pedigree <- data.table::fread(pedigree_file, sep = "auto", colClasses = "character")
    genos    <- data.table::fread(genotypes_file, sep = "auto")
  }, error = function(e) {
    stop("Error reading input files. Ensure paths are correct and files are TXT/TSV/CSV.")
  })

  #### Check required columns ####
  required_ped_cols <- c("id", "male_parent", "female_parent")
  missing_cols <- base::setdiff(required_ped_cols, base::names(pedigree))
  if (base::length(missing_cols) > 0)
    stop("Pedigree file missing required columns: ",
         base::paste(missing_cols, collapse = ", "))
  if (!"id" %in% base::names(genos))
    stop("Genotypes file must have an 'id' column")

  pedigree[, male_parent   := as.character(male_parent)]
  pedigree[, female_parent := as.character(female_parent)]
  original_pedigree <- data.table::copy(pedigree)

  #### Read founders list ####
  if (!is.null(founders_file)) {
    founders_raw <- tryCatch({
      data.table::fread(founders_file, header = FALSE, colClasses = "character")
    }, error = function(e) {
      stop("Could not read founders list. Ensure it is a plain text or CSV/TSV file.")
    })
    founder_ids <- unique(founders_raw[[1]])
  } else {
    founder_ids <- character(0)
  }

  #### Build genotype matrices ####
  genos_mat   <- base::as.matrix(genos, rownames = "id")
  genos_hom   <- data.table::copy(genos)
  marker_cols <- base::setdiff(base::names(genos_hom), "id")
  for (col in marker_cols)
    genos_hom[base::get(col) == 1, (col) := NA_integer_]
  genos_hom_mat <- base::as.matrix(genos_hom, rownames = "id")

  #### Identify trios missing from the genotype file ####
  valid_ids <- as.character(genos$id)
  has_geno <- pedigree[id %in% valid_ids &
                         (male_parent %in% valid_ids   | male_parent == "0") &
                         (female_parent %in% valid_ids | female_parent == "0")]
  no_geno_rows <- pedigree[!(id %in% valid_ids) |
                             (!(male_parent %in% valid_ids)   & male_parent != "0") |
                             (!(female_parent %in% valid_ids) & female_parent != "0")]
  if (base::nrow(no_geno_rows) > 0 && verbose)
    base::cat("Found", base::nrow(no_geno_rows),
              "trios with missing genotype data; flagged as no_genotype_data.\n")
  pedigree <- has_geno
  if (base::nrow(pedigree) == 0)
    stop("No valid trios remain after filtering for genotype availability.")

  #### Find best matching parent via homozygous mismatch ####
  find_best_parent <- function(prog_id, exclude_ids = base::character(0)) {
    candidates <- base::setdiff(base::rownames(genos_hom_mat),
                                c(prog_id, exclude_ids))
    if (base::length(candidates) == 0)
      return(base::list(id = NA_character_, error_pct = NA_real_))
    prog_hom <- genos_hom_mat[prog_id, ]
    errors <- base::sapply(candidates, function(cand_id) {
      cand_hom    <- genos_hom_mat[cand_id, ]
      comparisons <- base::sum(!base::is.na(cand_hom) & !base::is.na(prog_hom))
      if (comparisons == 0) return(NA_real_)
      (base::sum(cand_hom != prog_hom, na.rm = TRUE) / comparisons) * 100
    })
    if (base::all(base::is.na(errors)))
      return(base::list(id = NA_character_, error_pct = NA_real_))
    best_idx <- base::which.min(errors)
    base::list(id = candidates[best_idx],
               error_pct = base::round(errors[best_idx], 2))
  }

  #### Main trio evaluation loop ####
  results_list <- base::lapply(base::seq_len(base::nrow(pedigree)), function(i) {
    prog_id          <- pedigree$id[i]
    male_parent_id   <- pedigree$male_parent[i]
    female_parent_id <- pedigree$female_parent[i]
    correction_decision     <- "none"
    error_pct               <- NA_real_
    status                  <- "no_data"
    markers_tested          <- 0L
    male_parent_error_pct   <- NA_real_
    female_parent_error_pct <- NA_real_
    best_male_parent        <- NA_character_
    best_male_parent_pct    <- NA_real_
    best_female_parent      <- NA_character_
    best_female_parent_pct  <- NA_real_

    if (male_parent_id == "0" && female_parent_id == "0" &&
        prog_id %in% founder_ids) {
      status              <- "founders"
      correction_decision <- "none"
    } else {
      if (male_parent_id == "0" && female_parent_id == "0") {
        status              <- "missing_both_parents"
        correction_decision <- "none"
        best_m                 <- find_best_parent(prog_id, exclude_ids = character(0))
        best_male_parent       <- best_m$id
        best_male_parent_pct   <- best_m$error_pct
        best_f                 <- find_best_parent(prog_id, exclude_ids = c(best_m$id))
        best_female_parent     <- best_f$id
        best_female_parent_pct <- best_f$error_pct
      } else if (male_parent_id == "0" && female_parent_id != "0") {
        status              <- "missing_male_parent"
        correction_decision <- "none"
        best_m               <- find_best_parent(prog_id, exclude_ids = c(female_parent_id))
        best_male_parent     <- best_m$id
        best_male_parent_pct <- best_m$error_pct
      } else if (male_parent_id != "0" && female_parent_id == "0") {
        status              <- "missing_female_parent"
        correction_decision <- "none"
        best_f                 <- find_best_parent(prog_id, exclude_ids = c(male_parent_id))
        best_female_parent     <- best_f$id
        best_female_parent_pct <- best_f$error_pct
      } else {
        progeny_vec       <- genos_mat[prog_id, ]
        male_parent_vec   <- genos_mat[male_parent_id, ]
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
        if (markers_tested == 0) {
          status              <- "no_data"
          correction_decision <- "none"
        } else {
          error_pct <- (mismatches / markers_tested) * 100
          if (markers_tested < min_markers) {
            status <- "low_markers"
          } else if (error_pct <= trio_error_threshold) {
            status              <- "pass"
            correction_decision <- "none"
          } else {
            status <- "fail"
          }
          if (status %in% c("fail", "low_markers")) {
            progeny_hom       <- genos_hom_mat[prog_id, ]
            male_parent_hom   <- genos_hom_mat[male_parent_id, ]
            female_parent_hom <- genos_hom_mat[female_parent_id, ]
            male_comparisons <- base::sum(!base::is.na(male_parent_hom) &
                                            !base::is.na(progeny_hom))
            male_parent_error_pct <- if (male_comparisons == 0) NA_real_ else
              base::round((base::sum(male_parent_hom != progeny_hom, na.rm = TRUE) /
                             male_comparisons) * 100, 2)
            female_comparisons <- base::sum(!base::is.na(female_parent_hom) &
                                              !base::is.na(progeny_hom))
            female_parent_error_pct <- if (female_comparisons == 0) NA_real_ else
              base::round((base::sum(female_parent_hom != progeny_hom, na.rm = TRUE) /
                             female_comparisons) * 100, 2)
            male_acceptable   <- !is.na(male_parent_error_pct) &&
              male_parent_error_pct <= single_parent_error_threshold
            female_acceptable <- !is.na(female_parent_error_pct) &&
              female_parent_error_pct <= single_parent_error_threshold
            if (male_acceptable && female_acceptable) {
              correction_decision <- "keep_both"
            } else if (male_acceptable && !female_acceptable) {
              correction_decision    <- "remove_female_parent"
              best_f                 <- find_best_parent(prog_id, exclude_ids = c(male_parent_id))
              best_female_parent     <- best_f$id
              best_female_parent_pct <- best_f$error_pct
            } else if (!male_acceptable && female_acceptable) {
              correction_decision  <- "remove_male_parent"
              best_m               <- find_best_parent(prog_id, exclude_ids = c(female_parent_id))
              best_male_parent     <- best_m$id
              best_male_parent_pct <- best_m$error_pct
            } else {
              correction_decision    <- "remove_both"
              best_m                 <- find_best_parent(prog_id, exclude_ids = character(0))
              best_male_parent       <- best_m$id
              best_male_parent_pct   <- best_m$error_pct
              best_f                 <- find_best_parent(prog_id, exclude_ids = c(best_m$id))
              best_female_parent     <- best_f$id
              best_female_parent_pct <- best_f$error_pct
            }
            if (status == "low_markers")
              correction_decision <- paste0("low_markers_", correction_decision)
          }
        }
      }
    }

    data.table::data.table(
      id                              = prog_id,
      orig_male_parent                = male_parent_id,
      orig_female_parent              = female_parent_id,
      trio_mendelian_error_pct        = base::round(error_pct, 2),
      trio_markers_tested             = markers_tested,
      status                          = status,
      recommended_correction          = correction_decision,
      male_parent_hom_error_pct       = male_parent_error_pct,
      female_parent_hom_error_pct     = female_parent_error_pct,
      best_male_candidate             = best_male_parent,
      best_male_candidate_error_pct   = best_male_parent_pct,
      best_female_candidate           = best_female_parent,
      best_female_candidate_error_pct = best_female_parent_pct
    )
  })

  final_df <- data.table::rbindlist(results_list)

  #### Append no_genotype_data rows ####
  if (base::nrow(no_geno_rows) > 0) {
    no_geno_df <- data.table::data.table(
      id                              = no_geno_rows$id,
      orig_male_parent                = no_geno_rows$male_parent,
      orig_female_parent              = no_geno_rows$female_parent,
      trio_mendelian_error_pct        = NA_real_,
      trio_markers_tested             = 0L,
      status                          = "no_genotype_data",
      recommended_correction          = "none",
      male_parent_hom_error_pct       = NA_real_,
      female_parent_hom_error_pct     = NA_real_,
      best_male_candidate             = NA_character_,
      best_male_candidate_error_pct   = NA_real_,
      best_female_candidate           = NA_character_,
      best_female_candidate_error_pct = NA_real_
    )
    final_df <- data.table::rbindlist(list(final_df, no_geno_df))
  }

  #### Build corrected pedigree ####
  corrected_pedigree <- data.table::copy(original_pedigree)
  for (i in base::seq_len(base::nrow(final_df))) {
    prog_id  <- final_df$id[i]
    decision <- final_df$recommended_correction[i]
    row_idx  <- base::which(corrected_pedigree$id == prog_id)
    if (decision == "remove_male_parent") {
      data.table::set(corrected_pedigree, row_idx, "male_parent", "0")
    } else if (decision == "remove_female_parent") {
      data.table::set(corrected_pedigree, row_idx, "female_parent", "0")
    } else if (decision %in% c("remove_both",
                               "low_markers_remove_both",
                               "low_markers_remove_male_parent",
                               "low_markers_remove_female_parent")) {
      if (grepl("male", decision)) {
        data.table::set(corrected_pedigree, row_idx, "male_parent", "0")
      }
      if (grepl("female", decision)) {
        data.table::set(corrected_pedigree, row_idx, "female_parent", "0")
      }
      if (decision == "low_markers_remove_both" || decision == "remove_both") {
        data.table::set(corrected_pedigree, row_idx, "male_parent",   "0")
        data.table::set(corrected_pedigree, row_idx, "female_parent", "0")
      }
    }
  }

  #### Summary output ####
  if (verbose) {
    total_trios   <- base::nrow(final_df)
    status_counts <- base::table(final_df$status)
    base::cat("\n--- Trio Validation Summary ---\n")
    base::cat("Total trios in pedigree:", total_trios, "\n")
    for (s in base::names(status_counts))
      base::cat(base::sprintf("%-28s: %d (%.1f%%)\n", s,
                              status_counts[s],
                              (status_counts[s] / total_trios) * 100))
    base::cat("Error threshold:", trio_error_threshold, "%\n")
    base::cat("Homozygous threshold:", single_parent_error_threshold, "%\n")
    base::cat("Minimum markers required:", min_markers, "\n\n")
    corrections <- base::table(final_df$recommended_correction)
    base::cat("Correction summary:\n")
    for (decision in base::names(corrections))
      if (decision != "none")
        base::cat("  ", decision, ":", corrections[decision], "\n")
    base::cat("\n")
    base::print(final_df)
  }

  #### Plot results ####
  p <- NULL
  if (plot_results) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 is required for plot_results = TRUE. Please install it.", call. = FALSE)
    } else {
      plot_df <- final_df[!is.na(final_df$trio_mendelian_error_pct)]
      plot_df$plot_status <- dplyr::case_when(
        plot_df$recommended_correction %in% c("none",
                                              "keep_both",
                                              "low_markers_keep_both")                    ~ "pass",
        plot_df$recommended_correction %in% c("remove_male_parent",
                                              "remove_female_parent",
                                              "low_markers_remove_male_parent",
                                              "low_markers_remove_female_parent")         ~ "fail_one_parent",
        plot_df$recommended_correction %in% c("remove_both",
                                              "low_markers_remove_both")                  ~ "fail_both_parents",
        TRUE                                                                               ~ "other"
      )
      n_total <- nrow(plot_df)
      n_fail  <- sum(plot_df$trio_mendelian_error_pct > trio_error_threshold)
      n_pass  <- sum(plot_df$trio_mendelian_error_pct <= trio_error_threshold)
      threshold_label <- paste0(
        "Mendelian Error Threshold: ", trio_error_threshold, "%  |  ",
        "Lost: ", n_fail, " trios  |  ",
        "Kept: ", n_pass, " trios"
      )
      p <- ggplot2::ggplot(plot_df,
                           ggplot2::aes(x = trio_mendelian_error_pct, fill = plot_status)) +
        ggplot2::geom_histogram(binwidth = 1, color = "white", alpha = 0.9) +
        ggplot2::geom_vline(xintercept = trio_error_threshold,
                            linetype = "dashed", color = "black", linewidth = 1) +
        ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 5)) +
        ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 5)) +
        ggplot2::scale_fill_manual(
          values = c("pass"             = "#339900",
                     "fail_one_parent"  = "#F1C40F",
                     "fail_both_parents"= "#cc3333",
                     "other"            = "#BDC3C7"),
          labels = c("pass"             = "Pass",
                     "fail_one_parent"  = "Fail - One Parent",
                     "fail_both_parents"= "Fail - Both Parents",
                     "other"            = "Other")
        ) +
        ggplot2::labs(
          title    = "Trio Mendelian Error Distribution",
          subtitle = paste0("Trios with Genotype Data Tested: ", n_total, "\n \n", threshold_label),
          x        = "Mendelian Error (%)",
          y        = "Number of Trios",
          fill     = "Status"
        ) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(legend.position = "top")
      print(p)
    }
  }

  #### Compile and return named list ####
  output_list <- base::list(
    pass               = final_df[status == "pass"],
    fail               = final_df[status == "fail"],
    low_markers        = final_df[status == "low_markers"],
    no_genotype_data   = final_df[status == "no_genotype_data"],
    founders           = final_df[status == "founders"],
    missing_parents    = final_df[status %in% c("missing_both_parents",
                                                "missing_male_parent",
                                                "missing_female_parent")],
    full_results       = final_df,
    corrected_pedigree = corrected_pedigree,
    plot               = p
  )
  return(base::invisible(output_list))
}
