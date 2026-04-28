#' Validate Pedigree Trios Using Mendelian Error Analysis
#'
#' Validates parent-offspring trios by calculating Mendelian error rates from
#' SNP genotype data. Identifies incorrect parentage assignments and suggests
#' best-matching replacements. If a list of founders is supplied, trios that
#' are declared founders (both parents coded as 0) are preserved unchanged
#' with no recommendations. Trios removed due to missing genotype data are
#' retained in the output with a NO_GENOTYPE_DATA status.
#'
#' @param pedigree_file Character. Path to the pedigree file (TSV/CSV/TXT)
#'   with columns: ID, Male_Parent, Female_Parent.
#' @param genotypes_file Character. Path to the genotypes file (TSV/CSV/TXT)
#'   with an ID column followed by marker columns coded as 0, 1, 2.
#' @param founders_file Character, optional. Path to a one-column file
#'   listing the IDs of founder individuals. Founders with both parents
#'   coded as 0 are left unchanged with no recommendations. Defaults to NULL.
#' @param trio_error_threshold Numeric. Maximum Mendelian error percentage
#'   to classify a trio as PASS (default: 5.0). Must be between 0 and 100.
#' @param min_markers Integer. Minimum number of non-missing markers
#'   required to evaluate a trio (default: 10).
#' @param single_parent_error_threshold Numeric. Maximum homozygous-marker
#'   mismatch percentage for a parent to be considered acceptable during
#'   parent-level evaluation (default: 2.0). Must be between 0 and 100.
#' @param verbose Logical. If TRUE, prints progress messages, a summary
#'   table, and results to the console (default: TRUE).
#' @param write_results Logical. If TRUE, writes the validation results to
#'   output_filename and saves the plot as a .jpg file (default: TRUE).
#' @param output_filename Character. Path/name of the output results file
#'   (default: "__validation_report.txt"). The plot will be saved using the
#'   same base name with a .jpg extension.
#' @param plot_results Logical. If TRUE, prints a histogram of Trio Mendelian
#'   Error percentages with a threshold line (default: TRUE). If write_results
#'   is also TRUE, the plot is additionally saved as a .jpg file.
#' @param na_string Character. String for missing values in the output file.
#'   Use \code{"NA"} or \code{""} (default: \code{"NA"}).
#' @param corrected_pedigree_filename Character. Path/name of the output file
#'   for the corrected pedigree (default: "corrected_pedigree.txt"). Set to
#'   NULL to suppress writing the corrected pedigree.
#'
#' @return A data.table (returned invisibly) with one row per trio and
#'   the following columns:
#'   \describe{
#'     \item{ID}{Individual ID.}
#'     \item{Orig_Male_Parent}{Declared male parent ID.}
#'     \item{Orig_Female_Parent}{Declared female parent ID.}
#'     \item{Trio_Mendelian_Error_Pct}{Trio-level Mendelian error percentage.}
#'     \item{Trio_Markers_Tested}{Number of markers with non-missing genotypes.}
#'     \item{Status}{One of PASS, FAIL, LOW_MARKERS, NO_DATA, FOUNDERS,
#'           MISSING_MALE_PARENT, MISSING_FEMALE_PARENT, MISSING_BOTH_PARENTS,
#'           or NO_GENOTYPE_DATA.}
#'     \item{Recommended_Correction}{One of NONE, KEEP_BOTH,
#'           REMOVE_MALE_PARENT, REMOVE_FEMALE_PARENT, REMOVE_BOTH,
#'           LOW_MARKERS_KEEP_BOTH, LOW_MARKERS_REMOVE_MALE_PARENT,
#'           LOW_MARKERS_REMOVE_FEMALE_PARENT, LOW_MARKERS_REMOVE_BOTH.}
#'     \item{Male_Parent_Hom_Error_Pct}{Male parent homozygous-marker mismatch percentage.}
#'     \item{Female_Parent_Hom_Error_Pct}{Female parent homozygous-marker mismatch percentage.}
#'     \item{Best_Male_Candidate}{Best-matching male parent candidate ID.}
#'     \item{Best_Male_Candidate_Error_Pct}{Homozygous mismatch percentage for the best male candidate.}
#'     \item{Best_Female_Candidate}{Best-matching female parent candidate ID.}
#'     \item{Best_Female_Candidate_Error_Pct}{Homozygous mismatch percentage for the best female candidate.}
#'   }
#' @export
#'
#' @author Josué Chinchilla-Vargas
#'
#' @importFrom data.table fread fwrite copy data.table set rbindlist
validate_pedigree <- function(pedigree_file, genotypes_file,
                              founders_file = NULL,
                              trio_error_threshold = 5.0,
                              min_markers = 10,
                              single_parent_error_threshold = 2.0,
                              verbose = TRUE,
                              write_results = TRUE,
                              plot_results = TRUE,
                              na_string = "NA",
                              output_filename = "__validation_report.txt",
                              corrected_pedigree_filename = "corrected_pedigree.txt") {

  #### Input validation ####
  if (trio_error_threshold < 0 || trio_error_threshold > 100)
    stop("trio_error_threshold must be between 0 and 100")
  if (single_parent_error_threshold < 0 || single_parent_error_threshold > 100)
    stop("single_parent_error_threshold must be between 0 and 100")
  if (!na_string %in% c("NA", ""))
    stop("na_string must be either 'NA' or ''.")

  tryCatch({
    pedigree <- data.table::fread(pedigree_file, sep = "auto", colClasses = "character")
    genos    <- data.table::fread(genotypes_file, sep = "auto")
  }, error = function(e) {
    stop("Error reading input files. Ensure paths are correct and files are TXT/TSV/CSV.")
  })

  #### Check required columns ####
  required_ped_cols <- c("ID", "Male_Parent", "Female_Parent")
  missing_cols <- base::setdiff(required_ped_cols, base::names(pedigree))
  if (base::length(missing_cols) > 0)
    stop("Pedigree file missing required columns: ",
         base::paste(missing_cols, collapse = ", "))
  if (!"ID" %in% base::names(genos))
    stop("Genotypes file must have an 'ID' column")

  pedigree[, Male_Parent   := as.character(Male_Parent)]
  pedigree[, Female_Parent := as.character(Female_Parent)]
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
  genos_mat <- base::as.matrix(genos, rownames = "ID")
  genos_hom   <- data.table::copy(genos)
  marker_cols <- base::setdiff(base::names(genos_hom), "ID")
  for (col in marker_cols)
    genos_hom[base::get(col) == 1, (col) := NA_integer_]
  genos_hom_mat <- base::as.matrix(genos_hom, rownames = "ID")

  #### Identify trios missing from the genotype file ####
  valid_ids <- as.character(genos$ID)
  has_geno <- pedigree[ID %in% valid_ids &
                         (Male_Parent %in% valid_ids   | Male_Parent == "0") &
                         (Female_Parent %in% valid_ids | Female_Parent == "0")]
  no_geno_rows <- pedigree[!(ID %in% valid_ids) |
                             (!(Male_Parent %in% valid_ids)   & Male_Parent != "0") |
                             (!(Female_Parent %in% valid_ids) & Female_Parent != "0")]
  if (base::nrow(no_geno_rows) > 0 && verbose)
    base::cat("Found", base::nrow(no_geno_rows),
              "trios with missing genotype data; flagged as NO_GENOTYPE_DATA.\n")
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
    prog_id          <- pedigree$ID[i]
    male_parent_id   <- pedigree$Male_Parent[i]
    female_parent_id <- pedigree$Female_Parent[i]
    correction_decision     <- "NONE"
    error_pct               <- NA_real_
    status                  <- "NO_DATA"
    markers_tested          <- 0L
    male_parent_error_pct   <- NA_real_
    female_parent_error_pct <- NA_real_
    best_male_parent        <- NA_character_
    best_male_parent_pct    <- NA_real_
    best_female_parent      <- NA_character_
    best_female_parent_pct  <- NA_real_

    if (male_parent_id == "0" && female_parent_id == "0" &&
        prog_id %in% founder_ids) {
      status              <- "FOUNDERS"
      correction_decision <- "NONE"
    } else {
      if (male_parent_id == "0" && female_parent_id == "0") {
        status              <- "MISSING_BOTH_PARENTS"
        correction_decision <- "NONE"
        best_m                 <- find_best_parent(prog_id, exclude_ids = character(0))
        best_male_parent       <- best_m$id
        best_male_parent_pct   <- best_m$error_pct
        best_f                 <- find_best_parent(prog_id, exclude_ids = c(best_m$id))
        best_female_parent     <- best_f$id
        best_female_parent_pct <- best_f$error_pct
      } else if (male_parent_id == "0" && female_parent_id != "0") {
        status              <- "MISSING_MALE_PARENT"
        correction_decision <- "NONE"
        best_m               <- find_best_parent(prog_id, exclude_ids = c(female_parent_id))
        best_male_parent     <- best_m$id
        best_male_parent_pct <- best_m$error_pct
      } else if (male_parent_id != "0" && female_parent_id == "0") {
        status              <- "MISSING_FEMALE_PARENT"
        correction_decision <- "NONE"
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
          status              <- "NO_DATA"
          correction_decision <- "NONE"
        } else {
          error_pct <- (mismatches / markers_tested) * 100
          if (markers_tested < min_markers) {
            status <- "LOW_MARKERS"
          } else if (error_pct <= trio_error_threshold) {
            status              <- "PASS"
            correction_decision <- "NONE"
          } else {
            status <- "FAIL"
          }
          if (status %in% c("FAIL", "LOW_MARKERS")) {
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
            male_acceptable <- !is.na(male_parent_error_pct) &&
              male_parent_error_pct <= single_parent_error_threshold
            female_acceptable <- !is.na(female_parent_error_pct) &&
              female_parent_error_pct <= single_parent_error_threshold
            if (male_acceptable && female_acceptable) {
              correction_decision <- "KEEP_BOTH"
            } else if (male_acceptable && !female_acceptable) {
              correction_decision    <- "REMOVE_FEMALE_PARENT"
              best_f                 <- find_best_parent(prog_id, exclude_ids = c(male_parent_id))
              best_female_parent     <- best_f$id
              best_female_parent_pct <- best_f$error_pct
            } else if (!male_acceptable && female_acceptable) {
              correction_decision  <- "REMOVE_MALE_PARENT"
              best_m               <- find_best_parent(prog_id, exclude_ids = c(female_parent_id))
              best_male_parent     <- best_m$id
              best_male_parent_pct <- best_m$error_pct
            } else {
              correction_decision    <- "REMOVE_BOTH"
              best_m                 <- find_best_parent(prog_id, exclude_ids = character(0))
              best_male_parent       <- best_m$id
              best_male_parent_pct   <- best_m$error_pct
              best_f                 <- find_best_parent(prog_id, exclude_ids = c(best_m$id))
              best_female_parent     <- best_f$id
              best_female_parent_pct <- best_f$error_pct
            }
            if (status == "LOW_MARKERS")
              correction_decision <- paste0("LOW_MARKERS_", correction_decision)
          }
        }
      }
    }

    data.table::data.table(
      ID                              = prog_id,
      Orig_Male_Parent                = male_parent_id,
      Orig_Female_Parent              = female_parent_id,
      Trio_Mendelian_Error_Pct        = base::round(error_pct, 2),
      Trio_Markers_Tested             = markers_tested,
      Status                          = status,
      Recommended_Correction          = correction_decision,
      Male_Parent_Hom_Error_Pct       = male_parent_error_pct,
      Female_Parent_Hom_Error_Pct     = female_parent_error_pct,
      Best_Male_Candidate             = best_male_parent,
      Best_Male_Candidate_Error_Pct   = best_male_parent_pct,
      Best_Female_Candidate           = best_female_parent,
      Best_Female_Candidate_Error_Pct = best_female_parent_pct
    )
  })

  final_df <- data.table::rbindlist(results_list)

  #### Append NO_GENOTYPE_DATA rows ####
  if (base::nrow(no_geno_rows) > 0) {
    no_geno_df <- data.table::data.table(
      ID                              = no_geno_rows$ID,
      Orig_Male_Parent                = no_geno_rows$Male_Parent,
      Orig_Female_Parent              = no_geno_rows$Female_Parent,
      Trio_Mendelian_Error_Pct        = NA_real_,
      Trio_Markers_Tested             = 0L,
      Status                          = "NO_GENOTYPE_DATA",
      Recommended_Correction          = "NONE",
      Male_Parent_Hom_Error_Pct       = NA_real_,
      Female_Parent_Hom_Error_Pct     = NA_real_,
      Best_Male_Candidate             = NA_character_,
      Best_Male_Candidate_Error_Pct   = NA_real_,
      Best_Female_Candidate           = NA_character_,
      Best_Female_Candidate_Error_Pct = NA_real_
    )
    final_df <- data.table::rbindlist(list(final_df, no_geno_df))
  }

  #### Write corrected pedigree ####

    # pass NULL to suppress writing.
  if (!is.null(corrected_pedigree_filename)) {
    corrected_pedigree <- data.table::copy(original_pedigree)
    for (i in base::seq_len(base::nrow(final_df))) {
      prog_id  <- final_df$ID[i]
      decision <- final_df$Recommended_Correction[i]
      row_idx  <- base::which(corrected_pedigree$ID == prog_id)
      if (decision == "REMOVE_MALE_PARENT") {
        data.table::set(corrected_pedigree, row_idx, "Male_Parent", "0")
      } else if (decision == "REMOVE_FEMALE_PARENT") {
        data.table::set(corrected_pedigree, row_idx, "Female_Parent", "0")
      } else if (decision == "REMOVE_BOTH") {
        data.table::set(corrected_pedigree, row_idx, "Male_Parent", "0")
        data.table::set(corrected_pedigree, row_idx, "Female_Parent", "0")
      }
    }
    tryCatch({
      data.table::fwrite(corrected_pedigree, file = corrected_pedigree_filename,
                         sep = "\t", quote = FALSE)
      if (verbose) base::cat("Corrected pedigree written to:", corrected_pedigree_filename, "\n")
    }, error = function(e) {
      warning("Could not write corrected pedigree. Error: ", e$message, call. = FALSE)
    })
  }

  #### Summary output ####
  if (verbose) {
    total_trios   <- base::nrow(final_df)
    status_counts <- base::table(final_df$Status)
    base::cat("\n--- Trio Validation Summary ---\n")
    base::cat("Total trios in pedigree:", total_trios, "\n")
    for (s in base::names(status_counts))
      base::cat(base::sprintf("%-24s: %d (%.1f%%)\n", s,
                              status_counts[s],
                              (status_counts[s] / total_trios) * 100))
    base::cat("Error threshold:", trio_error_threshold, "%\n")
    base::cat("Homozygous threshold:", single_parent_error_threshold, "%\n")
    base::cat("Minimum markers required:", min_markers, "\n\n")
    corrections <- base::table(final_df$Recommended_Correction)
    base::cat("Correction summary:\n")
    for (decision in base::names(corrections))
      if (decision != "NONE")
        base::cat("  ", decision, ":", corrections[decision], "\n")
    base::cat("\n")
    base::print(final_df)
  }

  #### Write results ####
  if (write_results) {
    tryCatch({
      data.table::fwrite(final_df, file = output_filename,
                         sep = "\t", quote = FALSE, na = na_string)
      if (verbose) base::cat("Results written to:", output_filename, "\n")
    }, error = function(e) {
      warning("Could not write results. Error: ", e$message, call. = FALSE)
    })
  }

  #### Plot results ####
  if (plot_results) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 is required for plot_results = TRUE. Please install it.", call. = FALSE)
    } else {
      plot_df <- final_df[!is.na(final_df$Trio_Mendelian_Error_Pct)]
      plot_df$Plot_Status <- dplyr::case_when(
        plot_df$Recommended_Correction %in% c("NONE",
                                              "KEEP_BOTH",
                                              "LOW_MARKERS_KEEP_BOTH")              ~ "PASS",
        plot_df$Recommended_Correction %in% c("REMOVE_MALE_PARENT",
                                              "REMOVE_FEMALE_PARENT",
                                              "LOW_MARKERS_REMOVE_MALE_PARENT",
                                              "LOW_MARKERS_REMOVE_FEMALE_PARENT")   ~ "FAIL_ONE_PARENT",
        plot_df$Recommended_Correction %in% c("REMOVE_BOTH",
                                              "LOW_MARKERS_REMOVE_BOTH")            ~ "FAIL_BOTH_PARENTS",
        TRUE                                                                         ~ "OTHER"
      )
      n_total  <- nrow(plot_df)
      n_fail   <- sum(plot_df$Trio_Mendelian_Error_Pct > trio_error_threshold)
      n_pass   <- sum(plot_df$Trio_Mendelian_Error_Pct <= trio_error_threshold)
      pct_lost <- round((n_fail / n_total) * 100, 1)
      threshold_label <- paste0("Mendelian Error Threshold: ", trio_error_threshold, "%  |  ",
                                "Lost: ", n_fail, " trios  |  ",
                                "Kept: ", n_pass, " trios")
      p <- ggplot2::ggplot(plot_df,
                           ggplot2::aes(x = Trio_Mendelian_Error_Pct, fill = Plot_Status)) +
        ggplot2::geom_histogram(binwidth = 1, color = "white", alpha = 0.9) +
        ggplot2::geom_vline(xintercept = trio_error_threshold,
                            linetype = "dashed", color = "black", linewidth = 1) +
        ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 5)) +
        ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 5)) +
        ggplot2::scale_fill_manual(
          values = c(
            "PASS"              = "#339900",
            "FAIL_ONE_PARENT"   = "#F1C40F",
            "FAIL_BOTH_PARENTS" = "#cc3333",
            "OTHER"             = "#BDC3C7"
          ),
          labels = c(
            "PASS"              = "Pass",
            "FAIL_ONE_PARENT"   = "Fail - One Parent",
            "FAIL_BOTH_PARENTS" = "Fail - Both Parents",
            "OTHER"             = "Other"
          )
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
      if (write_results) {
        # Fix #3: use tools::file_path_sans_ext() to reliably build the .jpg
        # filename regardless of whether output_filename has an extension or not.
        plot_filename <- paste0(tools::file_path_sans_ext(output_filename), ".jpg")
        tryCatch({
          ggplot2::ggsave(plot_filename, plot = p,
                          device = "jpeg", width = 10, height = 6, dpi = 300)
          if (verbose) base::cat("Plot saved to:", plot_filename, "\n")
        }, error = function(e) {
          warning("Could not save plot. Error: ", e$message, call. = FALSE)
        })
      }
    }
  }

  return(base::invisible(final_df))
}
