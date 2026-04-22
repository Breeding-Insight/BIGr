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
#' @param write_txt Logical. If TRUE, writes validation results to
#'   output_filename (default: TRUE).
#' @param output_filename Character. Path/name of the output file
#'   (default: "pedigree_validation_results.txt").
#'
#' @return A data.table (returned invisibly) with one row per trio and
#'   the following columns:
#'   \describe{
#'     \item{ID}{Individual ID.}
#'     \item{Male_Parent}{Declared male parent ID.}
#'     \item{Female_Parent}{Declared female parent ID.}
#'     \item{Mendelian_Error_Pct}{Trio-level Mendelian error percentage.}
#'     \item{Markers_Tested}{Number of markers with non-missing genotypes.}
#'     \item{Status}{One of PASS, FAIL, LOW_MARKERS, NO_DATA, FOUNDERS,
#'           MISSING_MALE_PARENT, MISSING_FEMALE_PARENT, MISSING_BOTH_PARENTS,
#'           or NO_GENOTYPE_DATA.}
#'     \item{Correction_Decision}{One of NONE, KEEP_BOTH,
#'           REMOVE_MALE_PARENT, REMOVE_FEMALE_PARENT, REMOVE_BOTH,
#'           LOW_MARKERS_KEEP_BOTH, LOW_MARKERS_REMOVE_MALE_PARENT,
#'           LOW_MARKERS_REMOVE_FEMALE_PARENT, LOW_MARKERS_REMOVE_BOTH.}
#'     \item{Male_Parent_Hom_Error_Pct}{Male parent homozygous-marker mismatch percentage.}
#'     \item{Female_Parent_Hom_Error_Pct}{Female parent homozygous-marker mismatch percentage.}
#'     \item{Best_Male_Parent}{Best-matching male parent candidate ID.}
#'     \item{Best_Male_Parent_Error_Pct}{Homozygous mismatch percentage for the best male parent candidate.}
#'     \item{Best_Female_Parent}{Best-matching female parent candidate ID.}
#'     \item{Best_Female_Parent_Error_Pct}{Homozygous mismatch percentage for the best female parent candidate.}
#'   }
#' @export
#' @importFrom data.table fread fwrite copy  data.table set rbindlist
validate_pedigree <- function(pedigree_file, genotypes_file,
                              founders_file = NULL,
                              trio_error_threshold = 5.0,
                              min_markers = 10,
                              single_parent_error_threshold = 2.0,
                              verbose = TRUE,
                              write_txt = TRUE,
                              output_filename = "pedigree_validation_results.txt") {

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
  required_ped_cols <- c("ID", "Male_Parent", "Female_Parent")
  missing_cols <- base::setdiff(required_ped_cols, base::names(pedigree))
  if (base::length(missing_cols) > 0)
    stop("Pedigree file missing required columns: ",
         base::paste(missing_cols, collapse = ", "))
  if (!"ID" %in% base::names(genos))
    stop("Genotypes file must have an 'ID' column")

  # Ensure parent columns are character for consistent "0" comparisons
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

  #### Build genotype matrices ####
  genos_mat <- base::as.matrix(genos, rownames = "ID")

  # Homozygous-only matrix (het markers set to NA)
  genos_hom   <- data.table::copy(genos)
  marker_cols <- base::setdiff(base::names(genos_hom), "ID")
  for (col in marker_cols)
    genos_hom[base::get(col) == 1, (col) := NA_integer_]
  genos_hom_mat <- base::as.matrix(genos_hom, rownames = "ID")

  #### Helper: find best matching parent via homozygous mismatch ####
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

    best_idx <- base::which.min(errors)
    base::list(id = candidates[best_idx],
               error_pct = base::round(errors[best_idx], 2))
  }

  #### Main trio evaluation loop ####
  results_list <- base::lapply(base::seq_len(base::nrow(pedigree)), function(i) {

    prog_id          <- pedigree$ID[i]
    male_parent_id   <- pedigree$Male_Parent[i]
    female_parent_id <- pedigree$Female_Parent[i]

    # Default values
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

    ## Founder check — both parents "0" and ID in founders list
    if (male_parent_id == "0" && female_parent_id == "0" &&
        prog_id %in% founder_ids) {
      status              <- "FOUNDERS"
      correction_decision <- "NONE"

    } else {

      ## Missing parent(s) — recommendations only, "0"s preserved in pedigree
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

        ## Both parents present — Mendelian error calculation
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

          # LOW_MARKERS still computes parent mismatch/recommendations
          if (markers_tested < min_markers) {
            status <- "LOW_MARKERS"
          } else if (error_pct <= trio_error_threshold) {
            status              <- "PASS"
            correction_decision <- "NONE"
          } else {
            status <- "FAIL"
          }

          # Run parent-level evaluation for both FAIL and LOW_MARKERS
          if (status %in% c("FAIL", "LOW_MARKERS")) {

            # Homozygous mismatch per parent
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

            # Do not alter corrected pedigree for LOW_MARKERS rows
            if (status == "LOW_MARKERS")
              correction_decision <- paste0("LOW_MARKERS_", correction_decision)
          }
        }
      }
    }

    data.table::data.table(
      ID                           = prog_id,
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

  #### Append NO_GENOTYPE_DATA rows to the final report ####
  if (base::nrow(no_geno_rows) > 0) {
    no_geno_df <- data.table::data.table(
      ID                           = no_geno_rows$ID,
      Male_Parent                  = no_geno_rows$Male_Parent,
      Female_Parent                = no_geno_rows$Female_Parent,
      Mendelian_Error_Pct          = NA_real_,
      Markers_Tested               = 0L,
      Status                       = "NO_GENOTYPE_DATA",
      Correction_Decision          = "NONE",
      Male_Parent_Hom_Error_Pct    = NA_real_,
      Female_Parent_Hom_Error_Pct  = NA_real_,
      Best_Male_Parent             = NA_character_,
      Best_Male_Parent_Error_Pct   = NA_real_,
      Best_Female_Parent           = NA_character_,
      Best_Female_Parent_Error_Pct = NA_real_
    )
    final_df <- data.table::rbindlist(list(final_df, no_geno_df))
  }

  #### Write corrected pedigree ####
  corrected_pedigree <- data.table::copy(original_pedigree)
  for (i in base::seq_len(base::nrow(final_df))) {
    prog_id  <- final_df$ID[i]
    decision <- final_df$Correction_Decision[i]
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
    data.table::fwrite(corrected_pedigree, file = "corrected_pedigree.txt",
                       sep = "\t", quote = FALSE)
    if (verbose) base::cat("Corrected pedigree written to: corrected_pedigree.txt\n")
  }, error = function(e) {
    warning("Could not write corrected pedigree. Error: ", e$message, call. = FALSE)
  })

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

    corrections <- base::table(final_df$Correction_Decision)
    base::cat("Correction summary:\n")
    for (decision in base::names(corrections))
      if (decision != "NONE")
        base::cat("  ", decision, ":", corrections[decision], "\n")
    base::cat("\n")
    base::print(final_df)
  }

  #### Write results ####
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
