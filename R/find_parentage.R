#' Find Parentage Assignments for Progeny
#'
#' Assigns the most likely parent(s) to each progeny from SNP genotype data
#' using Mendelian error rates or homozygous mismatch rates. Parents or progeny
#' absent from the genotype file are removed with a warning.
#'
#' @param genotypes_file Path to a TSV/CSV/TXT file with an 'id' column
#'   followed by marker columns coded as 0, 1, 2 (allele dosage).
#' @param parents_file Path to a TSV/CSV/TXT file with an 'id' column and an
#'   optional 'sex' column ('M', 'F', or 'A'). If absent, all parents are
#'   treated as ambiguous.
#' @param progeny_file Path to a TSV/CSV/TXT file with an 'id' column.
#' @param method Character. One of \code{"best_male_parent"},
#'   \code{"best_female_parent"}, \code{"best_match"}, or
#'   \code{"best_pair"} (default). See Details.
#' @param min_markers Integer. Minimum markers required; fewer flags
#'   \code{low_markers} (default: \code{10}).
#' @param error_threshold Numeric. Maximum mismatch percentage; exceeded values
#'   flag \code{high_error} (default: \code{5.0}). Must be between 0 and 100.
#' @param show_ties Logical. If \code{TRUE}, tied best pairs are appended as
#'   suffix columns (e.g. \code{male_parent_2}) for \code{"best_pair"}.
#'   Default is \code{TRUE}.
#' @param allow_selfing Logical. If \code{FALSE}, pairs with identical male and
#'   female parent IDs are excluded. Default is \code{TRUE}.
#' @param verbose Logical. If \code{TRUE}, prints progress and summary to the
#'   console. Default is \code{TRUE}.
#' @param plot_results Logical. If \code{TRUE}, plots the Mendelian error
#'   distribution colored by status. Requires \code{ggplot2}. Default is \code{TRUE}.
#'
#' @return A named list (returned invisibly) with the following elements:
#' \describe{
#'   \item{pass}{Progeny with a confident parentage assignment.}
#'   \item{high_error}{Progeny whose best assignment exceeds the error threshold.}
#'   \item{low_markers}{Progeny with insufficient markers for a valid assignment.}
#'   \item{full_results}{Complete data.table with all progeny and all output columns.}
#'   \item{plot}{ggplot object. Use ggsave() to save if desired.}
#' }
#'
#' @author Josue Chinchilla-Vargas
#'
#' @importFrom data.table fread copy CJ rbindlist set data.table as.data.table
#' @export
find_parentage <- function(genotypes_file, parents_file, progeny_file,
                           method          = "best_pair",
                           min_markers     = 10,
                           error_threshold = 5.0,
                           show_ties       = TRUE,
                           allow_selfing   = TRUE,
                           verbose         = TRUE,
                           plot_results    = TRUE) {

  ## silence R CMD check NOTEs
  id <- sex <- male_parent <- female_parent <- NULL
  mendelian_error_pct <- plot_status <- status <- NULL

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

  valid_ids       <- genos$id
  removed_parents <- base::setdiff(all_parents$id, valid_ids)
  if (base::length(removed_parents) > 0) {
    warning("The following parent IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_parents, collapse = ", "), call. = FALSE)
    all_parents <- all_parents[id %in% valid_ids]
  }

  removed_progeny <- base::setdiff(progeny_candidates$id, valid_ids)
  if (base::length(removed_progeny) > 0) {
    warning("The following progeny IDs were not in the genotype file and will not be analyzed: ",
            paste(removed_progeny, collapse = ", "), call. = FALSE)
    progeny_candidates <- progeny_candidates[id %in% valid_ids]
  }

  if (!"sex" %in% base::colnames(all_parents)) {
    warning("No 'sex' column in parents file. All parents treated as ambiguous ('A').")
    all_parents[, sex := "A"]
  }

  all_parents[, sex := base::toupper(sex)]
  male_parent_candidates   <- all_parents[sex %in% c("M", "A", "NA"), .SD]
  female_parent_candidates <- all_parents[sex %in% c("F", "A", "NA")]

  if (base::nrow(male_parent_candidates) == 0 && method %in% c("best_male_parent", "best_pair"))
    warning("No valid male parent candidates remain after filtering.", call. = FALSE)
  if (base::nrow(female_parent_candidates) == 0 && method %in% c("best_female_parent", "best_pair"))
    warning("No valid female parent candidates remain after filtering.", call. = FALSE)
  if (base::nrow(progeny_candidates) == 0)
    stop("No valid progeny candidates remain after filtering.")

  #### Pre-compute genotype matrices once ####
  genos_mat   <- base::as.matrix(genos, rownames = "id")
  genos_hom   <- data.table::copy(genos)
  marker_cols <- base::setdiff(base::names(genos_hom), "id")
  for (col in marker_cols)
    genos_hom[base::get(col) == 1, (col) := NA_integer_]
  genos_hom_mat <- base::as.matrix(genos_hom, rownames = "id")

  #### Status helper ####
  assign_status <- function(markers, error_pct) {
    base::ifelse(markers < min_markers, "low_markers",
                 base::ifelse(error_pct > error_threshold, "high_error", "pass"))
  }

  #### Logic for Homozygous Matching Methods ####
  if (method %in% c("best_male_parent", "best_female_parent", "best_match")) {
    parent_ids <- base::switch(method,
                               "best_male_parent"   = male_parent_candidates$id,
                               "best_female_parent" = female_parent_candidates$id,
                               "best_match"         = base::union(male_parent_candidates$id,
                                                                  female_parent_candidates$id))
    parent_genos  <- genos_hom_mat[base::rownames(genos_hom_mat) %in% parent_ids,            , drop = FALSE]
    progeny_genos <- genos_hom_mat[base::rownames(genos_hom_mat) %in% progeny_candidates$id, , drop = FALSE]
    n_progeny  <- base::nrow(progeny_genos)
    results_dt <- data.table::data.table(
      id                  = base::rownames(progeny_genos),
      best_match          = NA_character_,
      mendelian_error_pct = NA_real_,
      markers_tested      = NA_integer_,
      status              = NA_character_
    )
    for (i in base::seq_len(n_progeny)) {
      progeny_vec      <- progeny_genos[i, ]
      mismatches       <- base::rowSums(parent_genos != progeny_vec, na.rm = TRUE)
      comparisons      <- base::rowSums(!base::is.na(parent_genos) & !base::is.na(progeny_vec))
      percent_mismatch <- (mismatches / comparisons) * 100
      percent_mismatch[base::is.nan(percent_mismatch)] <- NA
      best_idx <- base::which.min(percent_mismatch)
      if (base::length(best_idx) == 0) {
        data.table::set(results_dt, i, "markers_tested", 0L)
        data.table::set(results_dt, i, "status",         "low_markers")
        next
      }
      best_markers <- comparisons[best_idx]
      best_error   <- base::round(percent_mismatch[best_idx], 2)
      data.table::set(results_dt, i, "best_match",          base::rownames(parent_genos)[best_idx])
      data.table::set(results_dt, i, "mendelian_error_pct", best_error)
      data.table::set(results_dt, i, "markers_tested",      base::as.integer(best_markers))
      data.table::set(results_dt, i, "status",              assign_status(best_markers, best_error))
    }
    final_df <- results_dt
  }

  #### Logic for Best Pair Method ####
  if (method == "best_pair") {
    parent_pairs <- data.table::CJ(male_parent   = male_parent_candidates$id,
                                   female_parent = female_parent_candidates$id)
    if (!allow_selfing) {
      parent_pairs <- parent_pairs[male_parent != female_parent]
      if (verbose) base::cat("Selfing is disallowed. Pairs with identical parents are removed.\n")
    }
    if (base::nrow(parent_pairs) == 0) stop("No valid parent pairs to test.")

    male_parent_genos_mat   <- genos_mat[parent_pairs$male_parent,   , drop = FALSE]
    female_parent_genos_mat <- genos_mat[parent_pairs$female_parent, , drop = FALSE]
    progeny_ids <- progeny_candidates$id
    progeny_mat <- genos_mat[progeny_ids, , drop = FALSE]
    n_progeny   <- base::nrow(progeny_mat)
    n_pairs     <- base::nrow(parent_pairs)

    mismatch_mat <- base::matrix(
      base::vapply(base::seq_len(n_progeny), function(j) {
        progeny_vec <- progeny_mat[j, ]
        base::rowSums(
          (male_parent_genos_mat == 0 & female_parent_genos_mat == 0 & progeny_vec >  0) |
            (male_parent_genos_mat == 2 & female_parent_genos_mat == 2 & progeny_vec <  2) |
            ((male_parent_genos_mat == 0 & female_parent_genos_mat == 1) |
               (male_parent_genos_mat == 1 & female_parent_genos_mat == 0)) & (progeny_vec == 2) |
            ((male_parent_genos_mat == 2 & female_parent_genos_mat == 1) |
               (male_parent_genos_mat == 1 & female_parent_genos_mat == 2)) & (progeny_vec == 0) |
            ((male_parent_genos_mat == 0 & female_parent_genos_mat == 2) |
               (male_parent_genos_mat == 2 & female_parent_genos_mat == 0)) & (progeny_vec != 1),
          na.rm = TRUE
        )
      }, numeric(n_pairs)),
      nrow = n_pairs, ncol = n_progeny
    )

    comparison_mat <- base::matrix(
      base::vapply(base::seq_len(n_progeny), function(j) {
        progeny_vec <- progeny_mat[j, ]
        base::rowSums(!base::is.na(male_parent_genos_mat) &
                        !base::is.na(female_parent_genos_mat) &
                        !base::is.na(progeny_vec))
      }, numeric(n_pairs)),
      nrow = n_pairs, ncol = n_progeny
    )

    pct_mismatch_mat <- (mismatch_mat / comparison_mat) * 100
    pct_mismatch_mat[base::is.nan(pct_mismatch_mat)] <- NA

    results_dt <- data.table::data.table(
      id                  = progeny_ids,
      male_parent         = NA_character_,
      female_parent       = NA_character_,
      mendelian_error_pct = NA_character_,
      markers_tested      = NA_integer_,
      status              = NA_character_
    )

    results_list <- base::vector("list", n_progeny)
    for (j in base::seq_len(n_progeny)) {
      prog_id          <- progeny_ids[j]
      percent_mismatch <- pct_mismatch_mat[, j]
      comparisons      <- comparison_mat[,  j]
      min_mismatch_val <- base::min(percent_mismatch, na.rm = TRUE)

      if (base::is.infinite(min_mismatch_val)) {
        data.table::set(results_dt, j, "markers_tested", 0L)
        data.table::set(results_dt, j, "status",         "low_markers")
        next
      }

      best_indices <- base::which(percent_mismatch == min_mismatch_val)
      if (base::length(best_indices) > 1) {
        best_markers_per_pair <- comparisons[best_indices]
        max_markers           <- base::max(best_markers_per_pair)
        best_indices          <- best_indices[best_markers_per_pair == max_markers]
      }

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

      data.table::set(results_dt, j, "male_parent",         best_pairs$male_parent[1])
      data.table::set(results_dt, j, "female_parent",       best_pairs$female_parent[1])
      data.table::set(results_dt, j, "mendelian_error_pct", base::sprintf("%.2f", min_mismatch_val))
      data.table::set(results_dt, j, "markers_tested",      base::as.integer(best_markers))
      data.table::set(results_dt, j, "status",              a_status)

      if (show_ties && num_to_report > 1) {
        tie_row <- base::list(id = prog_id)
        for (k in base::seq(2, num_to_report)) {
          tie_row[[base::paste0("male_parent_",         k)]] <- best_pairs$male_parent[k]
          tie_row[[base::paste0("female_parent_",       k)]] <- best_pairs$female_parent[k]
          tie_row[[base::paste0("mendelian_error_pct_", k)]] <- min_mismatch_val
          tie_row[[base::paste0("markers_tested_",      k)]] <- comparisons[best_indices[k]]
        }
        results_list[[j]] <- data.table::as.data.table(tie_row)
      }
    }

    tie_rows <- data.table::rbindlist(
      base::Filter(Negate(base::is.null), results_list),
      fill      = TRUE,
      use.names = TRUE
    )
    if (base::nrow(tie_rows) > 0) {
      final_df <- merge(results_dt, tie_rows, by = "id", all.x = TRUE)
      for (col in base::names(final_df))
        data.table::set(final_df, which(final_df[[col]] == ""), col, NA_character_)
    } else {
      final_df <- results_dt
    }
  }

  #### Compile named list ####
  output_list <- list(
    pass         = final_df[status == "pass"],
    high_error   = final_df[status == "high_error"],
    low_markers  = final_df[status == "low_markers"],
    full_results = final_df
  )

  #### Verbose output ####
  if (verbose) {
    total_progeny <- base::nrow(final_df)
    base::cat("\n=== Parentage Assignment Report ===\n")
    base::cat("\nTotal progeny evaluated:", total_progeny, "\n")
    base::cat("Method:", method, "  |  ",
              "Error threshold:", error_threshold, "%  |  ",
              "Minimum markers:", min_markers, "\n")

    n_pass <- base::nrow(output_list$pass)
    if (n_pass > 0) {
      base::cat(base::sprintf("\n%d progeny passed (%.1f%%).\n",
                              n_pass, (n_pass / total_progeny) * 100))
    } else {
      base::cat("\nNo progeny passed.\n")
    }

    n_high <- base::nrow(output_list$high_error)
    if (n_high > 0) {
      base::cat(base::sprintf("\n%d progeny flagged high_error (%.1f%%):\n",
                              n_high, (n_high / total_progeny) * 100))
      base::print(output_list$high_error)
    } else {
      base::cat("\nNo progeny flagged for high error.\n")
    }

    n_low <- base::nrow(output_list$low_markers)
    if (n_low > 0) {
      base::cat(base::sprintf("\n%d progeny flagged low_markers (%.1f%%):\n",
                              n_low, (n_low / total_progeny) * 100))
      base::print(output_list$low_markers)
    } else {
      base::cat("\nNo progeny flagged for low marker count.\n")
    }

    base::cat("\nFull results are included in the returned list as $full_results.\n")
  }

  #### Plot Results ####
  p <- NULL
  if (plot_results) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 is required for plot_results = TRUE. Please install it.", call. = FALSE)
    } else {
      plot_df <- final_df[!is.na(final_df$mendelian_error_pct)]
      plot_df$mendelian_error_pct <- base::as.numeric(plot_df$mendelian_error_pct)
      plot_df$plot_status <- base::ifelse(
        plot_df$status == "pass",        "pass",
        base::ifelse(
          plot_df$status == "high_error",  "high_error",
          base::ifelse(
            plot_df$status == "low_markers", "low_markers",
            "other"
          )
        )
      )

      n_total <- base::nrow(plot_df)
      n_pass  <- base::sum(plot_df$status == "pass",        na.rm = TRUE)
      n_high  <- base::sum(plot_df$status == "high_error",  na.rm = TRUE)
      n_low   <- base::sum(plot_df$status == "low_markers", na.rm = TRUE)

      threshold_label <- base::paste0(
        "Error Threshold: ", error_threshold, "%  |  ",
        "Pass: ",        n_pass, "  |  ",
        "High Error: ",  n_high, "  |  ",
        "Low Markers: ", n_low
      )

      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = mendelian_error_pct, fill = plot_status)
      ) +
        ggplot2::geom_histogram(binwidth = 1, color = "white", alpha = 0.9) +
        ggplot2::geom_vline(xintercept = error_threshold,
                            linetype = "dashed", color = "black", linewidth = 1) +
        ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 5)) +
        ggplot2::scale_y_continuous(breaks = seq(0, 10000, by = 5)) +
        ggplot2::scale_fill_manual(
          values = c("pass"        = "#339900",
                     "high_error"  = "#cc3333",
                     "low_markers" = "#F1C40F",
                     "other"       = "#BDC3C7"),
          labels = c("pass"        = "Pass",
                     "high_error"  = "High Error",
                     "low_markers" = "Low Markers",
                     "other"       = "Other")
        ) +
        ggplot2::labs(
          title    = "Parentage Mendelian Error Distribution",
          subtitle = base::paste0("Progeny Tested: ", n_total, "\n \n", threshold_label),
          x        = "Mendelian Error (%)",
          y        = "Number of Progeny",
          fill     = "Status"
        ) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(legend.position = "top")

      base::print(p)
    }
  }

  output_list$plot <- p
  return(base::invisible(output_list))
}
