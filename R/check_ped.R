#' @title Check a pedigree file for structural issues
#'
#' @description
#' Reads a tab‑separated pedigree file with columns `id`,
#' `male_parent`, and `female_parent` and checks for five classes of
#' errors: missing required columns, exact duplicate rows,
#' repeated IDs with conflicting parents, inconsistent parent
#' sex roles, and pedigree cycles.  The input file is never
#' modified.  A report list is always assigned to the global
#' environment; when `correct = TRUE` a corrected pedigree is
#' assigned as well.
#'
#' @param ped.file Path to the pedigree text file.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical. If `TRUE` (default), prints findings to the
#'   console.
#' @param correct Logical. If `TRUE` (default), builds and assigns a
#'   corrected pedigree to the global environment. If `FALSE`, issues
#'   are only reported.
#' @param save_zip Logical. If `TRUE`, bundles report components into a
#'   `.zip` archive in the working directory.
#' @param save_corrected_zip Logical. If `TRUE` (requires `save_zip =
#'   TRUE` and `correct = TRUE`), includes the corrected pedigree in
#'   the `.zip` archive.
#'
#' @return An invisible named list of data frames describing detected
#'   issues:
#'   \item{exact_duplicates}{exact duplicate rows}
#'   \item{repeated_ids_diff}{rows where the same `id` has conflicting parents}
#'   \item{inconsistent_sex_roles}{rows involving IDs used as both parents}
#'   \item{missing_parents}{rows referencing parent IDs absent from `id`}
#'   \item{dependencies}{rows whose `id` is involved in a pedigree cycle}
#'
#' @importFrom dplyr %>% group_by filter ungroup distinct mutate
#'   summarize first bind_rows n_distinct n select row_number
#' @importFrom stats setNames
#' @importFrom utils read.table write.table zip capture.output
#' @importFrom tools file_path_sans_ext
#' @export
check_ped <- function(ped.file,
                      seed               = NULL,
                      verbose            = TRUE,
                      correct            = TRUE,
                      save_zip           = FALSE,
                      save_corrected_zip = FALSE) {
  #### setup ####
  if (!is.null(seed)) set.seed(seed)
  data <- utils::read.table(ped.file, header = TRUE)

  # ── Column alias remapping ──────────────────────────────────────────────────
  col_aliases <- list(
    id            = c("id", "ID", "animal", "Animal", "ind", "Ind"),
    male_parent   = c("male_parent", "sire",   "Sire",   "father", "Father", "pat", "Pat"),
    female_parent = c("female_parent", "dam",  "Dam",    "mother", "Mother", "mat", "Mat")
  )
  for (canonical in names(col_aliases)) {
    if (!canonical %in% colnames(data)) {
      match <- intersect(col_aliases[[canonical]], colnames(data))
      if (length(match) > 0) {
        colnames(data)[colnames(data) == match[1]] <- canonical
      }
    }
  }

  # Validate required column names before any processing
  required_cols <- c("id", "male_parent", "female_parent")
  missing_cols  <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      "Input file is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      ".\nExpected columns: id, male_parent, female_parent."
    )
  }
  data <- data %>%
    dplyr::mutate(
      id            = as.character(id),
      male_parent   = as.character(male_parent),
      female_parent = as.character(female_parent)
    )

  # Add row numbers before any processing
  data   <- data %>% dplyr::mutate(row_number = dplyr::row_number(), .before = id)
  errors <- list()

  #### check 1: exact duplicates ####
  exact_duplicates <- data[
    duplicated(data %>% dplyr::select(-row_number)) |
      duplicated(data %>% dplyr::select(-row_number), fromLast = TRUE),
  ]

  #### check 2: repeated IDs with conflicting male_parent/female_parent ####
  repeated_ids_diff <- data %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::filter(dplyr::n_distinct(male_parent) > 1 |
                    dplyr::n_distinct(female_parent) > 1) %>%
    dplyr::ungroup()

  #### check 3: inconsistent parent sex roles ####
  male_ids   <- unique(data$male_parent[data$male_parent   != "0"])
  female_ids <- unique(data$female_parent[data$female_parent != "0"])
  messy_ids  <- intersect(male_ids, female_ids)
  inconsistent_sex_roles <- data %>%
    dplyr::filter(male_parent %in% messy_ids | female_parent %in% messy_ids)

  #### check 4: missing parents ####
  all_ids <- unique(data$id)
  ref_ids <- unique(c(data$male_parent, data$female_parent))

  # Trim, remove NAs and zeros
  ref_ids <- trimws(ref_ids)
  ref_ids <- ref_ids[!is.na(ref_ids) & ref_ids != "" & ref_ids != "0"]
  missing_ids     <- setdiff(ref_ids, all_ids)
  missing_parents <- data %>%
    dplyr::filter(male_parent %in% missing_ids | female_parent %in% missing_ids)

  # Only build founder rows if there are actually missing parents to add
  if (length(missing_ids) > 0) {
    missing_founders <- data.frame(
      row_number    = NA_integer_,
      id            = missing_ids,
      male_parent   = "0",
      female_parent = "0",
      stringsAsFactors = FALSE
    )
  } else {
    missing_founders <- data.frame(
      row_number    = integer(0),
      id            = character(0),
      male_parent   = character(0),
      female_parent = character(0),
      stringsAsFactors = FALSE
    )
  }

  #### check 5: dependencies (cycles) ####
  detect_all_cycles <- function(data) {
    adj_list <- lapply(data$id, function(x) {
      row <- data[data$id == x, ]
      c(row$male_parent, row$female_parent)
    })
    names(adj_list) <- data$id

    dfs <- function(node, visited, rec_stack, path) {
      visited[node]   <- TRUE
      rec_stack[node] <- TRUE
      path   <- append(path, node)
      cycles <- list()
      for (neighbor in adj_list[[node]]) {
        if (neighbor %in% names(adj_list)) {
          if (!visited[neighbor]) {
            cycles <- append(cycles, dfs(neighbor, visited, rec_stack, path))
          } else if (rec_stack[neighbor]) {
            cycle_start <- match(neighbor, path)
            cycles <- append(cycles, list(path[cycle_start:length(path)]))
          }
        }
      }
      rec_stack[node] <- FALSE
      return(cycles)
    }

    visited   <- stats::setNames(rep(FALSE, length(adj_list)), names(adj_list))
    rec_stack <- stats::setNames(rep(FALSE, length(adj_list)), names(adj_list))
    all_cycles <- list()
    for (node in names(adj_list)) {
      if (!visited[node]) {
        node_cycles <- dfs(node, visited, rec_stack, character())
        if (length(node_cycles) > 0) all_cycles <- append(all_cycles, node_cycles)
      }
    }
    return(all_cycles)
  }
  cycles    <- detect_all_cycles(data)
  cycle_ids <- character(0)
  if (length(cycles) > 0) {
    for (cycle_group in cycles) {
      ids       <- unique(unlist(cycle_group))
      cycle_ids <- unique(c(cycle_ids, ids))
      errors    <- append(errors,
                          paste("Cycle detected involving IDs:",
                                paste(ids, collapse = " -> ")))
    }
  }
  dependencies <- data %>% dplyr::filter(id %in% cycle_ids)
  if (nrow(dependencies) == 0 && length(errors) > 0) {
    dependencies <- data.frame(Dependency = unique(unlist(errors)),
                               stringsAsFactors = FALSE)
  }

  #### compile report ####
  input_ped_report <- list(
    exact_duplicates       = exact_duplicates,
    repeated_ids_diff      = repeated_ids_diff,
    inconsistent_sex_roles = inconsistent_sex_roles,
    missing_parents        = missing_parents,
    dependencies           = dependencies
  )

  #### build corrected pedigree (no row_number column) ####
  if (correct) {
    corrected <- data %>% dplyr::select(-row_number)

    # Remove exact duplicates
    corrected <- dplyr::distinct(corrected)

    # Resolve conflicting IDs
    if (nrow(repeated_ids_diff) > 0) {
      corrected <- corrected %>%
        dplyr::group_by(id) %>%
        dplyr::summarize(
          male_parent   = if (dplyr::n_distinct(male_parent)   > 1) "0" else dplyr::first(male_parent),
          female_parent = if (dplyr::n_distinct(female_parent) > 1) "0" else dplyr::first(female_parent),
          .groups = "drop"
        )
    }

    # Add missing founders
    if (length(missing_ids) > 0) {
      corrected <- dplyr::bind_rows(
        corrected,
        missing_founders %>% dplyr::select(-row_number)
      )
    }
  }

  #### file names ####
  file_base      <- tools::file_path_sans_ext(basename(ped.file))
  corrected_name <- paste0(file_base, "_corrected")
  report_name    <- paste0(file_base, "_report")
  zip_name       <- paste0(file_base, "_report.zip")

  #### output ####
  if (verbose) {
    cat("\n=== Pedigree Quality Check Report ===\n")
    cat("\n--- Exact Duplicates ---\n")
    if (nrow(exact_duplicates) > 0) print(exact_duplicates) else cat("None found.\n")
    cat("\n--- Conflicting IDs ---\n")
    if (nrow(repeated_ids_diff) > 0) print(repeated_ids_diff) else cat("None found.\n")
    cat("\n--- Inconsistent Parent Sex Roles ---\n")
    if (nrow(inconsistent_sex_roles) > 0) print(inconsistent_sex_roles) else cat("None found.\n")
    cat("\n--- Missing Parents (rows that reference them) ---\n")
    if (nrow(missing_parents) > 0) print(missing_parents) else cat("None found.\n")
    cat("\n--- Cycles / Dependencies ---\n")
    if (nrow(dependencies) > 0) print(dependencies) else cat("None found.\n")
    if (correct) {
      cat(paste0("`correct = TRUE`: saving corrected pedigree as `", corrected_name, "`.\n"))
    } else {
      cat("`correct = FALSE`: no changes made to the pedigree.\n")
    }
  }

  # Section headers matching the console output labels
  section_labels <- c(
    exact_duplicates       = "Exact Duplicates",
    repeated_ids_diff      = "Conflicting IDs",
    inconsistent_sex_roles = "Inconsistent Parent Sex Roles",
    missing_parents        = "Missing Parents (rows that reference them)",
    dependencies           = "Cycles / Dependencies"
  )

  if (save_zip) {
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    zip_files <- character(0)

    for (component in names(section_labels)) {
      df        <- input_ped_report[[component]]
      file_path <- file.path(tmp_dir, paste0(file_base, "_", component, ".txt"))
      if (nrow(df) == 0) {
        writeLines(
          c(paste0("--- ", section_labels[[component]], " ---"),
            "None found."),
          con = file_path
        )
      } else {
        file_con <- file(file_path, open = "wt")
        writeLines(paste0("--- ", section_labels[[component]], " ---"),
                   con = file_con)
        writeLines(paste(colnames(df), collapse = "\t"),
                   con = file_con)
        close(file_con)
        utils::write.table(df, file = file_path, sep = "\t",
                           row.names = FALSE,
                           quote = FALSE,
                           append = TRUE,
                           col.names = FALSE)
      }
      zip_files <- c(zip_files, file_path)
    }

    if (save_corrected_zip && correct) {
      corrected_path <- file.path(tmp_dir, paste0(corrected_name, ".txt"))
      utils::write.table(corrected, file = corrected_path, sep = "\t",
                         row.names = FALSE,
                         quote = FALSE)
      zip_files <- c(zip_files, corrected_path)
    }

    zip_path <- file.path(getwd(), zip_name)
    invisible(capture.output(
      utils::zip(zipfile = zip_path, files = zip_files, flags = "-j")
    ))
    if (verbose) cat(paste0("\nZip archive saved to: ", zip_path, "\n"))
  }

  # Assign objects to the global environment
  assign(report_name, input_ped_report, envir = .GlobalEnv)
  if (correct) assign(corrected_name, corrected, envir = .GlobalEnv)

  invisible(input_ped_report)
}
