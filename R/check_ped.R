#' Check and optionally correct a pedigree file
#'
#' Reads a tab-separated pedigree file with columns `id`, `male_parent`, and
#' `female_parent` (any column order), performs structural and consistency checks,
#' reports findings, and optionally builds a partially corrected pedigree object.
#'
#' The input file is never modified. Findings are printed to the console when
#' `verbose = TRUE`, a report list is assigned to the global environment, and when
#' `correct = TRUE`, a corrected pedigree data frame is also assigned there.
#'
#' @section Checks performed:
#' \describe{
#'   \item{Missing required columns}{Stops immediately if `id`, `male_parent`, or
#'     `female_parent` are absent.}
#'   \item{Exact duplicate rows}{Rows identical across all three fields. Removed
#'     when `correct = TRUE`.}
#'   \item{Repeated IDs with conflicting parents}{Same `id` with differing parent
#'     assignments. Collapsed to one row when `correct = TRUE`; conflicting fields
#'     are set to `"0"`.}
#'   \item{Inconsistent sex roles}{IDs appearing as both `male_parent` and
#'     `female_parent`. Reported only; no automatic correction.}
#'   \item{Missing parents}{Parent IDs not listed as individuals. Reported, and
#'     added as founder rows (`"0"/"0"`) when `correct = TRUE`.}
#'   \item{Cycles}{Direct or indirect ancestry loops. Reported only; must be
#'     resolved manually.}
#' }
#'
#' @section Correction behavior:
#' When `correct = TRUE`, the function removes exact duplicates, collapses
#' conflicting repeated IDs (replacing inconsistent parent fields with `"0"`),
#' and appends missing parents as founders. Inconsistent sex roles and cycles are
#' not automatically corrected. Running the function again on the corrected
#' pedigree is recommended to reassess downstream issues.
#'
#' @param ped.file Path to the pedigree text file.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical. Print findings to the console? Default `TRUE`.
#' @param correct Logical. Build and save a corrected pedigree? Default `TRUE`.
#' @param save_zip Logical. Export report components as a `.zip` archive? Default `FALSE`.
#' @param save_corrected_zip Logical. Include corrected pedigree in the `.zip`?
#'   Only applies when `save_zip = TRUE` and `correct = TRUE`. Default `FALSE`.
#'
#' @return An invisible named list of data frames:
#' \describe{
#'   \item{exact_duplicates}{Exact duplicate rows.}
#'   \item{repeated_ids_diff}{Rows with the same `id` but conflicting parent assignments.}
#'   \item{inconsistent_sex_roles}{Rows involving IDs used as both male and female parent.}
#'   \item{missing_parents}{Rows referencing parent IDs absent from `id`.}
#'   \item{dependencies}{Rows whose `id` is involved in a detected cycle.}
#' }
#'
#' @author Josué Chinchilla-Vargas
#'
#' @importFrom dplyr %>% group_by filter ungroup distinct mutate summarize first bind_rows n_distinct n select row_number
#' @importFrom stats setNames
#' @importFrom utils read.table write.table zip
#' @importFrom tools file_path_sans_ext
#' @export
check_ped <- function(ped.file,
                      seed               = NULL,
                      verbose            = TRUE,
                      correct            = TRUE,  #add arguments for removing/filtering conflicting ids, repeated ids and founders are always fixed, dependencies are always manual.
                      save_zip           = FALSE, #remove all save zips, export each report and corrected_pedigree as a list
                      save_corrected_zip = FALSE) {

  #### setup ####
  if (!is.null(seed)) set.seed(seed)
  data <- utils::read.table(ped.file, header = TRUE)

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
    mutate(
      id            = as.character(id),
      male_parent   = as.character(male_parent),
      female_parent = as.character(female_parent)
    )
  # Add row numbers before any processing
  data <- data %>% mutate(row_number = row_number(), .before = id)
  errors <- list()
  # Add row numbers before any processing
  data <- data %>% mutate(row_number = row_number(), .before = id)
  errors <- list()

  #### check 1: exact duplicates ####
  exact_duplicates <- data[
    duplicated(data %>% select(-row_number)) |
      duplicated(data %>% select(-row_number), fromLast = TRUE),
  ]

  #### check 2: repeated IDs with conflicting male_parent/female_parent ####
  repeated_ids_diff <- data %>%
    group_by(id) %>%
    filter(n() > 1) %>%
    filter(n_distinct(male_parent) > 1 | n_distinct(female_parent) > 1) %>%
    ungroup()

  #### check 3: inconsistent parent sex roles ####
  male_ids   <- unique(data$male_parent[data$male_parent   != "0"])
  female_ids <- unique(data$female_parent[data$female_parent != "0"])
  messy_ids  <- intersect(male_ids, female_ids)
  inconsistent_sex_roles <- data %>%
    filter(male_parent %in% messy_ids | female_parent %in% messy_ids)

  #### check 4: missing parents ####
  all_ids     <- unique(data$id)
  ref_ids     <- unique(c(data$male_parent, data$female_parent))
  ref_ids     <- ref_ids[ref_ids != "0"]
  missing_ids <- setdiff(ref_ids, all_ids)
  missing_parents <- data %>%
    filter(male_parent %in% missing_ids | female_parent %in% missing_ids)

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
  dependencies <- data %>% filter(id %in% cycle_ids)
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
    corrected <- data %>% select(-row_number)
    # Remove exact duplicates
    corrected <- distinct(corrected)
    # Resolve conflicting IDs
    if (nrow(repeated_ids_diff) > 0) {
      corrected <- corrected %>%
        group_by(id) %>%
        summarize(
          male_parent   = if (n_distinct(male_parent)   > 1) "0" else first(male_parent),
          female_parent = if (n_distinct(female_parent) > 1) "0" else first(female_parent),
          .groups = "drop"
        )
    }
    # Add missing founders
    if (length(missing_ids) > 0) {
      corrected <- bind_rows(
        corrected,
        missing_founders %>% select(-row_number)
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
      cat(paste0("\n`correct = TRUE`: saving corrected pedigree as `", corrected_name, "`.\n"))
    } else {
      cat("\n`correct = FALSE`: no changes made to the pedigree.\n")
    }
  }

  # Always save report to global environment
  assign(report_name, input_ped_report, envir = .GlobalEnv)
  # Save corrected pedigree if correct = TRUE
  if (correct) assign(corrected_name, corrected, envir = .GlobalEnv)

  #### zip export ####
  if (save_zip) {
    tmp_dir   <- tempfile()
    dir.create(tmp_dir)
    zip_files <- character(0)

    # Section headers matching the console output labels
    section_labels <- c(
      exact_duplicates       = "Exact Duplicates",
      repeated_ids_diff      = "Conflicting IDs",
      inconsistent_sex_roles = "Inconsistent Parent Sex Roles",
      missing_parents        = "Missing Parents (rows that reference them)",
      dependencies           = "Cycles / Dependencies"
    )

    for (component in names(section_labels)) {
      df        <- input_ped_report[[component]]
      file_path <- file.path(tmp_dir, paste0(file_base, "_", component, ".txt"))

      if (nrow(df) == 0) {
        # Write header + "None found." for clean empty reports
        writeLines(
          c(paste0("--- ", section_labels[[component]], " ---"), "None found."),
          con = file_path
        )
      } else {
        # Write header, then column names, then data — avoids col.names warning
        file_con <- file(file_path, open = "wt")
        writeLines(paste0("--- ", section_labels[[component]], " ---"), con = file_con)
        writeLines(paste(colnames(df), collapse = "\t"), con = file_con)
        close(file_con)
        utils::write.table(df, file = file_path, sep = "\t", row.names = FALSE,
                           quote = FALSE, append = TRUE, col.names = FALSE)
      }

      zip_files <- c(zip_files, file_path)
    }

    # Optionally include corrected pedigree in zip
    if (save_corrected_zip && correct) {
      corrected_path <- file.path(tmp_dir, paste0(corrected_name, ".txt"))
      utils::write.table(corrected, file = corrected_path, sep = "\t",
                         row.names = FALSE, quote = FALSE)
      zip_files <- c(zip_files, corrected_path)
    }

    # Bundle all files into zip in the working directory
    zip_path <- file.path(getwd(), zip_name)
    utils::zip(zipfile = zip_path, files = zip_files, flags = "-j")
    if (verbose) cat(paste0("\nZip archive saved to: ", zip_path, "\n"))
  }

  invisible(input_ped_report)
}
