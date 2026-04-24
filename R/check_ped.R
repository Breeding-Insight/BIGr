#' Check a pedigree file for accuracy and report/correct common errors
#'
#' `check_ped` reads a 3-column pedigree file (tab-separated, columns labeled `id`, `male_parent`, `female_parent` in any order)
#' and performs quality checks, optionally correcting or flagging errors.
#'
#' The function checks for:
#' * Exact duplicate rows and removes them (keeping one copy)
#' * IDs that appear more than once with conflicting male_parent/female_parent assignments (sets both to "0")
#' * IDs that appear in both male_parent and female_parent columns
#' * Missing parents (IDs referenced as male_parent/female_parent but not in `id` column), adds them with both = "0"
#' * Direct and indirect pedigree dependencies (cycles), such as a parent being its own descendant
#'
#' After an initial run to clean exact duplicates and repeated IDs, you can run the function again to detect cycles more accurately.
#'
#' The function does **not** overwrite the input file. Instead, it prints findings to the console and optionally saves:
#' * Corrected pedigree as a dataframe in the global environment
#' * A report listing all detected issues
#'
#' @param ped.file Path to the pedigree text file.
#' @param seed Optional seed for reproducibility.
#' @param verbose Logical. If TRUE (default), prints findings and prompts for saving.
#' @param correct Logical. If TRUE (default), saves a corrected pedigree. If FALSE, only reports issues.
#' @param save_zip Logical. If TRUE, writes each report component as a .txt file and bundles them into a .zip archive in the working directory.
#' @param save_corrected_zip Logical. If TRUE (and save_zip = TRUE and correct = TRUE), also writes the corrected pedigree into the .zip archive.
#'
#' @return A named list of data.frames with offending rows (including row numbers):
#' * `exact_duplicates`
#' * `repeated_ids_diff`
#' * `inconsistent_sex_roles`
#' * `missing_parents`
#' * `dependencies`
#'
#' @importFrom dplyr %>% group_by filter ungroup distinct mutate summarize first bind_rows n_distinct n select
#' @importFrom janitor clean_names
#' @importFrom stats setNames
#' @importFrom utils read.table write.table zip
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
  data <- utils::read.table(ped.file, header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(
      id            = as.character(id),
      male_parent   = as.character(male_parent),
      female_parent = as.character(female_parent)
    )
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
  # For the corrected pedigree, build founder rows for truly absent parents
  missing_founders <- data.frame(
    row_number    = NA_integer_,
    id            = missing_ids,
    male_parent   = "0",
    female_parent = "0",
    stringsAsFactors = FALSE
  )

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
      cat(paste0("\nSave corrected pedigree as `", corrected_name, "`? (y/n): "))
      ans <- tolower(trimws(readline()))
      if (ans == "y") {
        assign(corrected_name, corrected, envir = .GlobalEnv)
        assign(report_name, input_ped_report, envir = .GlobalEnv)
        cat(paste0("Saved `", corrected_name, "` and `", report_name, "`.\n"))
      } else {
        cat("No corrected pedigree saved.\n")
        assign(report_name, input_ped_report, envir = .GlobalEnv)
      }
    } else {
      cat("\n`correct = FALSE`: no changes made to the pedigree.\n")
      assign(report_name, input_ped_report, envir = .GlobalEnv)
    }
  } else {
    # Silent mode: always save report; save corrected only if correct = TRUE
    if (correct) assign(corrected_name, corrected, envir = .GlobalEnv)
    assign(report_name, input_ped_report, envir = .GlobalEnv)
  }

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

