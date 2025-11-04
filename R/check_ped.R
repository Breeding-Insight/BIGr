#' Check a pedigree file for accuracy and report/correct common errors
#'
#' `check_ped` reads a 3-column pedigree file (tab-separated, columns labeled `id`, `sire`, `dam` in any order)
#' and performs quality checks, optionally correcting or flagging errors.
#'
#' The function checks for:
#' * Exact duplicate rows and removes them (keeping one copy)
#' * IDs that appear more than once with conflicting sire/dam assignments (sets sire/dam to "0")
#' * IDs that appear in both sire and dam columns
#' * Missing parents (IDs referenced as sire/dam but not in `id` column), adds them with sire/dam = "0"
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
#' @param verbose Logical. If TRUE (default), prints errors and prompts for interactive saving.
#'
#' @return A list of data.frames containing detected issues:
#' * `exact_duplicates`: rows that were exact duplicates
#' * `repeated_ids_diff`: IDs appearing more than once with conflicting sire/dam
#' * `messy_parents`: IDs appearing as both sire and dam
#' * `missing_parents`: parents added to the pedigree with 0 as sire/dam
#' * `dependencies`: detected cycles in the pedigree
#'
#' @examples
#' ped_file <- system.file("check_ped_test.txt", package = "BIGr")
#' ped_errors <- check_ped(ped.file = ped_file, seed = 101919)
#'
#' # Access messy parents
#' ped_errors$messy_parents
#'
#' # IDs with messy parents
#' messy_ids <- ped_errors$messy_parents$id
#' print(messy_ids)
#'
#' @import dplyr
#' @import janitor
#' @importFrom stats setNames
#' @importFrom utils read.table
#' @export
check_ped <- function(ped.file, seed = NULL, verbose = TRUE) {

  #### setup ####
  if (!is.null(seed)) set.seed(seed)

  # Read and clean data
  data <- utils::read.table(ped.file, header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(
      id = as.character(id),
      sire = as.character(sire),
      dam = as.character(dam)
    )

  original_data <- data
  errors <- list()
  missing_parents <- data.frame(id = character(), sire = character(), dam = character(), stringsAsFactors = FALSE)

  #### check 1: exact duplicates ####
  exact_duplicates <- data[duplicated(data), ]
  if (nrow(exact_duplicates) > 0) {
    data <- distinct(data)  # remove exact duplicates
  }

  #### check 2: repeated IDs with conflicting sire/dam ####
  repeated_ids <- data %>%
    group_by(id) %>%
    filter(n() > 1) %>%
    ungroup()

  # Only IDs with actual conflicting sire/dam
  conflicting_ids <- repeated_ids %>%
    group_by(id) %>%
    filter(n_distinct(sire) > 1 | n_distinct(dam) > 1) %>%
    ungroup()

  if (nrow(conflicting_ids) > 0) {
    # Keep one row per ID, set sire/dam to "0"
    data <- data %>%
      group_by(id) %>%
      summarize(
        sire = if(n_distinct(sire) > 1) "0" else first(sire),
        dam = if(n_distinct(dam) > 1) "0" else first(dam),
        .groups = "drop"
      )
  }

  repeated_ids_report <- conflicting_ids

  #### check 3: missing parents ####
  for (i in 1:nrow(data)) {
    id <- data$id[i]
    sire <- data$sire[i]
    dam <- data$dam[i]

    if (sire != "0" && sire != id && !sire %in% data$id) {
      missing_parents <- rbind(missing_parents, data.frame(id = sire, sire = "0", dam = "0", stringsAsFactors = FALSE))
    }
    if (dam != "0" && dam != id && !dam %in% data$id) {
      missing_parents <- rbind(missing_parents, data.frame(id = dam, sire = "0", dam = "0", stringsAsFactors = FALSE))
    }

    if (sire == id || dam == id) {
      errors <- append(errors, paste("Dependency: Individual", id, "cannot be its own parent"))
    }
  }

  missing_parents <- distinct(missing_parents)
  if (nrow(missing_parents) > 0) {
    data <- bind_rows(data, missing_parents)
  }

  #### check 4: messy parents ####
  sire_ids <- unique(data$sire[data$sire != "0"])
  dam_ids <- unique(data$dam[data$dam != "0"])
  messy_ids <- intersect(sire_ids, dam_ids)
  messy_parents <- data %>% filter(id %in% messy_ids)

  #### check 5: dependencies (cycles) ####
  detect_all_cycles <- function(data) {
    adj_list <- lapply(data$id, function(x) {
      row <- data[data$id == x, ]
      c(row$sire, row$dam)
    })
    names(adj_list) <- data$id

    dfs <- function(node, visited, rec_stack, path) {
      visited[node] <- TRUE
      rec_stack[node] <- TRUE
      path <- append(path, node)
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

    visited <- stats::setNames(rep(FALSE, length(adj_list)), names(adj_list))
    rec_stack <- stats::setNames(rep(FALSE, length(adj_list)), names(adj_list))
    all_cycles <- list()

    for (node in names(adj_list)) {
      if (!visited[node]) {
        node_cycles <- dfs(node, visited, rec_stack, character())
        if (length(node_cycles) > 0) {
          all_cycles <- append(all_cycles, node_cycles)
        }
      }
    }
    return(all_cycles)
  }

  cycles <- detect_all_cycles(data)
  if (length(cycles) > 0) {
    for (cycle_group in cycles) {
      cycle_ids <- unique(unlist(cycle_group))
      errors <- append(errors, paste("Cycle detected involving IDs:", paste(cycle_ids, collapse = " -> ")))
    }
  }

  #### compile findings ####
  input_ped_report <- list(
    exact_duplicates = exact_duplicates,
    repeated_ids_diff = repeated_ids_report,
    messy_parents = messy_parents,
    missing_parents = missing_parents,
    dependencies = data.frame(Dependency = unique(unlist(errors)))
  )

  #### file names ####
  file_base <- tools::file_path_sans_ext(basename(ped.file))
  corrected_name <- paste0(file_base, "_corrected")
  report_name <- paste0(file_base, "_report")

  #### output ####
  if (verbose) {
    cat("\n=== Pedigree Quality Check Report ===\n")

    if (nrow(exact_duplicates) > 0) {
      cat("\n Exact duplicate trios detected (only one copy will be kept in corrected pedigree):\n")
      print(exact_duplicates)
    } else cat("\nNo exact duplicate trios found.\n")

    if (nrow(repeated_ids_report) > 0) {
      cat("\nConflicting trios detected (sire/dam set to 0 in corrected pedigree):\n")
      print(repeated_ids_report)
    } else cat("\nNo conflicting repeated trios found.\n")

    if (nrow(missing_parents) > 0) {
      cat("\n Parents missing as IDs found in the pedigree (will be added as founders in corrected pedigree):\n")
      print(missing_parents)
    } else cat("\nNo missing parents found.\n")

    if (nrow(messy_parents) > 0) {
      cat("\n IDs found as both sire and dam (is selfing or hermaphrodytism possible?):\n")
      print(messy_parents)
    } else cat("\nNo IDs found as both sire and dam.\n")


    if (nrow(input_ped_report$dependencies) > 0) {
      cat("\nDependencies detected:\n")
      print(input_ped_report$dependencies)
    } else cat("\nNo dependencies detected.\n")

    #### interactive save ####
    cat(paste0("\nDo you want to save the corrected pedigree as dataframe `", corrected_name, "`? (y/n): "))
    ans <- tolower(trimws(readline()))
    if (ans == "y") {
      assign(corrected_name, data, envir = .GlobalEnv)
      assign("input_ped_report", input_ped_report, envir = .GlobalEnv)
      cat(paste0("Saved corrected pedigree as `", corrected_name, "` and report as `input_ped_report`.\n"))
    } else {
      cat("No corrected pedigree was saved.\n")
    }

  } else {
    # Silent automatic mode
    assign(corrected_name, data, envir = .GlobalEnv)
    assign(report_name, input_ped_report, envir = .GlobalEnv)
  }

  invisible(input_ped_report)
}
