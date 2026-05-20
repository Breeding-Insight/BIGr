#' Check a pedigree file for accuracy and report and correct common errors
#'
#' Reads a 3-column pedigree file (tab-separated, columns id, male_parent,
#' female_parent in any order) and performs quality checks, optionally
#' correcting detected errors.
#'
#' Exact duplicate rows and missing parents are always corrected.
#' Conflicting trios and inconsistent sex roles are corrected when their
#' respective arguments are TRUE. Cycles and dependencies are always
#' reported only and must be resolved manually.
#'
#' @param ped.file Path to the pedigree text file.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical. If TRUE (default), prints the report to the console.
#' @param correct_conflicting_trios Logical. If TRUE (default), sets conflicting
#'   male_parent and female_parent to 0 and collapses to one row per ID.
#' @param correct_inconsistent_sex_roles Logical. If TRUE (default), sets
#'   male_parent and female_parent to 0 for rows involving IDs found as both,
#'   then removes any resulting exact duplicates.
#'
#' @return An invisible named list of data frames:
#' \describe{
#'   \item{exact_duplicates}{Exact duplicate rows found in the input.}
#'   \item{conflicting_trios}{IDs with conflicting male_parent or female_parent assignments.}
#'   \item{inconsistent_sex_roles}{IDs appearing as both male_parent and female_parent.}
#'   \item{missing_parents}{Parent IDs absent from id, added as founders.}
#'   \item{dependencies}{Cycles detected in the pedigree. Must be resolved manually.}
#'   \item{corrected_pedigree}{Corrected pedigree table.}
#' }
#'
#' @examples
#' ped_file <- system.file("check_ped_test.txt", package = "BIGr")
#' ped_errors <- check_ped(ped.file = ped_file, seed = 101919, verbose = FALSE)
#'
#' ped_errors$inconsistent_sex_roles
#' ped_errors$corrected_pedigree
#'
#' conflicting_sex_ids <- ped_errors$inconsistent_sex_roles$id
#' print(conflicting_sex_ids)
#'
#' @importFrom dplyr %>% mutate filter group_by ungroup summarize distinct bind_rows select first n n_distinct if_else row_number
#' @importFrom stats setNames
#' @importFrom utils read.table
#' @export
check_ped <- function(ped.file,
                      seed                           = NULL,
                      verbose                        = TRUE,
                      correct_conflicting_trios      = TRUE,
                      correct_inconsistent_sex_roles = TRUE) {

  #### setup ####
  if (!is.null(seed)) set.seed(seed)

  data <- utils::read.table(ped.file, header = TRUE) %>%
    dplyr::mutate(
      id            = as.character(id),
      male_parent   = as.character(male_parent),
      female_parent = as.character(female_parent)
    )

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

  # Add row numbers before any processing so all reports reference original rows [1]
  data <- data %>% dplyr::mutate(row_number = dplyr::row_number(), .before = id)

  errors          <- list()
  missing_parents <- data.frame(
    row_number    = integer(),
    id            = character(),
    male_parent   = character(),
    female_parent = character(),
    stringsAsFactors = FALSE
  )

  #### check 1: exact duplicates (always fixed) ####
  exact_duplicates <- data[
    duplicated(data %>% dplyr::select(-row_number)) |
      duplicated(data %>% dplyr::select(-row_number), fromLast = TRUE),
  ]
  if (nrow(exact_duplicates) > 0) {
    data <- data %>%
      dplyr::select(-row_number) %>%
      dplyr::distinct() %>%
      dplyr::mutate(row_number = dplyr::row_number(), .before = id)
  }

  #### check 2: repeated IDs with conflicting male_parent/female_parent ####
  repeated_ids <- data %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()

  conflicting_ids <- repeated_ids %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::n_distinct(male_parent) > 1 | dplyr::n_distinct(female_parent) > 1) %>%
    dplyr::ungroup()

  if (correct_conflicting_trios && nrow(conflicting_ids) > 0) {
    # Set conflicting parents to "0" -- rows become exact duplicates, summarize collapses to one [1]
    data <- data %>%
      dplyr::group_by(id) %>%
      dplyr::summarize(
        row_number    = dplyr::first(row_number),
        male_parent   = if (dplyr::n_distinct(male_parent)   > 1) "0" else dplyr::first(male_parent),
        female_parent = if (dplyr::n_distinct(female_parent) > 1) "0" else dplyr::first(female_parent),
        .groups = "drop"
      ) %>%
      dplyr::select(row_number, id, male_parent, female_parent)
  }

  repeated_ids_report <- conflicting_ids

  #### check 3: missing parents (always fixed) ####
  for (i in seq_len(nrow(data))) {
    id            <- data$id[i]
    male_parent   <- data$male_parent[i]
    female_parent <- data$female_parent[i]

    if (male_parent != "0" && male_parent != id && !male_parent %in% data$id) {
      missing_parents <- rbind(
        missing_parents,
        data.frame(row_number = data$row_number[i], id = male_parent,
                   male_parent = "0", female_parent = "0",
                   stringsAsFactors = FALSE)
      )
    }
    if (female_parent != "0" && female_parent != id && !female_parent %in% data$id) {
      missing_parents <- rbind(
        missing_parents,
        data.frame(row_number = data$row_number[i], id = female_parent,
                   male_parent = "0", female_parent = "0",
                   stringsAsFactors = FALSE)
      )
    }

    if (male_parent == id || female_parent == id) {
      errors <- append(errors, paste("Dependency: Individual", id, "cannot be its own parent"))
    }
  }

  missing_parents <- dplyr::distinct(missing_parents)
  if (nrow(missing_parents) > 0) {
    data <- dplyr::bind_rows(data, missing_parents)
  }

  #### check 4: inconsistent sex roles ####
  male_ids            <- unique(data$male_parent[data$male_parent   != "0"])
  female_ids          <- unique(data$female_parent[data$female_parent != "0"])
  conflicting_sex_ids <- intersect(male_ids, female_ids)
  inconsistent_sex_roles <- data %>% dplyr::filter(id %in% conflicting_sex_ids)

  if (correct_inconsistent_sex_roles && length(conflicting_sex_ids) > 0) {
    # Zero out male_parent/female_parent wherever a conflicting ID appears,
    # then distinct() removes any rows that became exact duplicates [1]
    data <- data %>%
      dplyr::mutate(
        male_parent   = dplyr::if_else(male_parent   %in% conflicting_sex_ids, "0", male_parent),
        female_parent = dplyr::if_else(female_parent %in% conflicting_sex_ids, "0", female_parent)
      ) %>%
      dplyr::distinct(id, male_parent, female_parent, .keep_all = TRUE)
  }

  #### check 5: dependencies (cycles) -- reported only, never modified ####
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
      errors    <- append(errors,
                          paste("Cycle detected involving IDs:",
                                paste(cycle_ids, collapse = " -> ")))
    }
  }

  #### compile findings ####
  input_ped_report <- list(
    exact_duplicates       = exact_duplicates,
    conflicting_trios      = repeated_ids_report,
    inconsistent_sex_roles = inconsistent_sex_roles,
    missing_parents        = missing_parents,
    dependencies           = data.frame(Dependency = unique(unlist(errors)),
                                        stringsAsFactors = FALSE),
    corrected_pedigree     = data %>% dplyr::select(-row_number)
  )

  #### output ####
  if (verbose) {
    cat("\n=== Pedigree Quality Check Report ===\n")

    if (nrow(exact_duplicates) > 0) {
      cat("\nExact duplicate trios detected (removed in corrected pedigree):\n")
      print(exact_duplicates)
    } else cat("\nNo exact duplicate trios found.\n")

    if (nrow(repeated_ids_report) > 0) {
      cat("\nConflicting trios detected:\n")
      print(repeated_ids_report)
      if (correct_conflicting_trios) {
        cat("  -> parents set to 0 and collapsed to one row in corrected pedigree.\n")
      } else {
        cat("  -> correct_conflicting_trios = FALSE: left as-is in corrected pedigree.\n")
      }
    } else cat("\nNo conflicting trios found.\n")

    if (nrow(missing_parents) > 0) {
      cat("\nParents missing as IDs found in the pedigree (added as founders in corrected pedigree):\n")
      print(missing_parents)
    } else cat("\nNo missing parents found.\n")

    if (nrow(inconsistent_sex_roles) > 0) {
      cat("\nIDs found as both male_parent and female_parent (is selfing or hermaphrodytism possible?):\n")
      print(inconsistent_sex_roles)
      if (correct_inconsistent_sex_roles) {
        cat("  -> parent fields set to 0 for conflicting IDs; exact duplicates removed in corrected pedigree.\n")
      } else {
        cat("  -> correct_inconsistent_sex_roles = FALSE: left as-is in corrected pedigree.\n")
      }
    } else cat("\nNo IDs found as both male_parent and female_parent.\n")

    if (nrow(input_ped_report$dependencies) > 0) {
      cat("\nDependencies detected (must be resolved manually):\n")
      print(input_ped_report$dependencies)
    } else cat("\nNo dependencies detected.\n")

    cat("\nThe corrected pedigree is included in the returned list as corrected_pedigree.\n")
  }

  invisible(input_ped_report)
}
