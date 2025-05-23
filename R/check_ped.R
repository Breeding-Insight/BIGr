#' Evaluate Pedigree File for Accuracy
#'
#' Check a pedigree file for accuracy and output suspected errors
#'
#'check_ped takes a 3-column pedigree tab separated file with columns labeled as id sire dam in any order and checks for:
#'* Ids that appear more than once in the id column
#'* Ids that appear in both sire and dam columns
#'* Direct (e.g. parent is a offspring of his own daughter) and indirect (e.g. a great grandparent is son of its grandchild) dependencies within the pedigree.
#'* Individuals included in the pedigree as sire or dam but not on the id column and reports them back with unknown parents (0).
#'
#'When using check_ped, do a first run to check for repeated ids and parents that appear as sire and dam.
#'Once these errors are cleaned run the function again to check for dependencies as this will provide the most accurate report.
#'
#'Note: This function does not change the input file but prints any errors found in the console.
#'
#' @param ped.file path to pedigree text file. The pedigree file is a
#' 3-column pedigree tab separated file with columns labeled as id sire dam in any order
#' @param seed Optional seed for reproducibility
#' @param verbose Logical. If TRUE, print the errors to the console.
#' @return A list of data.frames of error types, and the output printed to the console
#' @examples
#' ##Get list with a dataframe for each error type
#' ped_file <- system.file("check_ped_test.txt", package="BIGr")
#' ped_errors <- check_ped(ped.file = ped_file,
#'                         seed = 101919)
#'
#' ##Access the "messy parents" dataframe result
#' ped_errors$messy_parents
#'
#' ##Get list of sample IDs with messy parents error
#' messy_parent_ids <- ped_errors$messy_parents$id
#' print(messy_parent_ids)
#' @import dplyr
#' @import janitor
#' @importFrom stats setNames
#' @importFrom utils read.table
#' @export
check_ped <- function(ped.file, seed = NULL, verbose = TRUE) {
  #### Function to check for hierarchical errors missing parents and repeated ids ####
  if(!is.null(seed)){
    set.seed(seed)
  }
  #### read in data ####
  data = utils::read.table(ped.file, header = T)
  data <- data %>%
    janitor::clean_names() %>%
    mutate(
      id = as.character(id),
      sire = as.character(sire),
      dam = as.character(dam)
      )
  #Missing parents dataframe initialize
  missing_parents <- data.frame(id = character(), sire = character(), dam = character(), stringsAsFactors = FALSE)
  errors <- list()
  # repeated id checks
  n_occur <- data.frame(table(data$id))
  repeated_ids = n_occur[n_occur$Freq > 1,] %>%
    rename(id = Var1)
  # Check for ids that appear as both sire and dam ###This is possible for plants so maybe do not control for this or do not delete these rows just print them
  messy_parents <- as.data.frame(intersect(data$sire, data$dam)) %>%
    rename(id = 1) %>%
    filter(id != 0)
  # Missing parents check
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
  # Remove duplicates
  missing_parents <- distinct(missing_parents)
  # Combine original data with missing parents
  corrected_data <- bind_rows(data, missing_parents)
  # Function to detect cycles in the pedigree graph and identify the nodes involved
  detect_all_cycles <- function(data) {
    # Create an adjacency list
    adj_list <- list()
    for (i in 1:nrow(data)) {
      adj_list[[data$id[i]]] <- c(data$sire[i], data$dam[i])
    }
    # Helper function to perform DFS and detect cycles
    dfs <- function(node, visited, rec_stack, path) {
      visited[node] <- TRUE
      rec_stack[node] <- TRUE
      path <- append(path, node)
      cycles <- list()
      for (neighbor in adj_list[[node]]) {
        if (neighbor != "0") {
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
    # Initialize visited and recursion stack
    visited <- stats::setNames(rep(FALSE, length(adj_list)), names(adj_list))
    rec_stack <- stats::setNames(rep(FALSE, length(adj_list)), names(adj_list))
    all_cycles <- list()
    # Check for cycles in the graph and return the nodes involved
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
  # Check for cycles in the corrected pedigree data
  cycles <- detect_all_cycles(corrected_data)
  if (length(cycles) > 0) {
    cycle_number <- 1
    for (cycle_group in cycles) {
      cycle_ids <- unique(unlist(cycle_group))
      errors <- append(errors, paste("Cycle detected involving nodes:", paste(cycle_ids, collapse = " -> ")))
    }
  }
  results <- list(missing_parents = missing_parents, dependencies = data.frame(Dependency = unlist(errors)), repeated_ids = repeated_ids, messy_parents = messy_parents)
  repeated_ids <- results$repeated_ids
  missing_parents <- results$missing_parents
  messy_parents <- results$messy_parents
  errors <- results$dependencies
  # Adding the dataframes as an output list
  output.results <- list()
  #### Print errors and cycles ####
  # Print repeated ids if any
  if (nrow(repeated_ids) > 0) {
    if (verbose) {
      cat("Repeated ids found:\n")
      message(repeated_ids)
    }
    output.results$repeated_ids <- repeated_ids

  } else {
    if (verbose) {
      cat("No repeated ids found.\n")
    }
  }
  #Print parents that appear as male and female
  if (nrow(messy_parents) > 0) {
    if (verbose) {
      cat("Ids found as male and female parent:\n")
      message(messy_parents)
    }
    output.results$messy_parents <- messy_parents

  } else {
    if (verbose) {
      cat("No ids found as male and female parent.\n")
    }
  }
  # Print missing parents if any
  if (nrow(missing_parents) > 0) {
    if (verbose) {
      cat("Missing parents found:\n")
      message(missing_parents)
    }
    output.results$missing_parents <- missing_parents

  } else {
    if (verbose) {
      cat("No missing parents found.\n")
    }
  }
  # Print errors if any
  if (nrow(errors) > 0) {
    if (verbose) {
      cat("Dependencies found:\n")
      message(unique(errors$Dependency))
    }
    output.results$dependencies <- data.frame(Dependency = unlist(errors))

  } else {
    if (verbose) {
      cat("No dependencies found.\n")
    }
  }

  return(results)
}
