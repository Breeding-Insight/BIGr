#' Evaluate Pedigree File for Accuracy
#'
#' clean_pedigree takes a 3-column pedigree tab separated file with columns labeled as id sire
#' dam in any order and checks for ids that appear more than once in the id
#' column, ids that appear as both sire and dam, and circular dependencies where
#' an id is its own sire or dam.
#'
#' This function also looks for any sires or dams that
#' are not included in the pedigree in the id column and adds them with unknown parents.
#'
#' @param ped.file path to pedigree text file
#' @return A cleaned pedigree .txt file
#' @import dplyr
#' @export
#### Function to check for hierarchical errors missing parents and repeated ids ####
check_ped <- function(ped.file) {
  set.seed(101919)


  #### read in data ####
  data = read.table(ped.file, header = T)

  #Missing parents dataframe initialize
  missing_parents <- data.frame(id = character(), sire = character(), dam = character(), stringsAsFactors = FALSE)
  errors <- list()

  # repeated id checks
  n_occur <- data.frame(table(data$id))
  repeated_ids = n_occur[n_occur$Freq > 1,] %>%
    rename(id = Var1)

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
      errors <- append(errors, paste("Error: Individual", id, "cannot be its own parent"))
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
    visited <- setNames(rep(FALSE, length(adj_list)), names(adj_list))
    rec_stack <- setNames(rep(FALSE, length(adj_list)), names(adj_list))
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

  results <- list(missing_parents = missing_parents, errors = errors, repeated_ids = repeated_ids)

  repeated_ids <- results$repeated_ids
  missing_parents <- results$missing_parents
  errors <- results$errors

  #### Print errors and cycles ####
  # Print repeated ids if any
  if (nrow(repeated_ids) > 0) {
    cat("Repeated ids found:\n")
    print(repeated_ids)
  } else {
    cat("No repeated ids found.\n")
  }

  # Print missing parents if any
  if (nrow(missing_parents) > 0) {
    cat("Missing parents found:\n")
    print(missing_parents)
  } else {
    cat("No missing parents found.\n")
  }

  # Print errors if any
  if (length(errors) > 0) {
    cat("Dependencies found:\n")
    for (error in unique(errors)) {
      cat(error, "\n")
    }
  } else {
    cat("No dependencies found.\n")
  }

  return(results)
}
