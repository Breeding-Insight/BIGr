#' Evaluate Pedigree File for Accuracy
#'
#' clean_pedigree takes a 3-column pedigree tab separated file with columns labeled as id sire
#' dam in any order and checks for ids that appear more than once in the id
#' column, ids that appear as both sire and dam and circular dependencies where
#' an id is its own sire or dam.
#'
#' This function also looks for any sires or dams that
#' are not included in the pedigree in the id column and adds them with unknown parents.
#'
#' @param input_ped path to pedigree text file
#' @param save.result save the cleaned pedigree file (True/False) (default = TRUE)
#' @return A cleaned pedigree .txt file
#' @import dplyr
#' @importFrom kinship2 fixParents
#' @importFrom kinship2 pedigree
#' @export
clean_pedigree <- function(input_ped, save.result = TRUE) {
  #Might need tidyverse and not just dplyr
  set.seed(101919)
  # Read in pedigree with cols id, dam, sire in any order
  raw_ped <- read.table(input_ped, header = TRUE) %>%
    mutate(id = as.factor(id),
           sire = as.factor(sire),
           dam = as.factor(dam))

  # Add sex column (male 0, female 1, unknown 2)
  sex_ped <- raw_ped %>%
    mutate(sex = ifelse(id %in% sire, 0, ifelse(id %in% dam, 1, 2)))

  # Check for ids that appear as both sire and dam
  messy_parents <- as.data.frame(intersect(sex_ped$sire, sex_ped$dam)) %>%
    rename(id = 1) %>%
    filter(id != 0)
  print("IDs appearing as sire and dam")
  print(messy_parents)

  # Change values in messy_parents to 0 in cols for sire and dam
  # Create new dataframe
  parents_fixed_ped <- sex_ped

  # Fix values in sire column
  parents_fixed_ped$sire[parents_fixed_ped$sire %in% messy_parents[, 1]] <- 0
  # Fix values in dam column
  parents_fixed_ped$dam[parents_fixed_ped$dam %in% messy_parents[, 1]] <- 0

  # Check for duplicate IDs
  n_occur <- data.frame(table(parents_fixed_ped$id)) %>%
    rename(id = 1, freq = 2) # Count occurrence of each ID

  double_ids <- n_occur[n_occur$freq > 1,] # Print ids with more than 1 occurrence

  print("IDs that appear more than once in the id column")
  print(double_ids)

  # Create new dataframe without doubled ids (all instances are removed to avoid data errors)
  nodup_ped <- parents_fixed_ped[!parents_fixed_ped$id %in% double_ids$id, ]

  # Check for individuals that appear as offspring and sire or offspring and dam of themselves
  circ_deps <- nodup_ped %>%
    mutate(id = as.character(id),
           sire = as.character(sire),
           dam = as.character(dam)) %>%
    filter(id == sire | id == dam)

  print("Circular dependencies")

  print(circ_deps)

  no_circ_deps <- anti_join(nodup_ped, circ_deps, by = "id")

  # Run FixPedigree to add parents that are not found as founders
  ready_ped <- with(no_circ_deps, kinship2::fixParents(id, sire, dam, sex, missid="0"))

  # Create a pedigree object
  final_ped <- with(ready_ped, kinship2::pedigree(id, dadid, momid, sex, missid="0"))

  # Write out clean pedigree
  if (save.result){
    output_file <- paste0("cleaned_",basename(input_ped))
    write.table(final_ped, file = output_file, row.names = FALSE, quote = FALSE)
  }
  return(final_ped)
}
