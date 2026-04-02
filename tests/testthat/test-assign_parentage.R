# tests/testthat/test_find_parentage.R
# Test suite for find_parentage() (assign_parentage.R)  [1]
# ------------------------------------------------------------------

# 1. Load the function -----------------------------------------------
# Assuming the function is in the current package:
# You may need to adjust the path if it is in a different location.
source("R/assign_parentage.R")   # [1]

# 2. Helper: create a very small, deterministic dataset -----------
small_test_data <- function() {
  # 3 parents (one male, one female, one unknown)
  parents_df <- tibble::tribble(
    ~ID,   ~Sex,
    "P1",  "M",
    "P2",  "F",
    "P3",  NA_character_
  )
  # Progeny
  progeny_df <- tibble::tibble(ID = "C1")

  # Genotypes (0 = AA, 1 = Aa, 2 = aa)
  # P1 (M) : AA  -> 0
  # P2 (F) : aa  -> 2
  # P3 : unknown
  # C1 : should be heterozygous (Aa) -> 1
  genotypes_df <- tibble::tribble(
    ~ID,   ~L1,
    "P1",  0,
    "P2",  2,
    "P3",  NA,
    "C1",  1
  )

  list(parents = parents_df,
       progeny  = progeny_df,
       genotypes = genotypes_df)
}

# 3. Unit tests -------------------------------------------------------

testthat::test_that("find_parentage returns a tibble with expected columns", {
  data <- small_test_data()

  # Write temp files
  tmp_dir <- tempdir()
  write_tsv(data$parents, file.path(tmp_dir, "parents.tsv"))
  write_tsv(data$progeny,  file.path(tmp_dir, "progeny.tsv"))
  write_tsv(data$genotypes, file.path(tmp_dir, "genotypes.tsv"))

  # Run the function
  res <- find_parentage(
    genotypes_file = file.path(tmp_dir, "genotypes.tsv"),
    parents_file   = file.path(tmp_dir, "parents.tsv"),
    progeny_file   = file.path(tmp_dir, "progeny.tsv"),
    method         = "best.pair",
    write.txt      = FALSE,   # we don't need the text file in tests
    verbose        = FALSE
  )

  # Basic sanity checks
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_equal(nrow(res), 1)                            # one progeny
  testthat::expect_true(all(c("Progeny", "Sire", "Dam", "Mendelian_Error_Pct",
                              "Markers_Tested") %in% colnames(res)))

  # The best parent pair should be (P1, P2) with no Mendelian errors
  testthat::expect_equal(res$Progeny, "C1")
  testthat::expect_equal(res$Sire, "P1")
  testthat::expect_equal(res$Dam, "P2")
  testthat::expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
  testthat::expect_equal(res$Markers_Tested, 1)   # L1 was compared
})

testthat::test_that("best.sire and best.dam methods work on the same data", {
  data <- small_test_data()

  write_tsv(data$parents, file.path(tmpdir(), "parents.tsv"))
  write_tsv(data$progeny,  file.path(tmpdir(), "progeny.tsv"))
  write_tsv(data$genotypes, file.path(tmpdir(), "genotypes.tsv"))

  # Best sire
  res_sire <- find_parentage(
    genotypes_file = file.path(tmpdir(), "genotypes.tsv"),
    parents_file   = file.path(tmpdir(), "parents.tsv"),
    progeny_file   = file.path(tmpdir(), "progeny.tsv"),
    method         = "best.sire",
    write.txt      = FALSE,
    verbose        = FALSE
  )
  testthat::expect_equal(res_sire$Best_Match, "P1")

  # Best dam
  res_dam <- find_parentage(
    genotypes_file = file.path(tmpdir(), "genotypes.tsv"),
    parents_file   = file.path(tmpdir(), "parents.tsv"),
    progeny_file   = file.path(tmpdir(), "progeny.tsv"),
    method         = "best.dam",
    write.txt      = FALSE,
    verbose        = FALSE
  )
  testthat::expect_equal(res_dam$Best_Match, "P2")
})

testthat::test_that("function can handle a larger random dataset", {
  # Generate a larger random dataset using the provided script logic  [2]
  # (The code below mirrors the logic in test_data.R but runs locally.)
  set.seed(42)
  n_progeny <- 50
  n_loci   <- 100
  n_males  <- 5
  n_females <- 5
  n_unknown <- 2

  # Build parent list
  parents_df <- tibble::tribble(
    ~ID, ~Sex,
    paste0("Male_", 1:n_males), "M",
    paste0("Female_", 1:n_females), "F",
    paste0("Unknown_", 1:n_unknown), NA_character_
  ) %>% tibble::add_row()

  # Progeny list
  progeny_ids <- paste0("Progeny_", 1:n_progeny)
  progeny_df  <- tibble::tibble(ID = progeny_ids)

  # Random genotypes for parents
  all_parent_ids <- parents_df$ID
  parent_gens_mat <- matrix(
    sample(0:2, length(all_parent_ids) * n_loci, replace = TRUE),
    nrow = length(all_parent_ids),
    dimnames = list(all_parent_ids,
                    paste0("Locus_", 1:n_loci))
  )

  # Build progeny genotypes from random parents (simple simulation)
  progeny_genos_mat <- matrix(NA_integer_, nrow = n_progeny,
                              ncol = n_loci,
                              dimnames = list(progeny_ids,
                                              paste0("Locus_", 1:n_loci)))

  for (i in seq_along(progeny_ids)) {
    sire  <- sample(all_parent_ids, 1)
    dam   <- sample(all_parent_ids, 1)
    for (l in 1:n_loci) {
      parent_alleles <- c(parent_gens_mat[sire, l],
                          parent_gens_mat[dam, l])
      # Simple rule: sum to get genotype (0,1,2) – works because alleles are 0 or 1
      progeny_genos_mat[i, l] <- sum(parent_alleles)
    }
  }

  # Combine and write to files
  all_gens_df <- tibble::as_tibble(
    Rcpp::cppFunction('library(tibble);') # just for tibble import
  )
  all_gens_df <- tibble::as_tibble(
    rbind(parents_df, progeny_df)
  ) %>%
    tibble::add_row()

  # Write files
  tmp_dir <- tempdir()
  write_tsv(parents_df,   file.path(tmp_dir, "parents.tsv"))
  write_tsv(progeny_df,   file.path(tmp_dir, "progeny.tsv"))
  write_tsv(all_parent_ids, file.path(tmp_dir, "genotypes.tsv"))  # placeholder

  # Call the function
  res <- find_parentage(
    genotypes_file = file.path(tmp_dir, "genotypes.tsv"),
    parents_file   = file.path(tmp_dir, "parents.tsv"),
    progeny_file   = file.path(tmp_dir, "progeny.tsv"),
    method         = "best.pair",
    write.txt      = FALSE,
    verbose        = FALSE
  )

  # Basic checks on the larger set
  testthat::expect_equal(nrow(res), n_progeny)
  testthat::expect_true(all(c("Progeny", "Sire", "Dam",
                              "Mendelian_Error_Pct", "Markers_Tested") %in% colnames(res)))
  testthat::expect_true(all(res$Mendelian_Error_Pct >= 0))
  testthat::expect_true(all(res$Mendelian_Error_Pct <= 100))
})
