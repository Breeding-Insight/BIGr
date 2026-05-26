# tests/testthat/test-find_parentage.R
# Run with: testthat::test_file("tests/testthat/test-find_parentage.R")
library(testthat)
library(data.table)

# ==============================================================================
# Helpers
# ==============================================================================
make_files <- function(genos, parents, progeny, dir = tempdir()) {
  geno_file    <- file.path(dir, paste0("genos_",   sample(1e6, 1), ".txt"))
  parent_file  <- file.path(dir, paste0("parents_", sample(1e6, 1), ".txt"))
  progeny_file <- file.path(dir, paste0("progeny_", sample(1e6, 1), ".txt"))
  data.table::fwrite(genos,   geno_file,    sep = "\t")
  data.table::fwrite(parents, parent_file,  sep = "\t")
  data.table::fwrite(progeny, progeny_file, sep = "\t")
  list(g = geno_file, p = parent_file, pr = progeny_file)
}

# Base toy data — lowercase id column
# S1 / D1: all 0  → child1 (all 0) is a perfect Mendelian child of S1 x D1
# S2 / D2: all 2  → child2 (all 2) is a perfect Mendelian child of S2 x D2
base_genos <- data.table::data.table(
  id = c("S1", "S2", "D1", "D2", "child1", "child2"),
  M1 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M2 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M3 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M4 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M5 = c(0L, 2L, 0L, 2L, 0L, 2L)
)
base_parents   <- data.table::data.table(id  = c("S1", "S2", "D1", "D2"),
                                         sex = c("M", "M", "F", "F"))
child1_progeny <- data.table::data.table(id = "child1")
child2_progeny <- data.table::data.table(id = "child2")
base_progeny   <- data.table::data.table(id = c("child1", "child2"))

# All-zero genotypes — every pair ties at 0% error
tied_genos <- data.table::data.table(
  id = c("S1", "S2", "D1", "D2", "child_tie"),
  M1 = c(0L, 0L, 0L, 0L, 0L),
  M2 = c(0L, 0L, 0L, 0L, 0L)
)
tied_parents <- data.table::data.table(id  = c("S1", "S2", "D1", "D2"),
                                       sex = c("M", "M", "F", "F"))
tied_progeny <- data.table::data.table(id = "child_tie")

# ==============================================================================
# 1. Input validation
# ==============================================================================
test_that("invalid method throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage(f$g, f$p, f$pr, method = "bad_method",
                   verbose = FALSE, plot_results = FALSE),
    regexp = "Method must be one of"
  )
})

test_that("min_markers < 1 throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage(f$g, f$p, f$pr, min_markers = 0,
                   verbose = FALSE, plot_results = FALSE),
    regexp = "min_markers"
  )
})

test_that("error_threshold out of range throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage(f$g, f$p, f$pr, error_threshold = 150,
                   verbose = FALSE, plot_results = FALSE),
    regexp = "error_threshold"
  )
  expect_error(
    find_parentage(f$g, f$p, f$pr, error_threshold = -1,
                   verbose = FALSE, plot_results = FALSE),
    regexp = "error_threshold"
  )
})

test_that("missing genotype file throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage("nonexistent.txt", f$p, f$pr,
                   verbose = FALSE, plot_results = FALSE),
    regexp = "Error reading input files"
  )
})

test_that("parent IDs absent from genotype file raise a warning and are dropped", {
  extra_parents <- rbind(base_parents,
                         data.table::data.table(id = "GHOST", sex = "M"))
  f <- make_files(base_genos, extra_parents, child1_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, plot_results = FALSE),
    regexp = "GHOST"
  )
})

test_that("progeny IDs absent from genotype file raise a warning and are dropped", {
  extra_progeny <- rbind(child1_progeny,
                         data.table::data.table(id = "GHOST_KID"))
  f <- make_files(base_genos, base_parents, extra_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, plot_results = FALSE),
    regexp = "GHOST_KID"
  )
})

test_that("no valid progeny candidates after filtering stops with an error", {
  ghost_progeny <- data.table::data.table(id = "NOBODY")
  f <- make_files(base_genos, base_parents, ghost_progeny)
  expect_warning(
    expect_error(
      find_parentage(f$g, f$p, f$pr, method = "best_pair",
                     verbose = FALSE, plot_results = FALSE),
      regexp = "No valid progeny"
    )
  )
})

test_that("missing sex column raises a warning and defaults to ambiguous", {
  parents_no_sex <- data.table::data.table(id = c("S1", "D1"))
  f <- make_files(base_genos, parents_no_sex, child1_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best_match",
                   verbose = FALSE, plot_results = FALSE),
    regexp = "sex"
  )
})

# ==============================================================================
# 2. Return structure — named list
# ==============================================================================
test_that("find_parentage returns an invisible named list with all required elements", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_type(out, "list")
  expect_named(out, c("pass", "high_error", "low_markers", "full_results", "plot"),
               ignore.order = TRUE)
})

test_that("full_results is a data.table", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_s3_class(out$full_results, "data.table")
})

test_that("full_results has one row per progeny", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_equal(nrow(out$full_results), 1L)
})

test_that("named list subsets together cover all full_results rows", {
  f   <- make_files(base_genos, base_parents, base_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  subset_total <- nrow(out$pass) + nrow(out$high_error) + nrow(out$low_markers)
  expect_equal(subset_total, nrow(out$full_results))
})

test_that("plot element is NULL when plot_results = FALSE", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        verbose = FALSE, plot_results = FALSE)

  expect_null(out$plot)
})

test_that("best_pair full_results has expected lowercase columns", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_true(all(c("id", "male_parent", "female_parent",
                    "mendelian_error_pct", "markers_tested",
                    "status") %in% names(out$full_results)))
})

test_that("best_male_parent full_results has expected lowercase columns", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_male_parent",
                        verbose = FALSE, plot_results = FALSE)

  expect_true(all(c("id", "best_match", "mendelian_error_pct",
                    "markers_tested", "status") %in% names(out$full_results)))
  expect_false("male_parent" %in% names(out$full_results))
})

test_that("best_female_parent full_results has expected lowercase columns", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_female_parent",
                        verbose = FALSE, plot_results = FALSE)

  expect_true(all(c("id", "best_match", "mendelian_error_pct",
                    "markers_tested", "status") %in% names(out$full_results)))
})

test_that("best_match full_results has expected lowercase columns", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_match",
                        verbose = FALSE, plot_results = FALSE)

  expect_true(all(c("id", "best_match", "mendelian_error_pct",
                    "markers_tested", "status") %in% names(out$full_results)))
})

# ==============================================================================
# 3. Biological correctness
# ==============================================================================
test_that("best_pair correctly identifies S1 x D1 for child1 with 0% error", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results$male_parent,   "S1")
  expect_equal(out$full_results$female_parent, "D1")
  expect_equal(as.numeric(out$full_results$mendelian_error_pct), 0)
})

test_that("best_pair correctly identifies S2 x D2 for child2 with 0% error", {
  f   <- make_files(base_genos, base_parents, child2_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results$male_parent,   "S2")
  expect_equal(out$full_results$female_parent, "D2")
  expect_equal(as.numeric(out$full_results$mendelian_error_pct), 0)
})

test_that("best_male_parent identifies S1 as best male for child1", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_male_parent",
                        verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results$best_match, "S1")
})

test_that("best_female_parent identifies D1 as best female for child1", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_female_parent",
                        verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results$best_match, "D1")
})

test_that("both child1 and child2 correctly assigned when run together", {
  f   <- make_files(base_genos, base_parents, base_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_equal(nrow(out$full_results), 2L)
  expect_equal(out$full_results[id == "child1"]$male_parent,   "S1")
  expect_equal(out$full_results[id == "child1"]$female_parent, "D1")
  expect_equal(out$full_results[id == "child2"]$male_parent,   "S2")
  expect_equal(out$full_results[id == "child2"]$female_parent, "D2")
})

test_that("markers_tested equals number of non-NA marker columns", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results$markers_tested, ncol(base_genos) - 1L)
})

test_that("mendelian_error_pct is between 0 and 100", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  pct <- as.numeric(out$full_results$mendelian_error_pct)
  expect_true(all(pct >= 0 & pct <= 100, na.rm = TRUE))
})

# ==============================================================================
# 4. status — lowercase values
# ==============================================================================
test_that("status = pass for perfect trio within thresholds", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        min_markers = 3, error_threshold = 5.0,
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results$status, "pass")
})

test_that("status = low_markers when min_markers exceeds available markers", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        min_markers = 99999, error_threshold = 5.0,
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results$status, "low_markers")
})

test_that("status = high_error when error rate exceeds threshold", {
  high_error_genos <- data.table::data.table(
    id = c("S1", "D1", "bad_child"),
    M1 = c(0L, 0L, 2L),
    M2 = c(0L, 0L, 2L),
    M3 = c(0L, 0L, 2L),
    M4 = c(0L, 0L, 2L),
    M5 = c(0L, 0L, 2L)
  )
  parents <- data.table::data.table(id = c("S1", "D1"), sex = c("M", "F"))
  progeny <- data.table::data.table(id = "bad_child")
  f       <- make_files(high_error_genos, parents, progeny)
  out     <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                            min_markers = 3, error_threshold = 5.0,
                            show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results$status, "high_error")
})

test_that("status column is present and lowercase in all methods", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  for (m in c("best_pair", "best_male_parent", "best_female_parent", "best_match")) {
    out <- find_parentage(f$g, f$p, f$pr, method = m,
                          show_ties = FALSE, verbose = FALSE, plot_results = FALSE)
    expect_true("status" %in% names(out$full_results),
                label = paste("status column present for method", m))
    expect_true(all(out$full_results$status %in%
                      c("pass", "high_error", "low_markers", NA)),
                label = paste("lowercase status values for method", m))
  }
})

test_that("low_markers is flagged for single-parent methods too", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  for (m in c("best_male_parent", "best_female_parent", "best_match")) {
    out <- find_parentage(f$g, f$p, f$pr, method = m,
                          min_markers = 99999, verbose = FALSE,
                          plot_results = FALSE)
    expect_equal(out$full_results$status, "low_markers",
                 label = paste("low_markers for method", m))
  }
})

test_that("pass list element contains only pass rows", {
  f   <- make_files(base_genos, base_parents, base_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  if (nrow(out$pass) > 0)
    expect_true(all(out$pass$status == "pass"))
})

test_that("high_error list element contains only high_error rows", {
  high_error_genos <- data.table::data.table(
    id = c("S1", "D1", "bad_child"),
    M1 = c(0L, 0L, 2L), M2 = c(0L, 0L, 2L), M3 = c(0L, 0L, 2L),
    M4 = c(0L, 0L, 2L), M5 = c(0L, 0L, 2L)
  )
  parents <- data.table::data.table(id = c("S1", "D1"), sex = c("M", "F"))
  progeny <- data.table::data.table(id = "bad_child")
  f       <- make_files(high_error_genos, parents, progeny)
  out     <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                            min_markers = 3, error_threshold = 5.0,
                            show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  if (nrow(out$high_error) > 0)
    expect_true(all(out$high_error$status == "high_error"))
})

# ==============================================================================
# 5. allow_selfing
# ==============================================================================
test_that("allow_selfing = FALSE removes self-pairs from candidates", {
  ambig_parents <- data.table::data.table(id = c("S1", "D1"), sex = c("A", "A"))
  f <- make_files(base_genos, ambig_parents, child1_progeny)
  out <- suppressWarnings(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   allow_selfing = FALSE, show_ties = FALSE,
                   verbose = FALSE, plot_results = FALSE)
  )
  r <- out$full_results
  if (!is.na(r$male_parent) && !is.na(r$female_parent))
    expect_false(r$male_parent == r$female_parent)
})

# ==============================================================================
# 6. show_ties
# ==============================================================================
test_that("show_ties = TRUE produces lowercase suffixed columns when ties exist", {
  f   <- make_files(tied_genos, tied_parents, tied_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = TRUE, verbose = FALSE, plot_results = FALSE)

  expect_true(any(grepl("^male_parent_",   names(out$full_results))))
  expect_true(any(grepl("^female_parent_", names(out$full_results))))
})

test_that("show_ties = FALSE warns about ties and returns single-result columns", {
  f <- make_files(tied_genos, tied_parents, tied_progeny)
  expect_warning(
    out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                          show_ties = FALSE, verbose = FALSE, plot_results = FALSE),
    regexp = "tied"
  )
  expect_true("male_parent"   %in% names(out$full_results))
  expect_true("female_parent" %in% names(out$full_results))
  expect_false(any(grepl("^male_parent_\\d",   names(out$full_results))))
  expect_false(any(grepl("^female_parent_\\d", names(out$full_results))))
})

test_that("base columns are always populated even when ties exist", {
  f   <- make_files(tied_genos, tied_parents, tied_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = TRUE, verbose = FALSE, plot_results = FALSE)

  expect_false(is.na(out$full_results$male_parent[1]))
  expect_false(is.na(out$full_results$female_parent[1]))
})

# ==============================================================================
# 7. verbose
# ==============================================================================
test_that("verbose = TRUE returns the result as a named list", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        verbose = TRUE, plot_results = FALSE)
  expect_type(out, "list")
})

test_that("verbose = FALSE suppresses console output", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_silent(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, plot_results = FALSE)
  )
})

# ==============================================================================
# 8. No write logic — find_parentage does not write files
# ==============================================================================
test_that("no output files are written to disk", {
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

  f <- make_files(base_genos, base_parents, child1_progeny, dir = tmp_dir)
  find_parentage(f$g, f$p, f$pr, method = "best_pair",
                 verbose = FALSE, plot_results = FALSE)

  written_files <- list.files(tmp_dir, pattern = "\\.txt$|\\.jpg$|\\.csv$")
  # Only the input files we created should be there
  expect_equal(length(written_files), 3L)
})

# ==============================================================================
# 9. Sex-based candidate filtering
# ==============================================================================
test_that("best_male_parent only assigns M or A parents", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_male_parent",
                        verbose = FALSE, plot_results = FALSE)

  valid_male_parents <- base_parents[sex %in% c("M", "A")]$id
  expect_true(out$full_results$best_match %in% valid_male_parents)
})

test_that("best_female_parent only assigns F or A parents", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_female_parent",
                        verbose = FALSE, plot_results = FALSE)

  valid_female_parents <- base_parents[sex %in% c("F", "A")]$id
  expect_true(out$full_results$best_match %in% valid_female_parents)
})

test_that("best_pair male slot contains only M or A parents", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  valid_males <- base_parents[sex %in% c("M", "A")]$id
  expect_true(out$full_results$male_parent %in% valid_males)
})

test_that("best_pair female slot contains only F or A parents", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)

  valid_females <- base_parents[sex %in% c("F", "A")]$id
  expect_true(out$full_results$female_parent %in% valid_females)
})

# ==============================================================================
# 10. Edge cases
# ==============================================================================
test_that("single progeny individual is handled correctly", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)
  expect_equal(nrow(out$full_results), 1L)
})

test_that("all-NA marker column does not cause an error", {
  na_genos <- data.table::copy(base_genos)
  na_genos[, M1 := NA_integer_]
  f <- make_files(na_genos, base_parents, child1_progeny)
  expect_no_error(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, plot_results = FALSE)
  )
})

test_that("id column in full_results contains the correct progeny IDs", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)
  expect_setequal(out$full_results$id, child1_progeny$id)
})

test_that("multiple progeny are all represented in full_results", {
  f   <- make_files(base_genos, base_parents, base_progeny)
  out <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, plot_results = FALSE)
  expect_setequal(out$full_results$id, base_progeny$id)
})

test_that("single parent pair (n_pairs = 1) does not cause dimension errors", {
  single_pair_parents <- data.table::data.table(id  = c("S1", "D1"),
                                                sex = c("M", "F"))
  f <- make_files(base_genos, single_pair_parents, child1_progeny)
  expect_no_error(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   show_ties = FALSE, verbose = FALSE, plot_results = FALSE)
  )
})

test_that("one row returned per progeny for single-parent methods", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  for (m in c("best_male_parent", "best_female_parent", "best_match")) {
    out <- find_parentage(f$g, f$p, f$pr, method = m,
                          verbose = FALSE, plot_results = FALSE)
    expect_equal(nrow(out$full_results), 1L,
                 label = paste("row count for method", m))
  }
})

# ==============================================================================
# 11. plot element
# ==============================================================================
test_that("plot element is a ggplot when plot_results = TRUE", {
  skip_if_not_installed("ggplot2")
  f   <- make_files(base_genos, base_parents, child1_progeny)
  out <- suppressWarnings(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, plot_results = TRUE)
  )
  expect_s3_class(out$plot, "ggplot")
})

# ==============================================================================
# 12. Return value is invisible
# ==============================================================================
test_that("find_parentage returns invisibly", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_invisible(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, plot_results = FALSE)
  )
})
