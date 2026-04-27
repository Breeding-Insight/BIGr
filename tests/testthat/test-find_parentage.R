# tests/testthat/test-find_parentage.R
# Run with: testthat::test_file("tests/testthat/test-find_parentage.R")
library(testthat)
library(data.table)

# ==============================================================================
# Helpers
# ==============================================================================

make_files <- function(genos, parents, progeny, dir = tempdir()) {
  geno_file    <- file.path(dir, paste0("genos_",   sample(1e6,1), ".txt"))
  parent_file  <- file.path(dir, paste0("parents_", sample(1e6,1), ".txt"))
  progeny_file <- file.path(dir, paste0("progeny_", sample(1e6,1), ".txt"))
  data.table::fwrite(genos,   geno_file,    sep = "\t")
  data.table::fwrite(parents, parent_file,  sep = "\t")
  data.table::fwrite(progeny, progeny_file, sep = "\t")
  list(g = geno_file, p = parent_file, pr = progeny_file)
}

# ==============================================================================
# Base toy data
# S1 / D1: all 0  → child1 (all 0) is a perfect Mendelian child of S1 x D1
# S2 / D2: all 2  → child2 (all 2) is a perfect Mendelian child of S2 x D2
# ==============================================================================

base_genos <- data.table::data.table(
  ID = c("S1","S2","D1","D2","child1","child2"),
  M1 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M2 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M3 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M4 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M5 = c(0L, 2L, 0L, 2L, 0L, 2L)
)
base_parents   <- data.table::data.table(ID  = c("S1","S2","D1","D2"),
                                         Sex = c("M","M","F","F"))
child1_progeny <- data.table::data.table(ID = "child1")
child2_progeny <- data.table::data.table(ID = "child2")
base_progeny   <- data.table::data.table(ID = c("child1","child2"))

# All-zero genotypes — every pair ties at 0% error
tied_genos <- data.table::data.table(
  ID = c("S1","S2","D1","D2","child_tie"),
  M1 = c(0L, 0L, 0L, 0L, 0L),
  M2 = c(0L, 0L, 0L, 0L, 0L)
)
tied_parents <- data.table::data.table(ID  = c("S1","S2","D1","D2"),
                                       Sex = c("M","M","F","F"))
tied_progeny <- data.table::data.table(ID = "child_tie")

# ==============================================================================
# 1. Input validation
# ==============================================================================

test_that("invalid method throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage(f$g, f$p, f$pr, method = "bad_method",
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
    regexp = "Method must be one of"
  )
})

test_that("min_markers < 1 throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage(f$g, f$p, f$pr, min_markers = 0,
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
    regexp = "min_markers"
  )
})

test_that("error_threshold out of range throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage(f$g, f$p, f$pr, error_threshold = 150,
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
    regexp = "error_threshold"
  )
  expect_error(
    find_parentage(f$g, f$p, f$pr, error_threshold = -1,
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
    regexp = "error_threshold"
  )
})

test_that("invalid na_string throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage(f$g, f$p, f$pr, na_string = "NULL",
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
    regexp = "na_string"
  )
})

test_that("missing genotype file throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage("nonexistent.txt", f$p, f$pr,
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
    regexp = "Error reading input files"
  )
})

test_that("parent IDs absent from genotype file raise a warning and are dropped", {
  extra_parents <- rbind(base_parents,
                         data.table::data.table(ID = "GHOST", Sex = "M"))
  f <- make_files(base_genos, extra_parents, child1_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
    regexp = "GHOST"
  )
})

test_that("progeny IDs absent from genotype file raise a warning and are dropped", {
  extra_progeny <- rbind(child1_progeny,
                         data.table::data.table(ID = "GHOST_KID"))
  f <- make_files(base_genos, base_parents, extra_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
    regexp = "GHOST_KID"
  )
})

test_that("no valid progeny candidates after filtering stops with an error", {
  ghost_progeny <- data.table::data.table(ID = "NOBODY")
  f <- make_files(base_genos, base_parents, ghost_progeny)
  expect_warning(
    expect_error(
      find_parentage(f$g, f$p, f$pr, method = "best_pair",
                     verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
      regexp = "No valid progeny"
    )
  )
})

test_that("missing Sex column raises a warning and defaults to ambiguous", {
  parents_no_sex <- data.table::data.table(ID = c("S1","D1"))
  f <- make_files(base_genos, parents_no_sex, child1_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best_match",
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE),
    regexp = "Sex"
  )
})

# ==============================================================================
# 2. Return structure
# ==============================================================================

test_that("best_pair returns a data.table with expected columns (no ties)", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny","Male_Parent","Female_Parent",
                    "Mendelian_Error_Pct","Markers_Tested",
                    "Assignment_Status") %in% names(res)))
  expect_equal(nrow(res), 1L)
})

test_that("best_male_parent returns expected columns", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_male_parent",
                        verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny","Best_Match","Mendelian_Error_Pct",
                    "Markers_Tested","Assignment_Status") %in% names(res)))
  expect_false("Male_Parent" %in% names(res))
  expect_equal(nrow(res), 1L)
})

test_that("best_female_parent returns expected columns", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_female_parent",
                        verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny","Best_Match","Mendelian_Error_Pct",
                    "Markers_Tested","Assignment_Status") %in% names(res)))
  expect_equal(nrow(res), 1L)
})

test_that("best_match returns expected columns", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_match",
                        verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny","Best_Match","Mendelian_Error_Pct",
                    "Markers_Tested","Assignment_Status") %in% names(res)))
  expect_equal(nrow(res), 1L)
})

test_that("one row returned per progeny for single-parent methods", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  for (m in c("best_male_parent","best_female_parent","best_match")) {
    res <- find_parentage(f$g, f$p, f$pr, method = m,
                          verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
    expect_equal(nrow(res), 1L, label = paste("row count for method", m))
  }
})

test_that("Markers_Tested equals the number of non-NA marker columns", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_equal(res$Markers_Tested, ncol(base_genos) - 1L)
})

test_that("Mendelian_Error_Pct is between 0 and 100", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  pct <- as.numeric(res$Mendelian_Error_Pct)
  expect_true(all(pct >= 0 & pct <= 100, na.rm = TRUE))
})

# ==============================================================================
# 3. Biological correctness
# ==============================================================================

test_that("best_pair correctly identifies S1 x D1 for child1 with 0% error", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_equal(res$Male_Parent,   "S1")
  expect_equal(res$Female_Parent, "D1")
  expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
})

test_that("best_pair correctly identifies S2 x D2 for child2 with 0% error", {
  f   <- make_files(base_genos, base_parents, child2_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_equal(res$Male_Parent,   "S2")
  expect_equal(res$Female_Parent, "D2")
  expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
})

test_that("best_male_parent identifies S1 as best male for child1", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_male_parent",
                        verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  expect_equal(res$Best_Match, "S1")
})

test_that("best_female_parent identifies D1 as best female for child1", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_female_parent",
                        verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  expect_equal(res$Best_Match, "D1")
})

test_that("Mendelian_Error_Pct is 0 for a perfect parent-progeny trio", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
})

test_that("both child1 and child2 are correctly assigned when run together", {
  f   <- make_files(base_genos, base_parents, base_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_equal(nrow(res), 2L)
  child1_row <- res[Progeny == "child1"]
  child2_row <- res[Progeny == "child2"]
  expect_equal(child1_row$Male_Parent,   "S1")
  expect_equal(child1_row$Female_Parent, "D1")
  expect_equal(child2_row$Male_Parent,   "S2")
  expect_equal(child2_row$Female_Parent, "D2")
})

# ==============================================================================
# 4. Assignment_Status — min_markers and error_threshold
# ==============================================================================

test_that("Assignment_Status = PASS for perfect trio within thresholds", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        min_markers = 3, error_threshold = 5.0,
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_equal(res$Assignment_Status, "PASS")
})

test_that("Assignment_Status = LOW_MARKERS when min_markers exceeds available markers", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        min_markers = 99999, error_threshold = 5.0,
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_equal(res$Assignment_Status, "LOW_MARKERS")
})

test_that("Assignment_Status = HIGH_ERROR when error rate exceeds threshold", {
  # Use wrong parents so error rate is high
  high_error_genos <- data.table::data.table(
    ID = c("S1","D1","bad_child"),
    M1 = c(0L, 0L, 2L),
    M2 = c(0L, 0L, 2L),
    M3 = c(0L, 0L, 2L),
    M4 = c(0L, 0L, 2L),
    M5 = c(0L, 0L, 2L)
  )
  parents  <- data.table::data.table(ID = c("S1","D1"), Sex = c("M","F"))
  progeny  <- data.table::data.table(ID = "bad_child")
  f        <- make_files(high_error_genos, parents, progeny)
  res      <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                             min_markers = 3, error_threshold = 5.0,
                             show_ties = FALSE, verbose = FALSE,
                             write_txt = FALSE, plot_results = FALSE)
  expect_equal(res$Assignment_Status, "HIGH_ERROR")
})

test_that("Assignment_Status column is present in all methods", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  for (m in c("best_pair","best_male_parent","best_female_parent","best_match")) {
    res <- find_parentage(f$g, f$p, f$pr, method = m,
                          show_ties = FALSE, verbose = FALSE,
                          write_txt = FALSE, plot_results = FALSE)
    expect_true("Assignment_Status" %in% names(res),
                label = paste("Assignment_Status present for method", m))
  }
})

test_that("LOW_MARKERS is flagged for single-parent methods too", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  for (m in c("best_male_parent","best_female_parent","best_match")) {
    res <- find_parentage(f$g, f$p, f$pr, method = m,
                          min_markers = 99999, verbose = FALSE,
                          write_txt = FALSE, plot_results = FALSE)
    expect_equal(res$Assignment_Status, "LOW_MARKERS",
                 label = paste("LOW_MARKERS for method", m))
  }
})

# ==============================================================================
# 5. allow_selfing
# ==============================================================================

test_that("allow_selfing = FALSE removes self-pairs from candidates", {
  ambig_parents <- data.table::data.table(ID = c("S1","D1"), Sex = c("A","A"))
  f <- make_files(base_genos, ambig_parents, child1_progeny)
  res <- suppressWarnings(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   allow_selfing = FALSE, show_ties = FALSE,
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  )
  if (!is.na(res$Male_Parent) && !is.na(res$Female_Parent))
    expect_false(res$Male_Parent == res$Female_Parent)
})

# ==============================================================================
# 6. show_ties
# ==============================================================================

test_that("show_ties = TRUE produces suffixed columns when ties exist", {
  f   <- make_files(tied_genos, tied_parents, tied_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = TRUE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_true(any(grepl("^Male_Parent_",   names(res))))
  expect_true(any(grepl("^Female_Parent_", names(res))))
})

test_that("show_ties = FALSE warns about ties and returns single-result columns", {
  f <- make_files(tied_genos, tied_parents, tied_progeny)
  expect_warning(
    res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                          show_ties = FALSE, verbose = FALSE,
                          write_txt = FALSE, plot_results = FALSE),
    regexp = "tied"
  )
  expect_true("Male_Parent"   %in% names(res))
  expect_true("Female_Parent" %in% names(res))
  expect_false(any(grepl("^Male_Parent_\\d",   names(res))))
  expect_false(any(grepl("^Female_Parent_\\d", names(res))))
})

test_that("base columns are always populated even when ties exist", {
  f   <- make_files(tied_genos, tied_parents, tied_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = TRUE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_false(is.na(res$Male_Parent[1]))
  expect_false(is.na(res$Female_Parent[1]))
})

# ==============================================================================
# 7. verbose / write_txt
# ==============================================================================

test_that("verbose = TRUE returns the result as data.table", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        verbose = TRUE, write_txt = FALSE, plot_results = FALSE)
  expect_s3_class(res, "data.table")
})

test_that("verbose = FALSE returns the result as data.table", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  expect_s3_class(res, "data.table")
})

test_that("write_txt = TRUE creates parentage_testing_results.txt", {
  old_wd <- getwd()
  tmp    <- tempdir()
  setwd(tmp)
  on.exit(setwd(old_wd), add = TRUE)
  f <- make_files(base_genos, base_parents, child1_progeny, dir = tmp)
  find_parentage(f$g, f$p, f$pr, method = "best_pair",
                 verbose = FALSE, write_txt = TRUE, plot_results = FALSE)
  expect_true(file.exists(file.path(tmp, "parentage_testing_results.txt")))
})

test_that("write_txt = FALSE does not create the output file", {
  old_wd   <- getwd()
  tmp      <- tempdir()
  setwd(tmp)
  on.exit(setwd(old_wd), add = TRUE)
  out_file <- file.path(tmp, "parentage_testing_results.txt")
  if (file.exists(out_file)) file.remove(out_file)
  f <- make_files(base_genos, base_parents, child1_progeny, dir = tmp)
  find_parentage(f$g, f$p, f$pr, method = "best_pair",
                 verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  expect_false(file.exists(out_file))
})

test_that("na_string = '' is accepted without error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_no_error(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   na_string = "", verbose = FALSE,
                   write_txt = FALSE, plot_results = FALSE)
  )
})

# ==============================================================================
# 8. Sex-based candidate filtering
# ==============================================================================

test_that("best_male_parent only assigns M or A parents", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_male_parent",
                        verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  valid_male_parents <- base_parents[Sex %in% c("M","A")]$ID
  expect_true(res$Best_Match %in% valid_male_parents)
})

test_that("best_female_parent only assigns F or A parents", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_female_parent",
                        verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  valid_female_parents <- base_parents[Sex %in% c("F","A")]$ID
  expect_true(res$Best_Match %in% valid_female_parents)
})

test_that("best_pair male slot contains only M or A parents", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  valid_males <- base_parents[Sex %in% c("M","A")]$ID
  expect_true(res$Male_Parent %in% valid_males)
})

test_that("best_pair female slot contains only F or A parents", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  valid_females <- base_parents[Sex %in% c("F","A")]$ID
  expect_true(res$Female_Parent %in% valid_females)
})

# ==============================================================================
# 9. Edge cases
# ==============================================================================

test_that("single progeny individual is handled correctly", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_equal(nrow(res), 1L)
})

test_that("all-NA marker column does not cause an error", {
  na_genos <- data.table::copy(base_genos)
  na_genos[, M1 := NA_integer_]
  f <- make_files(na_genos, base_parents, child1_progeny)
  expect_no_error(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, write_txt = FALSE, plot_results = FALSE)
  )
})

test_that("Progeny column contains the correct progeny IDs", {
  f   <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_setequal(res$Progeny, child1_progeny$ID)
})

test_that("multiple progeny are all represented in output", {
  f   <- make_files(base_genos, base_parents, base_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE,
                        write_txt = FALSE, plot_results = FALSE)
  expect_setequal(res$Progeny, base_progeny$ID)
})

test_that("single parent pair (n_pairs = 1) does not cause dimension errors", {
  single_pair_parents <- data.table::data.table(ID  = c("S1","D1"),
                                                Sex = c("M","F"))
  f   <- make_files(base_genos, single_pair_parents, child1_progeny)
  expect_no_error(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   show_ties = FALSE, verbose = FALSE,
                   write_txt = FALSE, plot_results = FALSE)
  )
})
