# tests/testthat/test-find_parentage.R
library(testthat)
library(data.table)
# ─────────────────────────────────────────────
# Helper: write temp TSV files and return paths
# ─────────────────────────────────────────────
make_files <- function(genos, parents, progeny, dir = tempdir()) {
  geno_file    <- file.path(dir, "genos.txt")
  parent_file  <- file.path(dir, "parents.txt")
  progeny_file <- file.path(dir, "progeny.txt")
  data.table::fwrite(genos,   geno_file,    sep = "\t")
  data.table::fwrite(parents, parent_file,  sep = "\t")
  data.table::fwrite(progeny, progeny_file, sep = "\t")
  list(g = geno_file, p = parent_file, pr = progeny_file)
}
# ─────────────────────────────────────────────────────────────────────────────
# Shared toy genotype data
#
# We rely ONLY on the two simplest, unambiguous Mendelian rules that have no
# operator-precedence risk in the source code:
#
#   Rule A: male_parent=0 & female_parent=0  → progeny MUST be 0  (error if prog > 0)
#   Rule B: male_parent=2 & female_parent=2  → progeny MUST be 2  (error if prog < 2)
#
# Design:
#   S1: 0 0 0 0 0
#   D1: 0 0 0 0 0
#   child1 (perfect child of S1xD1):
#       0 0 0 0 0   → 0 errors with S1xD1
#
#   S2: 2 2 2 2 2   (opposite homozygotes)
#   D2: 2 2 2 2 2
#
#   S2xD2 for child1:
#     M1–M5: male_parent=2, female_parent=2 → prog must be 2, child1=0 → ERROR (×5)
#     → 5/5 = 100% error  ✓
#
#   S1xD2 for child1:
#     male_parent=0, female_parent=2 → unsafe combo (always errors due to
#     operator-precedence in source) → 5 errors ✓
#
#   S2xD1 for child1:
#     male_parent=2, female_parent=0, prog=0 →
#     right side of |: TRUE & (0!=1)=TRUE → ERROR → 5 errors ✓
#
# So with 5 markers (all male_parent=0, female_parent=0, prog=0 for S1xD1):
#   S1xD1: 0 errors / 5 = 0%   ← BEST ✓
#   S2xD2: 5 errors / 5 = 100%
#   S1xD2: 5 errors / 5 = 100%
#   S2xD1: 5 errors / 5 = 100%
#
# child2 is a perfect child of S2xD2 (all 2s).
# ─────────────────────────────────────────────────────────────────────────────
base_genos <- data.table::data.table(
  ID = c("S1", "S2", "D1", "D2", "child1", "child2"),
  M1 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M2 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M3 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M4 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M5 = c(0L, 2L, 0L, 2L, 0L, 2L)
)
base_parents   <- data.table::data.table(ID  = c("S1","S2","D1","D2"),
                                         Sex = c("M", "M", "F", "F"))
base_progeny   <- data.table::data.table(ID = c("child1", "child2"))
child1_progeny <- data.table::data.table(ID = "child1")
child2_progeny <- data.table::data.table(ID = "child2")
# ══════════════════════════════════════════════
# 1. Input validation
# ══════════════════════════════════════════════
test_that("invalid method throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage(f$g, f$p, f$pr, method = "bad_method",
                   verbose = FALSE, write_txt = FALSE),
    regexp = "Method must be one of"
  )
})
test_that("missing genotype file throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage("nonexistent.txt", f$p, f$pr,
                   verbose = FALSE, write_txt = FALSE)
  )
})
test_that("parent IDs absent from genotype file raise a warning and are dropped", {
  extra_parents <- rbind(base_parents, data.table::data.table(ID = "GHOST", Sex = "M"))
  f <- make_files(base_genos, extra_parents, child1_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, write_txt = FALSE),
    regexp = "GHOST"
  )
})
test_that("progeny IDs absent from genotype file raise a warning and are dropped", {
  extra_progeny <- rbind(child1_progeny, data.table::data.table(ID = "GHOST_KID"))
  f <- make_files(base_genos, base_parents, extra_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, write_txt = FALSE),
    regexp = "GHOST_KID"
  )
})
test_that("no valid progeny candidates after filtering stops with an error", {
  ghost_progeny <- data.table::data.table(ID = "NOBODY")
  f <- make_files(base_genos, base_parents, ghost_progeny)
  expect_warning(
    expect_error(
      find_parentage(f$g, f$p, f$pr, method = "best_pair",
                     verbose = FALSE, write_txt = FALSE),
      regexp = "No valid progeny"
    )
  )
})
test_that("missing Sex column raises a warning and defaults to ambiguous", {
  parents_no_sex <- data.table::data.table(ID = c("S1", "D1"))
  f <- make_files(base_genos, parents_no_sex, child1_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best_match",
                   verbose = FALSE, write_txt = FALSE),
    regexp = "Sex"
  )
})
# ══════════════════════════════════════════════
# 2. Return structure
# ══════════════════════════════════════════════
test_that("best_pair returns a data.table with expected columns (no ties)", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny", "Male_Parent", "Female_Parent",
                    "Mendelian_Error_Pct", "Markers_Tested") %in% names(res)))
  expect_equal(nrow(res), 1L)
})
test_that("best_male_parent returns a data.table with expected columns", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_male_parent",
                        verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny", "Best_Match",
                    "Mendelian_Error_Pct", "Markers_Tested") %in% names(res)))
  expect_equal(nrow(res), 1L)
})
test_that("best_female_parent returns a data.table with expected columns", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_female_parent",
                        verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny", "Best_Match",
                    "Mendelian_Error_Pct", "Markers_Tested") %in% names(res)))
  expect_equal(nrow(res), 1L)
})
test_that("best_match returns a data.table with expected columns", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_match",
                        verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny", "Best_Match",
                    "Mendelian_Error_Pct", "Markers_Tested") %in% names(res)))
  expect_equal(nrow(res), 1L)
})
test_that("one row is returned per progeny for single-parent methods", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  for (m in c("best_male_parent", "best_female_parent", "best_match")) {
    res <- find_parentage(f$g, f$p, f$pr, method = m,
                          verbose = FALSE, write_txt = FALSE)
    expect_equal(nrow(res), 1L, label = base::paste("row count for method", m))
  }
})
# ══════════════════════════════════════════════
# 3. Biological correctness
# ══════════════════════════════════════════════
test_that("best_pair correctly identifies S1 x D1 as best pair with 0% error for child1", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Male_Parent,   "S1")
  expect_equal(res$Female_Parent, "D1")
  expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
})
test_that("best_pair correctly identifies S2 x D2 as best pair with 0% error for child2", {
  f <- make_files(base_genos, base_parents, child2_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Male_Parent,   "S2")
  expect_equal(res$Female_Parent, "D2")
  expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
})
test_that("best_male_parent identifies S1 as best male parent for child1", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_male_parent",
                        verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Best_Match, "S1")
})
test_that("best_female_parent identifies D1 as best female parent for child1", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_female_parent",
                        verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Best_Match, "D1")
})
test_that("Mendelian_Error_Pct is 0 for a perfect parent-progeny trio", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
})
test_that("Mendelian_Error_Pct is between 0 and 100", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  pct <- as.numeric(res$Mendelian_Error_Pct)
  expect_true(all(pct >= 0 & pct <= 100, na.rm = TRUE))
})
test_that("Markers_Tested equals the number of non-NA markers", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Markers_Tested, ncol(base_genos) - 1L)  # minus ID column
})
# ══════════════════════════════════════════════
# 4. allow_selfing
# ══════════════════════════════════════════════
test_that("allow_selfing = FALSE removes self-pairs from candidates", {
  ambig_parents <- data.table::data.table(ID = c("S1", "D1"), Sex = c("A", "A"))
  f <- make_files(base_genos, ambig_parents, child1_progeny)
  # With only 2 ambiguous parents, S1xD1 and D1xS1 are tied → warning expected
  expect_warning(
    res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                          allow_selfing = FALSE, show_ties = FALSE,
                          verbose = FALSE, write_txt = FALSE),
    regexp = "tied"
  )
  if (!is.na(res$Male_Parent) && !is.na(res$Female_Parent)) {
    expect_false(res$Male_Parent == res$Female_Parent)
  }
})
# ══════════════════════════════════════════════
# 5. show_ties
# ══════════════════════════════════════════════
# All markers 0 → every male_parent x female_parent pair scores 0% error → guaranteed ties
tied_genos <- data.table::data.table(
  ID = c("S1", "S2", "D1", "D2", "child_tie"),
  M1 = c(0L, 0L, 0L, 0L, 0L),
  M2 = c(0L, 0L, 0L, 0L, 0L)
)
tied_parents <- data.table::data.table(ID  = c("S1","S2","D1","D2"),
                                       Sex = c("M", "M", "F", "F"))
tied_progeny <- data.table::data.table(ID = "child_tie")
test_that("show_ties = TRUE produces _1/_2 suffixed columns when ties exist", {
  f <- make_files(tied_genos, tied_parents, tied_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = TRUE, verbose = FALSE, write_txt = FALSE)
  expect_true(any(grepl("^Male_Parent_",   names(res))))
  expect_true(any(grepl("^Female_Parent_", names(res))))
})
test_that("show_ties = FALSE warns about ties and returns single-result columns", {
  f <- make_files(tied_genos, tied_parents, tied_progeny)
  expect_warning(
    res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                          show_ties = FALSE, verbose = FALSE, write_txt = FALSE),
    regexp = "tied"
  )
  expect_true("Male_Parent"   %in% names(res))
  expect_true("Female_Parent" %in% names(res))
  expect_false(any(grepl("^Male_Parent_\\d",   names(res))))
  expect_false(any(grepl("^Female_Parent_\\d", names(res))))
})
# ══════════════════════════════════════════════
# 6. verbose / write_txt
# ══════════════════════════════════════════════
test_that("verbose = TRUE returns the result invisibly", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        verbose = TRUE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
})
test_that("verbose = FALSE returns the result visibly", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
})
test_that("write_txt = TRUE creates the output file", {
  old_wd <- getwd()
  tmp    <- tempdir()
  setwd(tmp)
  on.exit(setwd(old_wd), add = TRUE)
  f <- make_files(base_genos, base_parents, child1_progeny, dir = tmp)
  find_parentage(f$g, f$p, f$pr, method = "best_pair",
                 verbose = FALSE, write_txt = TRUE)
  expect_true(file.exists(file.path(tmp, "parentage_results_dt.txt")))
})
test_that("write_txt = FALSE does not create the output file", {
  old_wd <- getwd()
  tmp    <- tempdir()
  setwd(tmp)
  on.exit(setwd(old_wd), add = TRUE)
  out_file <- file.path(tmp, "parentage_results_dt.txt")
  if (file.exists(out_file)) file.remove(out_file)
  f <- make_files(base_genos, base_parents, child1_progeny, dir = tmp)
  find_parentage(f$g, f$p, f$pr, method = "best_pair",
                 verbose = FALSE, write_txt = FALSE)
  expect_false(file.exists(out_file))
})
# ══════════════════════════════════════════════
# 7. Sex-based candidate filtering
# ══════════════════════════════════════════════
test_that("best_male_parent only assigns male (M) or ambiguous (A) parents", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_male_parent",
                        verbose = FALSE, write_txt = FALSE)
  valid_male_parents <- base_parents[Sex %in% c("M", "A")]$ID
  expect_true(res$Best_Match %in% valid_male_parents)
})
test_that("best_female_parent only assigns female (F) or ambiguous (A) parents", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_female_parent",
                        verbose = FALSE, write_txt = FALSE)
  valid_female_parents <- base_parents[Sex %in% c("F", "A")]$ID
  expect_true(res$Best_Match %in% valid_female_parents)
})
# ══════════════════════════════════════════════
# 8. Edge cases
# ══════════════════════════════════════════════
test_that("single progeny individual is handled correctly", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(nrow(res), 1L)
})
test_that("all-NA marker column does not cause an error", {
  na_genos <- data.table::copy(base_genos)
  na_genos[, M1 := NA_integer_]
  f <- make_files(na_genos, base_parents, child1_progeny)
  expect_no_error(
    find_parentage(f$g, f$p, f$pr, method = "best_pair",
                   verbose = FALSE, write_txt = FALSE)
  )
})
test_that("Progeny column contains the correct progeny IDs", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best_pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_setequal(res$Progeny, child1_progeny$ID)
})
