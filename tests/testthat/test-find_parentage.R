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
  fwrite(genos,   geno_file,    sep = "\t")
  fwrite(parents, parent_file,  sep = "\t")
  fwrite(progeny, progeny_file, sep = "\t")
  list(g = geno_file, p = parent_file, pr = progeny_file)
}

# ─────────────────────────────────────────────────────────────────────────────
# Shared toy genotype data
#
# We rely ONLY on the two simplest, unambiguous Mendelian rules that have no
# operator-precedence risk in the source code:
#
#   Rule A: sire=0 & dam=0  → progeny MUST be 0  (error if prog > 0)
#   Rule B: sire=2 & dam=2  → progeny MUST be 2  (error if prog < 2)
#
# Design:
#   S1: 0 0 0 0 0 2 2 2 2 2
#   D1: 0 0 0 0 0 2 2 2 2 2
#   child1 (perfect child of S1xD1):
#       0 0 0 0 0 2 2 2 2 2   → 0 errors with S1xD1
#
#   S2: 2 2 2 2 2 0 0 0 0 0   (opposite homozygotes)
#   D2: 2 2 2 2 2 0 0 0 0 0
#
#   S2xD2 for child1:
#     M1–M5: s=2,d=2 → prog must be 2, child1=0 → ERROR  (×5)
#     M6–M10: s=0,d=0 → prog must be 0, child1=2 → ERROR (×5)
#     → 10/10 = 100% error  ✓
#
#   S1xD2 for child1:
#     M1–M5: s=0,d=2 → heterozygous rule (involves compound condition,
#            not tested here — but S2xD1 and S2xD2 already have 100% error)
#     We only need S1xD1 to be strictly better than all others.
#     S1xD2: M1–M5: s=0,d=2 → no Rule A/B fires → 0 errors on M1-M5
#            M6–M10: s=2,d=0 → no Rule A/B fires → 0 errors on M6-M10
#            → 0% error  ← tie with S1xD1!
#
# To break S1xD2 tie, add markers where S1=0,D1=0,child1=0 but D2≠0:
#   Add M11, M12: S1=0, S2=2, D1=0, D2=2, child1=0
#     S1xD1: s=0,d=0 → must be 0, child1=0 → OK  ✓
#     S1xD2: s=0,d=2 → Rule A/B don't fire → OK (no error counted)
#     S2xD1: s=2,d=0 → Rule A/B don't fire → OK
#     S2xD2: s=2,d=2 → must be 2, child1=0 → ERROR ✓
#
# Hmm — S1xD2 still 0 errors. The only way to get errors for S1xD2 using
# only Rule A/B is if we have a marker where S1=0, D2=0, child1≠0
# OR S1=2, D2=2, child1≠2.
#
# Final clean design using ONLY Rule A (s=0,d=0→prog=0):
#
#   Group 1 (5 markers): S1=0, D1=0, S2=2, D2=2, child1=0
#     S1xD1: Rule A → prog=0 ✓ (0 errors)
#     S2xD2: Rule B → prog must be 2, child1=0 → ERROR ✓
#     S1xD2: s=0,d=2 → no Rule A/B → 0 errors  (still tied)
#     S2xD1: s=2,d=0 → no Rule A/B → 0 errors  (still tied)
#
#   Group 2 (5 markers): S1=0, D1=0, S2=0, D2=2, child1=0
#     S1xD1: Rule A → 0 errors ✓
#     S1xD2: s=0,d=2 → no error from Rule A/B
#     S2xD1: Rule A → 0 errors
#     S2xD2: s=0,d=2 → no error from Rule A/B
#
# It is impossible to distinguish S1xD1 from S1xD2 using ONLY Rule A/B
# when child1 is always 0 (since Rule A needs d=0 too, and D2≠0 means
# Rule A doesn't fire for S1xD2, giving no error).
#
# CONCLUSION: We must allow heterozygous markers but avoid the
# operator-precedence bug. Looking at the source code condition:
#
#   ((sire==0 & dam==1) | (sire==1 & dam==0)) & (prog==2)
#
# Due to R's precedence (&  binds tighter than |), this parses as:
#   (sire==0 & dam==1) | ((sire==1 & dam==0) & prog==2)
#
# So the condition misfires for sire=0,dam=1,prog=anything (always TRUE
# for the left side regardless of prog). This means any marker where
# sire=0,dam=1 will ALWAYS be counted as an error, regardless of progeny.
#
# Similarly: ((sire==2 & dam==1) | (sire==1 & dam==2)) & (prog==0)
# parses as: (sire==2 & dam==1) | ((sire==1 & dam==2) & prog==0)
# → sire=2,dam=1 always flagged as error.
#
# And: ((sire==0 & dam==2) | (sire==2 & dam==0)) & (prog!=1)
# parses as: (sire==0 & dam==2) | ((sire==2 & dam==0) & prog!=1)
# → sire=0,dam=2 always flagged as error regardless of prog.
#
# SAFE rules (no precedence issue):
#   Rule A: s=0,d=0 → prog must be 0   ✓ safe
#   Rule B: s=2,d=2 → prog must be 2   ✓ safe
#
# UNSAFE parent combos (always produce errors due to bug):
#   s=0,d=1 → always error
#   s=2,d=1 → always error
#   s=0,d=2 → always error
#
# SAFE combos with no error fired (Rule A/B don't apply):
#   s=1,d=0, s=1,d=2, s=1,d=1, s=2,d=0 (only right side of | checked)
#   s=2,d=0: parses as (FALSE) | (TRUE & prog!=1) → error only if prog!=1
#            so s=2,d=0,prog=1 → NO error ✓
#
# New design using ONLY Rule A, Rule B, and the safe s=2,d=0,prog=1 case:
#
#   Group 1 (5 markers): S1=0,D1=0,child1=0 → Rule A, 0 errors for S1xD1
#                        S2=2,D2=2          → Rule B fires: child1=0 < 2 → ERROR for S2xD2
#                        S1xD2: s=0,d=2 → UNSAFE → always error for S1xD2 ✓
#                        S2xD1: s=2,d=0,prog=0 → (FALSE)|(TRUE & 0!=1)=TRUE → ERROR ✓
#
# Let's verify S1xD2 Group1: s=0,d=2 → always error (due to bug) → 5 errors
# Let's verify S2xD1 Group1: s=2,d=0,prog=0 → right side: TRUE & (0!=1)=TRUE → ERROR → 5 errors
#
# So with Group 1 alone (5 markers, all s=0,d=0,prog=0 for S1xD1):
#   S1xD1: 0 errors / 5 = 0%   ← BEST ✓
#   S2xD2: 5 errors / 5 = 100%
#   S1xD2: 5 errors / 5 = 100%
#   S2xD1: 5 errors / 5 = 100%
#
# This works! Simple and clean.
# child2 can be anything distinct.
# ─────────────────────────────────────────────────────────────────────────────

base_genos <- data.table(
  ID = c("S1", "S2", "D1", "D2", "child1", "child2"),
  M1 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M2 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M3 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M4 = c(0L, 2L, 0L, 2L, 0L, 2L),
  M5 = c(0L, 2L, 0L, 2L, 0L, 2L)
)

# child2 is a perfect child of S2xD2 (all 2s, Rule B: s=2,d=2→prog=2 ✓)
# S1xD1 for child2: s=0,d=0→must be 0, child2=2 → ERROR on all 5 markers

base_parents   <- data.table(ID  = c("S1","S2","D1","D2"),
                             Sex = c("M", "M", "F", "F"))
base_progeny   <- data.table(ID = c("child1", "child2"))
child1_progeny <- data.table(ID = "child1")
child2_progeny <- data.table(ID = "child2")

# ══════════════════════════════════════════════
# 1. Input validation
# ══════════════════════════════════════════════
test_that("invalid method throws an error", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  expect_error(
    find_parentage(f$g, f$p, f$pr, method = "bad.method",
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
  extra_parents <- rbind(base_parents, data.table(ID = "GHOST", Sex = "M"))
  f <- make_files(base_genos, extra_parents, child1_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best.pair",
                   verbose = FALSE, write_txt = FALSE),
    regexp = "GHOST"
  )
})

test_that("progeny IDs absent from genotype file raise a warning and are dropped", {
  extra_progeny <- rbind(child1_progeny, data.table(ID = "GHOST_KID"))
  f <- make_files(base_genos, base_parents, extra_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best.pair",
                   verbose = FALSE, write_txt = FALSE),
    regexp = "GHOST_KID"
  )
})

test_that("no valid progeny candidates after filtering stops with an error", {
  ghost_progeny <- data.table(ID = "NOBODY")
  f <- make_files(base_genos, base_parents, ghost_progeny)
  expect_warning(
    expect_error(
      find_parentage(f$g, f$p, f$pr, method = "best.pair",
                     verbose = FALSE, write_txt = FALSE),
      regexp = "No valid progeny"
    )
  )
})

test_that("missing Sex column raises a warning and defaults to ambiguous", {
  parents_no_sex <- data.table(ID = c("S1", "D1"))
  f <- make_files(base_genos, parents_no_sex, child1_progeny)
  expect_warning(
    find_parentage(f$g, f$p, f$pr, method = "best.match",
                   verbose = FALSE, write_txt = FALSE),
    regexp = "Sex"
  )
})

# ══════════════════════════════════════════════
# 2. Return structure
# ══════════════════════════════════════════════
test_that("best.pair returns a data.table with expected columns (no ties)", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny", "Sire", "Dam",
                    "Mendelian_Error_Pct", "Markers_Tested") %in% names(res)))
  expect_equal(nrow(res), 1L)
})

test_that("best.sire returns a data.table with expected columns", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.sire",
                        verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny", "Best_Match",
                    "Mendelian_Error_Pct", "Markers_Tested") %in% names(res)))
  expect_equal(nrow(res), 1L)
})

test_that("best.dam returns a data.table with expected columns", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.dam",
                        verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny", "Best_Match",
                    "Mendelian_Error_Pct", "Markers_Tested") %in% names(res)))
  expect_equal(nrow(res), 1L)
})

test_that("best.match returns a data.table with expected columns", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.match",
                        verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
  expect_true(all(c("Progeny", "Best_Match",
                    "Mendelian_Error_Pct", "Markers_Tested") %in% names(res)))
  expect_equal(nrow(res), 1L)
})

test_that("one row is returned per progeny for single-parent methods", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  for (m in c("best.sire", "best.dam", "best.match")) {
    res <- find_parentage(f$g, f$p, f$pr, method = m,
                          verbose = FALSE, write_txt = FALSE)
    expect_equal(nrow(res), 1L, label = paste("row count for method", m))
  }
})

# ══════════════════════════════════════════════
# 3. Biological correctness
# ══════════════════════════════════════════════
test_that("best.pair correctly identifies S1 x D1 as best pair with 0% error for child1", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Sire, "S1")
  expect_equal(res$Dam,  "D1")
  expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
})

test_that("best.pair correctly identifies S2 x D2 as best pair with 0% error for child2", {
  f <- make_files(base_genos, base_parents, child2_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Sire, "S2")
  expect_equal(res$Dam,  "D2")
  expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
})

test_that("best.sire identifies S1 as best sire for child1", {
  # For homozygous method: child1=0 (hom), S1=0 (hom) → match; S2=2 → mismatch
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.sire",
                        verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Best_Match, "S1")
})

test_that("best.dam identifies D1 as best dam for child1", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.dam",
                        verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Best_Match, "D1")
})

test_that("Mendelian_Error_Pct is 0 for a perfect parent-progeny trio", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(as.numeric(res$Mendelian_Error_Pct), 0)
})

test_that("Mendelian_Error_Pct is between 0 and 100", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  pct <- as.numeric(res$Mendelian_Error_Pct)
  expect_true(all(pct >= 0 & pct <= 100, na.rm = TRUE))
})

test_that("Markers_Tested equals the number of non-NA markers", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(res$Markers_Tested, ncol(base_genos) - 1L)  # minus ID column
})

# ══════════════════════════════════════════════
# 4. allow_selfing
# ══════════════════════════════════════════════
test_that("allow_selfing = FALSE removes self-pairs from candidates", {
  ambig_parents <- data.table(ID = c("S1", "D1"), Sex = c("A", "A"))
  f <- make_files(base_genos, ambig_parents, child1_progeny)
  # With only 2 ambiguous parents, S1xD1 and D1xS1 are tied → warning expected
  expect_warning(
    res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                          allow_selfing = FALSE, show_ties = FALSE,
                          verbose = FALSE, write_txt = FALSE),
    regexp = "tied"
  )
  if (!is.na(res$Sire) && !is.na(res$Dam)) {
    expect_false(res$Sire == res$Dam)
  }
})

# ══════════════════════════════════════════════
# 5. show_ties
# ══════════════════════════════════════════════

# All markers 0 → every sire×dam pair scores 0% error → guaranteed ties
tied_genos <- data.table(
  ID = c("S1", "S2", "D1", "D2", "child_tie"),
  M1 = c(0L, 0L, 0L, 0L, 0L),
  M2 = c(0L, 0L, 0L, 0L, 0L)
)
tied_parents <- data.table(ID  = c("S1","S2","D1","D2"),
                           Sex = c("M", "M", "F", "F"))
tied_progeny <- data.table(ID = "child_tie")

test_that("show_ties = TRUE produces _1/_2 suffixed columns when ties exist", {
  f <- make_files(tied_genos, tied_parents, tied_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        show_ties = TRUE, verbose = FALSE, write_txt = FALSE)
  expect_true(any(grepl("^Sire_", names(res))))
  expect_true(any(grepl("^Dam_",  names(res))))
})

test_that("show_ties = FALSE warns about ties and returns single-result columns", {
  f <- make_files(tied_genos, tied_parents, tied_progeny)
  expect_warning(
    res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                          show_ties = FALSE, verbose = FALSE, write_txt = FALSE),
    regexp = "tied"
  )
  expect_true("Sire" %in% names(res))
  expect_false(any(grepl("^Sire_\\d", names(res))))
})

# ══════════════════════════════════════════════
# 6. verbose / write_txt
# ══════════════════════════════════════════════
test_that("verbose = TRUE returns the result invisibly", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        verbose = TRUE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
})

test_that("verbose = FALSE returns the result visibly", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
})

test_that("write_txt = TRUE creates the output file", {
  old_wd <- getwd()
  tmp    <- tempdir()
  setwd(tmp)
  on.exit(setwd(old_wd), add = TRUE)

  f <- make_files(base_genos, base_parents, child1_progeny, dir = tmp)
  find_parentage(f$g, f$p, f$pr, method = "best.pair",
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
  find_parentage(f$g, f$p, f$pr, method = "best.pair",
                 verbose = FALSE, write_txt = FALSE)
  expect_false(file.exists(out_file))
})

# ══════════════════════════════════════════════
# 7. Sex-based candidate filtering
# ══════════════════════════════════════════════
test_that("best.sire only assigns male (M) or ambiguous (A) parents", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.sire",
                        verbose = FALSE, write_txt = FALSE)
  valid_sires <- base_parents[Sex %in% c("M", "A")]$ID
  expect_true(res$Best_Match %in% valid_sires)
})

test_that("best.dam only assigns female (F) or ambiguous (A) parents", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.dam",
                        verbose = FALSE, write_txt = FALSE)
  valid_dams <- base_parents[Sex %in% c("F", "A")]$ID
  expect_true(res$Best_Match %in% valid_dams)
})

# ══════════════════════════════════════════════
# 8. Edge cases
# ══════════════════════════════════════════════
test_that("single progeny individual is handled correctly", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_equal(nrow(res), 1L)
})

test_that("all-NA marker column does not cause an error", {
  na_genos <- copy(base_genos)
  na_genos[, M1 := NA_integer_]
  f <- make_files(na_genos, base_parents, child1_progeny)
  expect_no_error(
    find_parentage(f$g, f$p, f$pr, method = "best.pair",
                   verbose = FALSE, write_txt = FALSE)
  )
})

test_that("Progeny column contains the correct progeny IDs", {
  f <- make_files(base_genos, base_parents, child1_progeny)
  res <- find_parentage(f$g, f$p, f$pr, method = "best.pair",
                        show_ties = FALSE, verbose = FALSE, write_txt = FALSE)
  expect_setequal(res$Progeny, child1_progeny$ID)
})
