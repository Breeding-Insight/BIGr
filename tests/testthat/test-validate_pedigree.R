# tests/testthat/test-validate_pedigree.R
# Run with: testthat::test_file("tests/testthat/test-validate_pedigree.R")

library(testthat)
library(data.table)

# ==============================================================================
# Helpers
# ==============================================================================

make_genos <- function() {
  n_markers    <- 20
  marker_names <- paste0("M", seq_len(n_markers))
  dt <- data.table(
    id = c("IND_A", "IND_B", "IND_C", "IND_D", "IND_E"),
    rbind(
      rep(0L, n_markers),  # IND_A — all ref homozygous
      rep(2L, n_markers),  # IND_B — all alt homozygous
      rep(1L, n_markers),  # IND_C — all het: valid child of IND_A x IND_B
      rep(0L, n_markers),  # IND_D — impossible child of IND_B x IND_A
      rep(0L, n_markers)   # IND_E — all ref
    )
  )
  setnames(dt, c("id", marker_names))
  dt
}

make_pedigree <- function() {
  # IND_C: perfect Mendelian child of IND_A x IND_B -> pass
  # IND_D: declared parents swapped    -> fail
  data.table(
    id            = c("IND_C", "IND_D"),
    male_parent   = c("IND_A", "IND_B"),
    female_parent = c("IND_B", "IND_A")
  )
}

write_temp_files <- function(genos = make_genos(), ped = make_pedigree()) {
  ped_file   <- tempfile(fileext = ".txt")
  genos_file <- tempfile(fileext = ".txt")
  fwrite(ped,   ped_file,   sep = "\t")
  fwrite(genos, genos_file, sep = "\t")
  list(ped = ped_file, genos = genos_file)
}

# ==============================================================================
# 1. Input validation
# ==============================================================================

test_that("trio_error_threshold out of range raises an error", {
  f <- write_temp_files()
  expect_error(
    validate_pedigree(f$ped, f$genos, trio_error_threshold = 150,
                      verbose = FALSE, plot_results = FALSE),
    regexp = "trio_error_threshold"
  )
  expect_error(
    validate_pedigree(f$ped, f$genos, trio_error_threshold = -1,
                      verbose = FALSE, plot_results = FALSE),
    regexp = "trio_error_threshold"
  )
})

test_that("single_parent_error_threshold out of range raises an error", {
  f <- write_temp_files()
  expect_error(
    validate_pedigree(f$ped, f$genos, single_parent_error_threshold = 101,
                      verbose = FALSE, plot_results = FALSE),
    regexp = "single_parent_error_threshold"
  )
  expect_error(
    validate_pedigree(f$ped, f$genos, single_parent_error_threshold = -5,
                      verbose = FALSE, plot_results = FALSE),
    regexp = "single_parent_error_threshold"
  )
})

test_that("boundary values 0 and 100 are accepted for trio_error_threshold", {
  f <- write_temp_files()
  expect_no_error(
    validate_pedigree(f$ped, f$genos, trio_error_threshold = 0,
                      verbose = FALSE, plot_results = FALSE)
  )
  expect_no_error(
    validate_pedigree(f$ped, f$genos, trio_error_threshold = 100,
                      verbose = FALSE, plot_results = FALSE)
  )
})

test_that("nonexistent pedigree file throws 'Error reading input files'", {
  f <- write_temp_files()
  expect_error(
    validate_pedigree("nonexistent.txt", f$genos,
                      verbose = FALSE, plot_results = FALSE),
    regexp = "Error reading input files"
  )
})

test_that("nonexistent genotypes file throws 'Error reading input files'", {
  f <- write_temp_files()
  expect_error(
    validate_pedigree(f$ped, "nonexistent.txt",
                      verbose = FALSE, plot_results = FALSE),
    regexp = "Error reading input files"
  )
})

test_that("missing required pedigree column raises an error", {
  bad_ped <- data.table(id = "IND_C", parent1 = "IND_A", female_parent = "IND_B")
  f <- write_temp_files(ped = bad_ped)
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE),
    regexp = "missing required columns"
  )
})

test_that("missing id column in genotypes raises an error", {
  bad_genos <- copy(make_genos())
  setnames(bad_genos, "id", "SampleID")
  f <- write_temp_files(genos = bad_genos)
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE),
    regexp = "id"
  )
})

test_that("all trios with no genotype data stops with an error", {
  ped <- data.table(id = "GHOST", male_parent = "IND_A", female_parent = "IND_B")
  f   <- write_temp_files(ped = ped)
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE),
    regexp = "No valid trios remain"
  )
})

test_that("unreadable founders file raises an error", {
  f <- write_temp_files()
  expect_error(
    validate_pedigree(f$ped, f$genos, founders_file = "nonexistent_founders.txt",
                      verbose = FALSE, plot_results = FALSE)
  )
})

# ==============================================================================
# 2. Return structure — named list
# ==============================================================================

test_that("returns an invisible named list with all required elements", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_type(out, "list")
  expect_named(out, c("pass", "fail", "low_markers", "no_genotype_data",
                      "founders", "missing_parents", "full_results",
                      "corrected_pedigree", "plot"),
               ignore.order = TRUE)
})

test_that("validate_pedigree returns invisibly", {
  f <- write_temp_files()
  expect_invisible(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  )
})

test_that("full_results is a data.table", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_s3_class(out$full_results, "data.table")
})

test_that("full_results has one row per pedigree entry", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_equal(nrow(out$full_results), 2L)
})

test_that("full_results has all expected lowercase columns", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expected_cols <- c(
    "id", "orig_male_parent", "orig_female_parent",
    "trio_mendelian_error_pct", "trio_markers_tested", "status",
    "recommended_correction",
    "male_parent_hom_error_pct", "female_parent_hom_error_pct",
    "best_male_candidate",   "best_male_candidate_error_pct",
    "best_female_candidate", "best_female_candidate_error_pct"
  )
  expect_true(all(expected_cols %in% names(out$full_results)))
})

test_that("corrected_pedigree is a data.table with lowercase columns", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_s3_class(out$corrected_pedigree, "data.table")
  expect_true(all(c("id", "male_parent", "female_parent") %in%
                    names(out$corrected_pedigree)))
})

test_that("corrected_pedigree has same number of rows as original pedigree", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_equal(nrow(out$corrected_pedigree), nrow(make_pedigree()))
})

test_that("plot element is NULL when plot_results = FALSE", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_null(out$plot)
})

test_that("plot element is a ggplot when plot_results = TRUE", {
  skip_if_not_installed("ggplot2")
  f   <- write_temp_files()
  out <- suppressWarnings(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = TRUE)
  )
  expect_s3_class(out$plot, "ggplot")
})

test_that("named list subsets sum correctly to full_results row count", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  subset_total <- nrow(out$pass) + nrow(out$fail) + nrow(out$low_markers) +
    nrow(out$no_genotype_data) + nrow(out$founders) + nrow(out$missing_parents)
  expect_equal(subset_total, nrow(out$full_results))
})

# ==============================================================================
# 3. pass / fail / low_markers statuses
# ==============================================================================

test_that("pass trio is correctly identified with 0% error", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_C"]
  expect_equal(nrow(r), 1L)
  expect_equal(r$status, "pass")
  expect_equal(r$trio_mendelian_error_pct, 0)
  expect_equal(r$recommended_correction, "none")
})

test_that("fail trio is correctly identified with error above threshold", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_D"]
  expect_equal(nrow(r), 1L)
  expect_equal(r$status, "fail")
  expect_gt(r$trio_mendelian_error_pct, 5.0)
})

test_that("fail trio has a non-NA recommended_correction", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_D"]
  expect_false(is.na(r$recommended_correction))
  expect_false(r$recommended_correction == "none")
})

test_that("fail trio with one acceptable parent gets remove_* correction", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_D"]
  expect_true(r$recommended_correction %in%
                c("remove_male_parent", "remove_female_parent", "remove_both",
                  "keep_both"))
})

test_that("trio_mendelian_error_pct is 0 for a perfect Mendelian trio", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results[id == "IND_C"]$trio_mendelian_error_pct, 0)
})

test_that("trio_mendelian_error_pct is between 0 and 100 for all trios", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  pct <- out$full_results$trio_mendelian_error_pct
  expect_true(all(pct >= 0 & pct <= 100, na.rm = TRUE))
})

test_that("trio_markers_tested equals number of markers for complete data", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results[id == "IND_C"]$trio_markers_tested, 20L)
})

test_that("low_markers status assigned when markers_tested < min_markers", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE,
                           plot_results = FALSE, min_markers = 25L)
  expect_true(all(out$full_results$status == "low_markers"))
  expect_true(all(grepl("^low_markers_", out$full_results$recommended_correction)))
})

test_that("NA markers reduce trio_markers_tested and do not cause errors", {
  genos <- make_genos()
  genos[id == "IND_C", M1 := NA_integer_]
  genos[id == "IND_C", M2 := NA_integer_]
  f   <- write_temp_files(genos = genos)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results[id == "IND_C"]$trio_markers_tested, 18L)
  expect_equal(out$full_results[id == "IND_C"]$status, "pass")
})

test_that("pass list element contains only pass rows", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  if (nrow(out$pass) > 0)
    expect_true(all(out$pass$status == "pass"))
})

test_that("fail list element contains only fail rows", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  if (nrow(out$fail) > 0)
    expect_true(all(out$fail$status == "fail"))
})

test_that("low_markers list element contains only low_markers rows", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE,
                           plot_results = FALSE, min_markers = 25L)
  if (nrow(out$low_markers) > 0)
    expect_true(all(out$low_markers$status == "low_markers"))
})

test_that("raising trio_error_threshold turns fail rows into pass rows", {
  f      <- write_temp_files()
  strict <- validate_pedigree(f$ped, f$genos, trio_error_threshold = 5.0,
                              verbose = FALSE, plot_results = FALSE)
  lenient <- validate_pedigree(f$ped, f$genos, trio_error_threshold = 100.0,
                               verbose = FALSE, plot_results = FALSE)
  expect_gte(nrow(lenient$pass), nrow(strict$pass))
})

# ==============================================================================
# 4. missing parent statuses
# ==============================================================================

test_that("missing_male_parent status and recommendation are correct", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "0",
                          female_parent = "IND_B"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_E"]
  expect_equal(r$status,                 "missing_male_parent")
  expect_equal(r$recommended_correction, "none")
  expect_false(is.na(r$best_male_candidate))
  expect_true(is.na(r$best_female_candidate))
})

test_that("missing_female_parent status and recommendation are correct", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "IND_A",
                          female_parent = "0"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_E"]
  expect_equal(r$status,                 "missing_female_parent")
  expect_equal(r$recommended_correction, "none")
  expect_true(is.na(r$best_male_candidate))
  expect_false(is.na(r$best_female_candidate))
})

test_that("missing_both_parents status and recommendations are correct", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "0",
                          female_parent = "0"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_E"]
  expect_equal(r$status,                 "missing_both_parents")
  expect_equal(r$recommended_correction, "none")
  expect_false(is.na(r$best_male_candidate))
  expect_false(is.na(r$best_female_candidate))
})

test_that("best_male_candidate for missing_male_parent excludes the known female parent", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "0",
                          female_parent = "IND_B"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_E"]
  expect_false(r$best_male_candidate == "IND_B")
})

test_that("best_female_candidate for missing_female_parent excludes the known male parent", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "IND_A",
                          female_parent = "0"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_E"]
  expect_false(r$best_female_candidate == "IND_A")
})

test_that("missing_parents list element contains only missing_* rows", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "0",
                          female_parent = "0"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  if (nrow(out$missing_parents) > 0)
    expect_true(all(grepl("^missing_", out$missing_parents$status)))
})

# ==============================================================================
# 5. founders status
# ==============================================================================

test_that("founders status is assigned when ID is in founders list with 0 0 parents", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_A", male_parent = "0",
                          female_parent = "0"))
  founders_file <- tempfile(fileext = ".txt")
  fwrite(data.table(id = "IND_A"), founders_file,
         sep = "\t", quote = FALSE, col.names = FALSE)
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, founders_file = founders_file,
                           verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_A"]
  expect_equal(r$status,                 "founders")
  expect_equal(r$recommended_correction, "none")
  expect_true(is.na(r$best_male_candidate))
  expect_true(is.na(r$best_female_candidate))
})

test_that("founders list element contains only founders rows", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_A", male_parent = "0",
                          female_parent = "0"))
  founders_file <- tempfile(fileext = ".txt")
  fwrite(data.table(id = "IND_A"), founders_file,
         sep = "\t", quote = FALSE, col.names = FALSE)
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, founders_file = founders_file,
                           verbose = FALSE, plot_results = FALSE)
  if (nrow(out$founders) > 0)
    expect_true(all(out$founders$status == "founders"))
})

test_that("non-founder rows still evaluated normally when founders file is supplied", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_A", male_parent = "0",
                          female_parent = "0"))
  founders_file <- tempfile(fileext = ".txt")
  fwrite(data.table(id = "IND_A"), founders_file,
         sep = "\t", quote = FALSE, col.names = FALSE)
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, founders_file = founders_file,
                           verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results[id == "IND_C"]$status, "pass")
  expect_equal(out$full_results[id == "IND_D"]$status, "fail")
})

test_that("0 0 parents NOT in founders list get missing_both_parents", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "0",
                          female_parent = "0"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results[id == "IND_E"]$status, "missing_both_parents")
})

test_that("founder row does not appear in pass, fail, or missing_parents", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_A", male_parent = "0",
                          female_parent = "0"))
  founders_file <- tempfile(fileext = ".txt")
  fwrite(data.table(id = "IND_A"), founders_file,
         sep = "\t", quote = FALSE, col.names = FALSE)
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, founders_file = founders_file,
                           verbose = FALSE, plot_results = FALSE)
  expect_false("IND_A" %in% out$pass$id)
  expect_false("IND_A" %in% out$fail$id)
  expect_false("IND_A" %in% out$missing_parents$id)
})

# ==============================================================================
# 6. no_genotype_data status
# ==============================================================================

test_that("no_genotype_data is flagged for progeny absent from genotype file", {
  ped <- rbind(make_pedigree(),
               data.table(id = "GHOST", male_parent = "IND_A",
                          female_parent = "IND_B"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "GHOST"]
  expect_equal(nrow(r),                  1L)
  expect_equal(r$status,                 "no_genotype_data")
  expect_equal(r$recommended_correction, "none")
})

test_that("no_genotype_data rows have NA/0 for all analysis columns", {
  ped <- rbind(make_pedigree(),
               data.table(id = "GHOST", male_parent = "IND_A",
                          female_parent = "IND_B"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "GHOST"]
  expect_true(is.na(r$trio_mendelian_error_pct))
  expect_equal(r$trio_markers_tested, 0L)
  expect_true(is.na(r$best_male_candidate))
  expect_true(is.na(r$best_female_candidate))
})

test_that("no_genotype_data flagged when a declared parent is absent from genotype file", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_C_GHOST", male_parent = "GHOST_DAD",
                          female_parent = "IND_B"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results[id == "IND_C_GHOST"]$status, "no_genotype_data")
})

test_that("no_genotype_data list element contains only no_genotype_data rows", {
  ped <- rbind(make_pedigree(),
               data.table(id = "GHOST", male_parent = "IND_A",
                          female_parent = "IND_B"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  if (nrow(out$no_genotype_data) > 0)
    expect_true(all(out$no_genotype_data$status == "no_genotype_data"))
})

test_that("valid trios still evaluated correctly when ghost rows are present", {
  ped <- rbind(make_pedigree(),
               data.table(id = "GHOST", male_parent = "IND_A",
                          female_parent = "IND_B"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results[id == "IND_C"]$status, "pass")
  expect_equal(out$full_results[id == "IND_D"]$status, "fail")
})

# ==============================================================================
# 7. corrected_pedigree contents
# ==============================================================================

test_that("corrected_pedigree: pass parents are unchanged", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  corr <- out$corrected_pedigree
  expect_equal(as.character(corr[id == "IND_C"]$male_parent),   "IND_A")
  expect_equal(as.character(corr[id == "IND_C"]$female_parent), "IND_B")
})

test_that("corrected_pedigree: removed parent set to 0 for remove_male_parent", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_D"]
  corr <- out$corrected_pedigree
  if (r$recommended_correction == "remove_male_parent") {
    expect_equal(corr[id == "IND_D"]$male_parent, "0")
  }
})

test_that("corrected_pedigree: removed parent set to 0 for remove_female_parent", {
  # construct a trio where IND_A (all 0) is correct male and female is wrong
  genos <- make_genos()
  ped   <- data.table(id = "IND_E", male_parent = "IND_A",
                      female_parent = "IND_B")
  # IND_E is all ref (0); IND_A is all ref (0); IND_B is all alt (2)
  # IND_E as child of IND_A x IND_B is impossible → remove_female_parent
  f   <- write_temp_files(genos = genos, ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  corr <- out$corrected_pedigree
  r    <- out$full_results[id == "IND_E"]
  if (r$recommended_correction == "remove_female_parent") {
    expect_equal(corr[id == "IND_E"]$female_parent, "0")
    expect_false(corr[id == "IND_E"]$male_parent == "0")
  }
})

test_that("corrected_pedigree preserves id column values", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_setequal(out$corrected_pedigree$id, make_pedigree()$id)
})

# ==============================================================================
# 8. No write logic — function does not write files
# ==============================================================================

test_that("no output files are written to disk", {
  f       <- write_temp_files()
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)
  validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  written_files <- list.files(tmp_dir)
  expect_length(written_files, 0)
})

# ==============================================================================
# 9. verbose
# ==============================================================================

test_that("verbose = FALSE suppresses console output", {
  f <- write_temp_files()
  expect_silent(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  )
})

test_that("verbose = TRUE returns valid named list without error", {
  f <- write_temp_files()
  invisible(capture.output(
    out <- validate_pedigree(f$ped, f$genos, verbose = TRUE, plot_results = FALSE)
  ))
  expect_type(out, "list")
  expect_named(out, c("pass", "fail", "low_markers", "no_genotype_data",
                      "founders", "missing_parents", "full_results",
                      "corrected_pedigree", "plot"),
               ignore.order = TRUE)
})

# ==============================================================================
# 10. Mendelian error correctness
# ==============================================================================

test_that("0x0 parents produce 0% error for dosage-0 child", {
  genos <- data.table(
    id = c("S", "D", "C"),
    M1 = c(0L, 0L, 0L), M2 = c(0L, 0L, 0L), M3 = c(0L, 0L, 0L),
    M4 = c(0L, 0L, 0L), M5 = c(0L, 0L, 0L)
  )
  ped <- data.table(id = "C", male_parent = "S", female_parent = "D")
  f   <- write_temp_files(genos = genos, ped = ped)
  out <- validate_pedigree(f$ped, f$genos, trio_error_threshold = 5.0,
                           min_markers = 1, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results$trio_mendelian_error_pct, 0)
})

test_that("2x2 parents produce 0% error for dosage-2 child", {
  genos <- data.table(
    id = c("S", "D", "C"),
    M1 = c(2L, 2L, 2L), M2 = c(2L, 2L, 2L), M3 = c(2L, 2L, 2L),
    M4 = c(2L, 2L, 2L), M5 = c(2L, 2L, 2L)
  )
  ped <- data.table(id = "C", male_parent = "S", female_parent = "D")
  f   <- write_temp_files(genos = genos, ped = ped)
  out <- validate_pedigree(f$ped, f$genos, trio_error_threshold = 5.0,
                           min_markers = 1, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results$trio_mendelian_error_pct, 0)
})

test_that("0x0 parents produce 100% error for dosage-2 child", {
  genos <- data.table(
    id = c("S", "D", "C"),
    M1 = c(0L, 0L, 2L), M2 = c(0L, 0L, 2L), M3 = c(0L, 0L, 2L),
    M4 = c(0L, 0L, 2L), M5 = c(0L, 0L, 2L)
  )
  ped <- data.table(id = "C", male_parent = "S", female_parent = "D")
  f   <- write_temp_files(genos = genos, ped = ped)
  out <- validate_pedigree(f$ped, f$genos, trio_error_threshold = 5.0,
                           min_markers = 1, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results$trio_mendelian_error_pct, 100)
})

test_that("0x2 parents produce 0% error for dosage-1 child", {
  genos <- data.table(
    id = c("S", "D", "C"),
    M1 = c(0L, 2L, 1L), M2 = c(0L, 2L, 1L), M3 = c(0L, 2L, 1L),
    M4 = c(0L, 2L, 1L), M5 = c(0L, 2L, 1L)
  )
  ped <- data.table(id = "C", male_parent = "S", female_parent = "D")
  f   <- write_temp_files(genos = genos, ped = ped)
  out <- validate_pedigree(f$ped, f$genos, trio_error_threshold = 5.0,
                           min_markers = 1, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results$trio_mendelian_error_pct, 0)
})

test_that("0x2 parents produce 100% error for dosage-0 child", {
  genos <- data.table(
    id = c("S", "D", "C"),
    M1 = c(0L, 2L, 0L), M2 = c(0L, 2L, 0L), M3 = c(0L, 2L, 0L),
    M4 = c(0L, 2L, 0L), M5 = c(0L, 2L, 0L)
  )
  ped <- data.table(id = "C", male_parent = "S", female_parent = "D")
  f   <- write_temp_files(genos = genos, ped = ped)
  out <- validate_pedigree(f$ped, f$genos, trio_error_threshold = 5.0,
                           min_markers = 1, verbose = FALSE, plot_results = FALSE)
  expect_equal(out$full_results$trio_mendelian_error_pct, 100)
})
# ==============================================================================
# 11. In-memory input — data.frame / data.table accepted directly
# ==============================================================================

test_that("validate_pedigree accepts a data.table pedigree directly", {
  f   <- write_temp_files()
  ped <- make_pedigree()
  out <- validate_pedigree(ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_s3_class(out$full_results, "data.table")
  expect_equal(nrow(out$full_results), 2L)
})

test_that("validate_pedigree accepts a data.table genotypes object directly", {
  f     <- write_temp_files()
  genos <- make_genos()
  out   <- validate_pedigree(f$ped, genos, verbose = FALSE, plot_results = FALSE)
  expect_s3_class(out$full_results, "data.table")
  expect_equal(nrow(out$full_results), 2L)
})

test_that("validate_pedigree accepts both inputs as data.tables directly", {
  ped   <- make_pedigree()
  genos <- make_genos()
  out   <- validate_pedigree(ped, genos, verbose = FALSE, plot_results = FALSE)
  expect_s3_class(out$full_results, "data.table")
  expect_equal(nrow(out$full_results), 2L)
})

test_that("validate_pedigree accepts a data.frame pedigree directly", {
  f   <- write_temp_files()
  ped <- as.data.frame(make_pedigree())
  out <- validate_pedigree(ped, f$genos, verbose = FALSE, plot_results = FALSE)
  expect_s3_class(out$full_results, "data.table")
  expect_equal(nrow(out$full_results), 2L)
})

test_that("in-memory and file-path inputs produce identical results for validate_pedigree", {
  f     <- write_temp_files()
  ped   <- make_pedigree()
  genos <- make_genos()
  out_file   <- validate_pedigree(f$ped, f$genos,
                                  verbose = FALSE, plot_results = FALSE)
  out_memory <- validate_pedigree(ped, genos,
                                  verbose = FALSE, plot_results = FALSE)
  expect_equal(out_file$full_results$status,
               out_memory$full_results$status)
  expect_equal(out_file$full_results$trio_mendelian_error_pct,
               out_memory$full_results$trio_mendelian_error_pct)
})

test_that("invalid input type raises an error for validate_pedigree", {
  f <- write_temp_files()
  expect_error(
    validate_pedigree(list(id = "IND_C"), f$genos,
                      verbose = FALSE, plot_results = FALSE),
    regexp = "Error reading input files"
  )
})
