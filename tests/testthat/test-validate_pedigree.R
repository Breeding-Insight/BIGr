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
      rep(0L, n_markers),   # IND_A — all ref homozygous
      rep(2L, n_markers),   # IND_B — all alt homozygous
      rep(1L, n_markers),   # IND_C — all het (valid child of IND_A x IND_B)
      rep(0L, n_markers),   # IND_D — impossible child of IND_B x IND_A
      rep(0L, n_markers)    # IND_E — all ref
    )
  )
  setnames(dt, c("id", marker_names))
  dt
}

make_pedigree <- function() {
  # IND_C: perfect Mendelian child of IND_A x IND_B -> pass
  # IND_D: declared parents swapped -> fail
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
  expect_error(validate_pedigree(f$ped, f$genos,
                                 trio_error_threshold = 150,
                                 verbose = FALSE, plot_results = FALSE))
  expect_error(validate_pedigree(f$ped, f$genos,
                                 trio_error_threshold = -1,
                                 verbose = FALSE, plot_results = FALSE))
})

test_that("single_parent_error_threshold out of range raises an error", {
  f <- write_temp_files()
  expect_error(validate_pedigree(f$ped, f$genos,
                                 single_parent_error_threshold = 101,
                                 verbose = FALSE, plot_results = FALSE))
  expect_error(validate_pedigree(f$ped, f$genos,
                                 single_parent_error_threshold = -5,
                                 verbose = FALSE, plot_results = FALSE))
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

test_that("full_results is a data.table", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  expect_s3_class(out$full_results, "data.table")
})

test_that("corrected_pedigree is a data.table with lowercase columns", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  expect_s3_class(out$corrected_pedigree, "data.table")
  expect_true(all(c("id", "male_parent", "female_parent") %in%
                    names(out$corrected_pedigree)))
})

test_that("plot element is NULL when plot_results = FALSE", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  expect_null(out$plot)
})

test_that("named list subsets sum correctly to full_results row count", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  subset_total <- nrow(out$pass) + nrow(out$fail) + nrow(out$low_markers) +
    nrow(out$no_genotype_data) + nrow(out$founders) + nrow(out$missing_parents)
  expect_equal(subset_total, nrow(out$full_results))
})

# ==============================================================================
# 3. pass / fail / low_markers / no_data statuses — lowercase
# ==============================================================================
test_that("pass trio is correctly identified", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  pass_row <- out$full_results[id == "IND_C"]
  expect_equal(nrow(pass_row), 1L)
  expect_equal(pass_row$status, "pass")
  expect_equal(pass_row$trio_mendelian_error_pct, 0)
  expect_equal(pass_row$recommended_correction, "none")
})

test_that("fail trio is correctly identified", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  fail_row <- out$full_results[id == "IND_D"]
  expect_equal(nrow(fail_row), 1L)
  expect_equal(fail_row$status, "fail")
  expect_gt(fail_row$trio_mendelian_error_pct, 5.0)
})

test_that("fail trio has remove_male_parent decision and best candidate populated", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  fail_row <- out$full_results[id == "IND_D"]
  expect_equal(fail_row$recommended_correction, "remove_male_parent")
  expect_false(is.na(fail_row$best_male_candidate))
  expect_true(is.na(fail_row$best_female_candidate))
})

test_that("trio_mendelian_error_pct is 0 for a perfect trio", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results[id == "IND_C"]$trio_mendelian_error_pct, 0)
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

  expect_true(all(out$pass$status == "pass"))
})

test_that("fail list element contains only fail rows", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  expect_true(all(out$fail$status == "fail"))
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
               data.table(id = "IND_E", male_parent = "0", female_parent = "0"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_E"]

  expect_equal(r$status,                 "missing_both_parents")
  expect_equal(r$recommended_correction, "none")
  expect_false(is.na(r$best_male_candidate))
  expect_false(is.na(r$best_female_candidate))
})

test_that("best_male_candidate for missing_male_parent is not the known female parent", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "0",
                          female_parent = "IND_B"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  r   <- out$full_results[id == "IND_E"]

  expect_false(r$best_male_candidate == "IND_B")
})

test_that("missing_parents list element contains only missing_* rows", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "0", female_parent = "0"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  expect_true(all(grepl("^missing_", out$missing_parents$status)))
})

# ==============================================================================
# 5. founders status
# ==============================================================================
test_that("founders status is assigned when ID in founders list with 0 0 parents", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_A", male_parent = "0", female_parent = "0"))
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
               data.table(id = "IND_A", male_parent = "0", female_parent = "0"))
  founders_file <- tempfile(fileext = ".txt")
  fwrite(data.table(id = "IND_A"), founders_file,
         sep = "\t", quote = FALSE, col.names = FALSE)
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, founders_file = founders_file,
                           verbose = FALSE, plot_results = FALSE)

  expect_true(all(out$founders$status == "founders"))
})

test_that("non-founder rows still evaluated normally when founders file is supplied", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_A", male_parent = "0", female_parent = "0"))
  founders_file <- tempfile(fileext = ".txt")
  fwrite(data.table(id = "IND_A"), founders_file,
         sep = "\t", quote = FALSE, col.names = FALSE)
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, founders_file = founders_file,
                           verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results[id == "IND_C"]$status, "pass")
})

test_that("0 0 parents NOT in founders list get missing_both_parents", {
  ped <- rbind(make_pedigree(),
               data.table(id = "IND_E", male_parent = "0", female_parent = "0"))
  f   <- write_temp_files(ped = ped)
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  expect_equal(out$full_results[id == "IND_E"]$status, "missing_both_parents")
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

test_that("no_genotype_data flagged when declared parent absent from genotype file", {
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

test_that("corrected_pedigree: bad parent set to 0 for fail trio with remove_male_parent", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)
  corr <- out$corrected_pedigree

  expect_equal(corr[id == "IND_D"]$male_parent,   "0")
  expect_equal(corr[id == "IND_D"]$female_parent, "IND_A")
})

test_that("corrected_pedigree has same number of rows as original pedigree", {
  f   <- write_temp_files()
  out <- validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = FALSE)

  expect_equal(nrow(out$corrected_pedigree), nrow(make_pedigree()))
})

# ==============================================================================
# 8. No write logic — functions do not write files
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
# 9. verbose and plot_results
# ==============================================================================
test_that("verbose = FALSE suppresses console output", {
  f <- write_temp_files()
  expect_silent(validate_pedigree(f$ped, f$genos, verbose = FALSE,
                                  plot_results = FALSE))
})

test_that("plot element is a ggplot when plot_results = TRUE", {
  skip_if_not_installed("ggplot2")
  f   <- write_temp_files()
  out <- suppressWarnings(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, plot_results = TRUE)
  )
  expect_s3_class(out$plot, "ggplot")
})

# ==============================================================================
# 10. Return value is invisible
# ==============================================================================
test_that("validate_pedigree returns invisibly", {
  f <- write_temp_files()
  expect_invisible(validate_pedigree(f$ped, f$genos, verbose = FALSE,
                                     plot_results = FALSE))
})
