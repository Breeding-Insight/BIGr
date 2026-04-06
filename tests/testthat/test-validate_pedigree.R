#### Tests for validate_pedigree() ####
# Run with: testthat::test_file("test-validate_pedigree.R")
# Requires:  data.table, testthat

library(testthat)
library(data.table)

#### Helpers: Minimal test data ####

# Genotypes: 20 markers, coded 0/1/2
# IND_A and IND_B are parents; IND_C is a valid progeny (Mendelian-consistent)
# IND_D is a progeny whose sire is wrong (high Mendelian error)

make_genos <- function() {
  set.seed(42)
  n_markers <- 20
  marker_names <- paste0("M", seq_len(n_markers))

  # Parent A: all homozygous ref (0)
  pa <- rep(0L, n_markers)
  # Parent B: all homozygous alt (2)
  pb <- rep(2L, n_markers)
  # Valid progeny: all het (1) — perfectly Mendelian from A x B
  pc <- rep(1L, n_markers)
  # Bad progeny: all homozygous ref (0) — impossible if sire is B (2) and dam is A (0)
  # 0/0 x 2/2 -> must be 1; so all-0 is 100% error
  pd <- rep(0L, n_markers)
  # Parent C2: homozygous ref (0) — correct sire for pd
  pe <- rep(0L, n_markers)

  dt <- data.table(
    ID = c("IND_A", "IND_B", "IND_C", "IND_D", "IND_E"),
    rbind(pa, pb, pc, pd, pe)
  )
  setnames(dt, c("ID", marker_names))
  dt
}

make_pedigree <- function() {
  data.table(
    Progeny = c("IND_C", "IND_D"),
    Sire    = c("IND_A", "IND_B"),   # IND_D sire is wrong (IND_B)
    Dam     = c("IND_B", "IND_A")
  )
}

# Write temp files and return paths
write_temp_files <- function(genos = make_genos(), ped = make_pedigree()) {
  ped_file   <- tempfile(fileext = ".txt")
  genos_file <- tempfile(fileext = ".txt")
  fwrite(ped,   ped_file,   sep = "\t")
  fwrite(genos, genos_file, sep = "\t")
  list(ped = ped_file, genos = genos_file)
}

#### Test suite ####

test_that("PASS trio is correctly identified", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  pass_row <- res[Progeny == "IND_C"]
  expect_equal(nrow(pass_row), 1L)
  expect_equal(pass_row$Status, "PASS")
  expect_equal(pass_row$Mendelian_Error_Pct, 0)
  expect_equal(pass_row$Correction_Decision, "NONE")
})

test_that("FAIL trio is correctly identified", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  fail_row <- res[Progeny == "IND_D"]
  expect_equal(nrow(fail_row), 1L)
  expect_equal(fail_row$Status, "FAIL")
  expect_gt(fail_row$Mendelian_Error_Pct, 5.0)
})

test_that("FAIL trio has correct correction decision (REMOVE_SIRE)", {
  # IND_D: sire IND_B (all-2) is wrong; dam IND_A (all-0) matches IND_D (all-0) perfectly
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  fail_row <- res[Progeny == "IND_D"]
  expect_equal(fail_row$Correction_Decision, "REMOVE_SIRE")
  expect_false(is.na(fail_row$Best_Sire))
  expect_true(is.na(fail_row$Best_Dam))          # dam was fine, no replacement needed
})

test_that("Mendelian_Error_Pct is 0 for perfect trio", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  expect_equal(res[Progeny == "IND_C"]$Mendelian_Error_Pct, 0)
})

test_that("Markers_Tested equals number of markers for complete data", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  expect_equal(res[Progeny == "IND_C"]$Markers_Tested, 20L)
})

test_that("Returns a data.table invisibly", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), 2L)
})

test_that("Result has all expected columns", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  expected_cols <- c("Progeny", "Sire", "Dam", "Mendelian_Error_Pct",
                     "Markers_Tested", "Status", "Correction_Decision",
                     "Sire_Hom_Error_Pct", "Dam_Hom_Error_Pct",
                     "Best_Sire", "Best_Sire_Error_Pct",
                     "Best_Dam",  "Best_Dam_Error_Pct")
  expect_true(all(expected_cols %in% names(res)))
})

test_that("write_txt writes output file with correct name", {
  f        <- write_temp_files()
  out_file <- tempfile(fileext = ".txt")
  validate_pedigree(f$ped, f$genos, verbose = FALSE,
                    write_txt = TRUE, output_filename = out_file)

  expect_true(file.exists(out_file))
  written <- fread(out_file)
  expect_equal(nrow(written), 2L)
})

test_that("write_txt = FALSE does not create default output file", {
  f        <- write_temp_files()
  out_file <- tempfile(fileext = ".txt")
  validate_pedigree(f$ped, f$genos, verbose = FALSE,
                    write_txt = FALSE, output_filename = out_file)

  expect_false(file.exists(out_file))
})

test_that("corrected_pedigree.txt is always written with zeros for bad parents", {
  f       <- write_temp_files()
  tmp_dir <- tempfile()                  # unique dir per test — no cross-test pollution
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

  validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  corr <- fread(file.path(tmp_dir, "corrected_pedigree.txt"))
  # IND_D sire should be zeroed out (written as character "0")
  expect_equal(corr[Progeny == "IND_D"]$Sire, "0")
  # IND_D dam should be unchanged
  expect_equal(corr[Progeny == "IND_D"]$Dam, "IND_A")
  # IND_C should be completely unchanged
  expect_equal(corr[Progeny == "IND_C"]$Sire, "IND_A")
  expect_equal(corr[Progeny == "IND_C"]$Dam,  "IND_B")
})

test_that("fill_pedigree = TRUE writes filled_pedigree.txt with replacement IDs", {
  f       <- write_temp_files()
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

  validate_pedigree(f$ped, f$genos, verbose = FALSE,
                    write_txt = FALSE, fill_pedigree = TRUE)

  filled <- fread(file.path(tmp_dir, "filled_pedigree.txt"))
  # IND_D sire should be replaced with a valid ID (not 0, not the wrong IND_B)
  new_sire <- filled[Progeny == "IND_D"]$Sire
  expect_false(new_sire == "IND_B")
  expect_false(new_sire == "0")
})

test_that("fill_pedigree = FALSE does not write filled_pedigree.txt", {
  f       <- write_temp_files()
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

  validate_pedigree(f$ped, f$genos, verbose = FALSE,
                    write_txt = FALSE, fill_pedigree = FALSE)

  expect_false(file.exists(file.path(tmp_dir, "filled_pedigree.txt")))
})

test_that("Trios with missing genotype data are removed with a message", {
  ped <- data.table(Progeny = "GHOST", Sire = "IND_A", Dam = "IND_B")
  f   <- write_temp_files(ped = ped)

  # No valid trios remain -> should stop
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE),
    "No valid trios remain"
  )
})

test_that("LOW_MARKERS status assigned when markers_tested < min_markers", {
  # Set min_markers higher than the 20 in our test data
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE,
                           write_txt = FALSE, min_markers = 25L)

  expect_true(all(res$Status == "LOW_MARKERS"))
  expect_true(all(res$Correction_Decision == "NONE"))
})

test_that("error_threshold out of range raises an error", {
  f <- write_temp_files()
  expect_error(validate_pedigree(f$ped, f$genos, error_threshold = 150,
                                 verbose = FALSE, write_txt = FALSE))
  expect_error(validate_pedigree(f$ped, f$genos, error_threshold = -1,
                                 verbose = FALSE, write_txt = FALSE))
})

test_that("homozygous_threshold out of range raises an error", {
  f <- write_temp_files()
  expect_error(validate_pedigree(f$ped, f$genos, homozygous_threshold = 101,
                                 verbose = FALSE, write_txt = FALSE))
  expect_error(validate_pedigree(f$ped, f$genos, homozygous_threshold = -5,
                                 verbose = FALSE, write_txt = FALSE))
})

test_that("missing required pedigree column raises an error", {
  bad_ped <- data.table(Progeny = "IND_C", Parent1 = "IND_A", Dam = "IND_B")
  f       <- write_temp_files(ped = bad_ped)
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE),
    "missing required columns"
  )
})

test_that("missing ID column in genotypes raises an error", {
  bad_genos <- copy(make_genos())
  setnames(bad_genos, "ID", "SampleID")
  f <- write_temp_files(genos = bad_genos)
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE),
    "ID"
  )
})

test_that("NA markers do not cause errors and are handled gracefully", {
  genos <- make_genos()
  # Introduce NAs in a few markers for IND_C
  genos[ID == "IND_C", M1 := NA_integer_]
  genos[ID == "IND_C", M2 := NA_integer_]
  f   <- write_temp_files(genos = genos)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  expect_equal(res[Progeny == "IND_C"]$Markers_Tested, 18L)  # 2 NAs excluded
  expect_equal(res[Progeny == "IND_C"]$Status, "PASS")
})
