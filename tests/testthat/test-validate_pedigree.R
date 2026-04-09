#### Tests for validate_pedigree() ####
# Run with: testthat::test_file("test-validate_pedigree.R")
library(testthat)
library(data.table)

#### Helpers ####
make_genos <- function() {
  n_markers    <- 20
  marker_names <- paste0("M", seq_len(n_markers))
  pa <- rep(0L, n_markers)
  pb <- rep(2L, n_markers)
  pc <- rep(1L, n_markers)
  pd <- rep(0L, n_markers)
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
    ID            = c("IND_C", "IND_D"),   # <-- changed from Progeny
    Male_Parent   = c("IND_A", "IND_B"),
    Female_Parent = c("IND_B", "IND_A")
  )
}

write_temp_files <- function(genos = make_genos(), ped = make_pedigree()) {
  ped_file   <- tempfile(fileext = ".txt")
  genos_file <- tempfile(fileext = ".txt")
  fwrite(ped,   ped_file,   sep = "\t")
  fwrite(genos, genos_file, sep = "\t")
  list(ped = ped_file, genos = genos_file)
}

#### Tests ####

test_that("PASS trio is correctly identified", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  pass_row <- res[ID == "IND_C"]                                  # <-- changed
  expect_equal(nrow(pass_row), 1L)
  expect_equal(pass_row$Status, "PASS")
  expect_equal(pass_row$Mendelian_Error_Pct, 0)
  expect_equal(pass_row$Correction_Decision, "NONE")
})

test_that("FAIL trio is correctly identified", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  fail_row <- res[ID == "IND_D"]                                  # <-- changed
  expect_equal(nrow(fail_row), 1L)
  expect_equal(fail_row$Status, "FAIL")
  expect_gt(fail_row$Mendelian_Error_Pct, 5.0)
})

test_that("FAIL trio has correct correction decision (REMOVE_MALE_PARENT)", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  fail_row <- res[ID == "IND_D"]                                  # <-- changed
  expect_equal(fail_row$Correction_Decision, "REMOVE_MALE_PARENT")
  expect_false(is.na(fail_row$Best_Male_Parent))
  expect_true(is.na(fail_row$Best_Female_Parent))
})

test_that("Mendelian_Error_Pct is 0 for perfect trio", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(res[ID == "IND_C"]$Mendelian_Error_Pct, 0)        # <-- changed
})

test_that("Markers_Tested equals number of markers for complete data", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(res[ID == "IND_C"]$Markers_Tested, 20L)           # <-- changed
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
  expected_cols <- c(
    "ID", "Male_Parent", "Female_Parent",                        # <-- changed
    "Mendelian_Error_Pct", "Markers_Tested", "Status",
    "Correction_Decision",
    "Male_Parent_Hom_Error_Pct", "Female_Parent_Hom_Error_Pct",
    "Best_Male_Parent",   "Best_Male_Parent_Error_Pct",
    "Best_Female_Parent", "Best_Female_Parent_Error_Pct"
  )
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
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

  validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)

  corr <- fread(file.path(tmp_dir, "corrected_pedigree.txt"))
  expect_equal(corr[ID == "IND_D"]$Male_Parent,   "0")           # <-- changed
  expect_equal(corr[ID == "IND_D"]$Female_Parent, "IND_A")       # <-- changed
  expect_equal(corr[ID == "IND_C"]$Male_Parent,   "IND_A")       # <-- changed
  expect_equal(corr[ID == "IND_C"]$Female_Parent, "IND_B")       # <-- changed
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

  filled       <- fread(file.path(tmp_dir, "filled_pedigree.txt"))
  new_male_par <- filled[ID == "IND_D"]$Male_Parent               # <-- changed
  expect_false(new_male_par == "IND_B")
  expect_false(new_male_par == "0")
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

test_that("Trios with missing genotype data are removed and error is thrown", {
  ped <- data.table(
    ID            = "GHOST",                                       # <-- changed
    Male_Parent   = "IND_A",
    Female_Parent = "IND_B"
  )
  f <- write_temp_files(ped = ped)
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE),
    "No valid trios remain"
  )
})

test_that("LOW_MARKERS status assigned when markers_tested < min_markers", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE,
                           write_txt = FALSE, min_markers = 25L)
  expect_true(all(res$Status == "LOW_MARKERS"))
  expect_true(all(res$Correction_Decision == "NONE"))
})

test_that("trio_error_threshold out of range raises an error", {
  f <- write_temp_files()
  expect_error(validate_pedigree(f$ped, f$genos,
                                 trio_error_threshold = 150,
                                 verbose = FALSE, write_txt = FALSE))
  expect_error(validate_pedigree(f$ped, f$genos,
                                 trio_error_threshold = -1,
                                 verbose = FALSE, write_txt = FALSE))
})

test_that("single_parent_error_threshold out of range raises an error", {
  f <- write_temp_files()
  expect_error(validate_pedigree(f$ped, f$genos,
                                 single_parent_error_threshold = 101,
                                 verbose = FALSE, write_txt = FALSE))
  expect_error(validate_pedigree(f$ped, f$genos,
                                 single_parent_error_threshold = -5,
                                 verbose = FALSE, write_txt = FALSE))
})

test_that("missing required pedigree column raises an error", {
  bad_ped <- data.table(
    ID            = "IND_C",                                       # <-- changed
    Parent1       = "IND_A",
    Female_Parent = "IND_B"
  )
  f <- write_temp_files(ped = bad_ped)
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
  genos[ID == "IND_C", M1 := NA_integer_]
  genos[ID == "IND_C", M2 := NA_integer_]
  f   <- write_temp_files(genos = genos)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(res[ID == "IND_C"]$Markers_Tested, 18L)            # <-- changed
  expect_equal(res[ID == "IND_C"]$Status, "PASS")                 # <-- changed
})
