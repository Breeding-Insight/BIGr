# tests/testthat/test-validate_pedigree.R
# Run with: testthat::test_file("tests/testthat/test-validate_pedigree.R")
library(testthat)
library(data.table)

# ==============================================================================
# Helpers
# ==============================================================================

make_genos <- function() {
  # IND_A: all 0, IND_B: all 2, IND_C: all 1 (het), IND_D: all 0, IND_E: all 0
  n_markers    <- 20
  marker_names <- paste0("M", seq_len(n_markers))
  dt <- data.table(
    ID = c("IND_A", "IND_B", "IND_C", "IND_D", "IND_E"),
    rbind(
      rep(0L, n_markers),   # IND_A — all ref homozygous
      rep(2L, n_markers),   # IND_B — all alt homozygous
      rep(1L, n_markers),   # IND_C — all het (valid child of IND_A x IND_B)
      rep(0L, n_markers),   # IND_D — all ref (impossible child of IND_B x IND_A)
      rep(0L, n_markers)    # IND_E — all ref
    )
  )
  setnames(dt, c("ID", marker_names))
  dt
}

make_pedigree <- function() {
  # IND_C: perfect Mendelian child of IND_A x IND_B -> PASS
  # IND_D: declared parents swapped -> FAIL
  data.table(
    ID            = c("IND_C", "IND_D"),
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

# ==============================================================================
# 1. Input validation
# ==============================================================================

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
  bad_ped <- data.table(ID = "IND_C", Parent1 = "IND_A", Female_Parent = "IND_B")
  f <- write_temp_files(ped = bad_ped)
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE),
    regexp = "missing required columns"
  )
})

test_that("missing ID column in genotypes raises an error", {
  bad_genos <- copy(make_genos())
  setnames(bad_genos, "ID", "SampleID")
  f <- write_temp_files(genos = bad_genos)
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE),
    regexp = "ID"
  )
})

test_that("all trios with no genotype data stop with an error", {
  ped <- data.table(ID = "GHOST", Male_Parent = "IND_A", Female_Parent = "IND_B")
  f   <- write_temp_files(ped = ped)
  expect_error(
    validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE),
    regexp = "No valid trios remain"
  )
})

# ==============================================================================
# 2. Return structure
# ==============================================================================

test_that("returns a data.table", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_s3_class(res, "data.table")
})

test_that("result has one row per pedigree entry", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(nrow(res), 2L)
})

test_that("result has all expected columns", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expected_cols <- c(
    "ID", "Male_Parent", "Female_Parent",
    "Mendelian_Error_Pct", "Markers_Tested", "Status",
    "Correction_Decision",
    "Male_Parent_Hom_Error_Pct", "Female_Parent_Hom_Error_Pct",
    "Best_Male_Parent",   "Best_Male_Parent_Error_Pct",
    "Best_Female_Parent", "Best_Female_Parent_Error_Pct"
  )
  expect_true(all(expected_cols %in% names(res)))
})

# ==============================================================================
# 3. PASS / FAIL / LOW_MARKERS / NO_DATA statuses
# ==============================================================================

test_that("PASS trio is correctly identified", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  pass_row <- res[ID == "IND_C"]
  expect_equal(nrow(pass_row), 1L)
  expect_equal(pass_row$Status, "PASS")
  expect_equal(pass_row$Mendelian_Error_Pct, 0)
  expect_equal(pass_row$Correction_Decision, "NONE")
})

test_that("FAIL trio is correctly identified", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  fail_row <- res[ID == "IND_D"]
  expect_equal(nrow(fail_row), 1L)
  expect_equal(fail_row$Status, "FAIL")
  expect_gt(fail_row$Mendelian_Error_Pct, 5.0)
})

test_that("FAIL trio has REMOVE_MALE_PARENT decision with best match populated", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  fail_row <- res[ID == "IND_D"]
  expect_equal(fail_row$Correction_Decision, "REMOVE_MALE_PARENT")
  expect_false(is.na(fail_row$Best_Male_Parent))
  expect_true(is.na(fail_row$Best_Female_Parent))
})

test_that("Mendelian_Error_Pct is 0 for a perfect trio", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(res[ID == "IND_C"]$Mendelian_Error_Pct, 0)
})

test_that("Markers_Tested equals number of markers for complete data", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(res[ID == "IND_C"]$Markers_Tested, 20L)
})

test_that("LOW_MARKERS status assigned when markers_tested < min_markers", {
  f   <- write_temp_files()
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE,
                           write_txt = FALSE, min_markers = 25L)
  expect_true(all(res$Status == "LOW_MARKERS"))
  expect_true(all(grepl("^LOW_MARKERS_", res$Correction_Decision)))
})

test_that("NA markers reduce Markers_Tested and do not cause errors", {
  genos <- make_genos()
  genos[ID == "IND_C", M1 := NA_integer_]
  genos[ID == "IND_C", M2 := NA_integer_]
  f   <- write_temp_files(genos = genos)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(res[ID == "IND_C"]$Markers_Tested, 18L)
  expect_equal(res[ID == "IND_C"]$Status, "PASS")
})

# ==============================================================================
# 4. Missing parent statuses (MISSING_MALE_PARENT / MISSING_FEMALE_PARENT /
#    MISSING_BOTH_PARENTS)
# Note: each test includes make_pedigree() rows so has_geno is never empty,
#       and filters res by the specific ID under test [3][4]
# ==============================================================================

test_that("MISSING_MALE_PARENT status and recommendation are correct", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_E", Male_Parent = "0", Female_Parent = "IND_B")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  r   <- res[ID == "IND_E"]
  expect_equal(r$Status,              "MISSING_MALE_PARENT")
  expect_equal(r$Correction_Decision, "NONE")
  expect_false(is.na(r$Best_Male_Parent))
  expect_true(is.na(r$Best_Female_Parent))
})

test_that("MISSING_FEMALE_PARENT status and recommendation are correct", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_E", Male_Parent = "IND_A", Female_Parent = "0")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  r   <- res[ID == "IND_E"]
  expect_equal(r$Status,              "MISSING_FEMALE_PARENT")
  expect_equal(r$Correction_Decision, "NONE")
  expect_true(is.na(r$Best_Male_Parent))
  expect_false(is.na(r$Best_Female_Parent))
})

test_that("MISSING_BOTH_PARENTS status and recommendations are correct", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_E", Male_Parent = "0", Female_Parent = "0")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  r   <- res[ID == "IND_E"]
  expect_equal(r$Status,              "MISSING_BOTH_PARENTS")
  expect_equal(r$Correction_Decision, "NONE")
  expect_false(is.na(r$Best_Male_Parent))
  expect_false(is.na(r$Best_Female_Parent))
})

test_that("MISSING_* rows preserve 0s in corrected pedigree", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_E", Male_Parent = "0", Female_Parent = "IND_B")
  )
  f       <- write_temp_files(ped = ped)
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)
  validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  corr <- fread(file.path(tmp_dir, "corrected_pedigree.txt"),
                colClasses = "character")
  expect_equal(corr[ID == "IND_E"]$Male_Parent, "0")
})

test_that("Best_Male_Parent for MISSING_MALE_PARENT is excluded from being the known female", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_E", Male_Parent = "0", Female_Parent = "IND_B")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  r   <- res[ID == "IND_E"]
  # The known female parent should not be suggested as the best male parent
  expect_false(r$Best_Male_Parent == "IND_B")
})

# ==============================================================================
# 5. FOUNDERS status
# ==============================================================================

test_that("FOUNDERS status is assigned when ID in founders list with 0 0 parents", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_A", Male_Parent = "0", Female_Parent = "0")
  )
  founders_file <- tempfile(fileext = ".txt")
  fwrite(data.table(ID = "IND_A"), founders_file,
         sep = "\t", quote = FALSE, col.names = FALSE)
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos,
                           founders_file = founders_file,
                           verbose       = FALSE,
                           write_txt     = FALSE)
  founder_row <- res[ID == "IND_A"]
  expect_equal(founder_row$Status,              "FOUNDERS")
  expect_equal(founder_row$Correction_Decision, "NONE")
  expect_true(is.na(founder_row$Best_Male_Parent))
  expect_true(is.na(founder_row$Best_Female_Parent))
})

test_that("non-founder rows are still evaluated normally when founders file is supplied", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_A", Male_Parent = "0", Female_Parent = "0")
  )
  founders_file <- tempfile(fileext = ".txt")
  fwrite(data.table(ID = "IND_A"), founders_file,
         sep = "\t", quote = FALSE, col.names = FALSE)
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos,
                           founders_file = founders_file,
                           verbose       = FALSE,
                           write_txt     = FALSE)
  # IND_C has real parents — should still get PASS
  expect_equal(res[ID == "IND_C"]$Status, "PASS")
})

test_that("0 0 parents NOT in founders list get MISSING_BOTH_PARENTS", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_E", Male_Parent = "0", Female_Parent = "0")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(res[ID == "IND_E"]$Status, "MISSING_BOTH_PARENTS")
})

test_that("0 0 parents with no founders file gets MISSING_BOTH_PARENTS not FOUNDERS", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_E", Male_Parent = "0", Female_Parent = "0")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos,
                           founders_file = NULL,
                           verbose       = FALSE,
                           write_txt     = FALSE)
  expect_equal(res[ID == "IND_E"]$Status, "MISSING_BOTH_PARENTS")
})

# ==============================================================================
# 6. NO_GENOTYPE_DATA status
# ==============================================================================

test_that("NO_GENOTYPE_DATA is flagged for progeny absent from genotype file", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "GHOST", Male_Parent = "IND_A", Female_Parent = "IND_B")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  ghost_row <- res[ID == "GHOST"]
  expect_equal(nrow(ghost_row),               1L)
  expect_equal(ghost_row$Status,              "NO_GENOTYPE_DATA")
  expect_equal(ghost_row$Correction_Decision, "NONE")
})

test_that("NO_GENOTYPE_DATA rows have NA for all analysis columns", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "GHOST", Male_Parent = "IND_A", Female_Parent = "IND_B")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  ghost_row <- res[ID == "GHOST"]
  expect_true(is.na(ghost_row$Mendelian_Error_Pct))
  expect_equal(ghost_row$Markers_Tested, 0L)
  expect_true(is.na(ghost_row$Best_Male_Parent))
  expect_true(is.na(ghost_row$Best_Female_Parent))
})

test_that("NO_GENOTYPE_DATA flagged when declared parent is absent from genotype file", {
  ped <- rbind(
    make_pedigree(),   # ensures has_geno is not empty
    data.table(ID = "IND_C_GHOST", Male_Parent = "GHOST_DAD", Female_Parent = "IND_B")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(res[ID == "IND_C_GHOST"]$Status, "NO_GENOTYPE_DATA")
})

test_that("valid trios still evaluated correctly when ghost rows are present", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "GHOST", Male_Parent = "IND_A", Female_Parent = "IND_B")
  )
  f   <- write_temp_files(ped = ped)
  res <- validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_equal(res[ID == "IND_C"]$Status, "PASS")
  expect_equal(res[ID == "IND_D"]$Status, "FAIL")
})

# ==============================================================================
# 7. Corrected pedigree output
# ==============================================================================

test_that("corrected_pedigree.txt is written and PASS parents are unchanged", {
  f       <- write_temp_files()
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)
  validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  expect_true(file.exists(file.path(tmp_dir, "corrected_pedigree.txt")))
  corr <- fread(file.path(tmp_dir, "corrected_pedigree.txt"))
  # IND_C passes — parents must be unchanged
  expect_equal(as.character(corr[ID == "IND_C"]$Male_Parent),   "IND_A")
  expect_equal(as.character(corr[ID == "IND_C"]$Female_Parent), "IND_B")
})

test_that("corrected_pedigree.txt sets bad parent to 0 for FAIL trio", {
  f       <- write_temp_files()
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)
  validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  corr <- fread(file.path(tmp_dir, "corrected_pedigree.txt"),
                colClasses = "character")
  # IND_D fails with REMOVE_MALE_PARENT — male should become "0"
  expect_equal(corr[ID == "IND_D"]$Male_Parent,   "0")
  expect_equal(corr[ID == "IND_D"]$Female_Parent, "IND_A")
})

test_that("corrected_pedigree.txt preserves 0s for MISSING_* rows", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "IND_E", Male_Parent = "0", Female_Parent = "IND_B")
  )
  f       <- write_temp_files(ped = ped)
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd  <- getwd()
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)
  validate_pedigree(f$ped, f$genos, verbose = FALSE, write_txt = FALSE)
  corr <- fread(file.path(tmp_dir, "corrected_pedigree.txt"),
                colClasses = "character")
  # MISSING_MALE_PARENT — male stays "0", female unchanged
  expect_equal(corr[ID == "IND_E"]$Male_Parent,   "0")
  expect_equal(corr[ID == "IND_E"]$Female_Parent, "IND_B")
})

# ==============================================================================
# 8. write_txt / output file
# ==============================================================================

test_that("write_txt = TRUE writes output file with correct row count", {
  f        <- write_temp_files()
  out_file <- tempfile(fileext = ".txt")
  validate_pedigree(f$ped, f$genos, verbose = FALSE,
                    write_txt = TRUE, output_filename = out_file)
  expect_true(file.exists(out_file))
  written <- fread(out_file)
  expect_equal(nrow(written), 2L)
})

test_that("write_txt = FALSE does not create output file", {
  f        <- write_temp_files()
  out_file <- tempfile(fileext = ".txt")
  validate_pedigree(f$ped, f$genos, verbose = FALSE,
                    write_txt = FALSE, output_filename = out_file)
  expect_false(file.exists(out_file))
})

test_that("output file contains correct number of rows when ghost rows present", {
  ped <- rbind(
    make_pedigree(),
    data.table(ID = "GHOST", Male_Parent = "IND_A", Female_Parent = "IND_B")
  )
  f        <- write_temp_files(ped = ped)
  out_file <- tempfile(fileext = ".txt")
  validate_pedigree(f$ped, f$genos, verbose = FALSE,
                    write_txt = TRUE, output_filename = out_file)
  written <- fread(out_file)
  # 2 valid trios + 1 ghost = 3 rows total
  expect_equal(nrow(written), 3L)
})
