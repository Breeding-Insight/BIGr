context("check_ped – Pedigree Quality Checks")

# ---------------------------------------------------------------------------
# Helper: write a data.frame to a temp tab-separated file and return the path
# ---------------------------------------------------------------------------
write_ped <- function(df) {
  f <- tempfile(fileext = ".txt")
  utils::write.table(df, f, sep = "\t", row.names = FALSE, quote = FALSE)
  f
}

# ===========================================================================
# 1. Return structure
# Tests that check_ped returns a named list of exactly 5 data frame
# components covering all five issue categories.
# ===========================================================================
test_that("check_ped returns a named list of length 5", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("0", "A", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), seed = 1, verbose = FALSE, correct = FALSE)

  expect_type(out, "list")
  expect_length(out, 5)
  expect_named(out, c(
    "exact_duplicates",
    "repeated_ids_diff",
    "inconsistent_sex_roles",
    "missing_parents",
    "dependencies"
  ))
})

# ===========================================================================
# 2. Component types
# Tests that every element in the returned list is a data.frame, ensuring
# downstream code can safely call nrow(), colnames(), etc. on any component.
# ===========================================================================
test_that("check_ped report components are data.frames", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("0", "A", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE, correct = FALSE)

  expect_true(is.data.frame(out$exact_duplicates))
  expect_true(is.data.frame(out$repeated_ids_diff))
  expect_true(is.data.frame(out$inconsistent_sex_roles))
  expect_true(is.data.frame(out$missing_parents))
  expect_true(is.data.frame(out$dependencies))
})

# ===========================================================================
# 3. Clean pedigree → empty report
# Tests that a structurally valid pedigree with no issues produces zero rows
# in every report component (no false positives).
# ===========================================================================
test_that("clean pedigree produces no issues", {
  ped <- data.frame(
    id            = c("G1", "G2", "P1"),
    male_parent   = c("0",  "0",  "G1"),
    female_parent = c("0",  "0",  "G2")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE, correct = FALSE)

  expect_equal(nrow(out$exact_duplicates),       0)
  expect_equal(nrow(out$repeated_ids_diff),       0)
  expect_equal(nrow(out$inconsistent_sex_roles),  0)
  expect_equal(nrow(out$missing_parents),         0)
  expect_equal(nrow(out$dependencies),            0)
})

# ===========================================================================
# 4. Exact duplicate detection and correction
# Tests that rows sharing the same (id, male_parent, female_parent) triplet
# are flagged, and that correct = TRUE retains only one copy of each.
# ===========================================================================
test_that("check_ped detects exact duplicates", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("0", "0", "A"),
    female_parent = c("0", "0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE, correct = FALSE)

  expect_equal(nrow(out$exact_duplicates), 2)
  expect_true(all(out$exact_duplicates$id == "A"))
})

test_that("correct = TRUE removes exact duplicates", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("0", "0", "A"),
    female_parent = c("0", "0", "0")
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  corr_name  <- paste0(file_base, "_corrected")
  rep_name   <- paste0(file_base, "_report")
  on.exit({
    if (exists(corr_name, envir = .GlobalEnv)) rm(list = corr_name, envir = .GlobalEnv)
    if (exists(rep_name,  envir = .GlobalEnv)) rm(list = rep_name,  envir = .GlobalEnv)
  })

  check_ped(f, verbose = FALSE, correct = TRUE)
  corrected <- get(corr_name, envir = .GlobalEnv)

  # Only one "A" row should remain
  expect_equal(sum(corrected$id == "A"), 1)
})

# ===========================================================================
# 5. Repeated IDs with conflicting parent assignments
# Tests that the same ID appearing on multiple rows with different parent
# values is flagged, and that correct = TRUE blanks only the conflicting
# parent field while preserving consistent values.
# ===========================================================================
test_that("check_ped detects repeated IDs with conflicting parents", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("X", "Y", "A"),
    female_parent = c("M", "M", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE, correct = FALSE)

  expect_equal(nrow(out$repeated_ids_diff), 2)
  expect_true(all(out$repeated_ids_diff$id == "A"))
})

test_that("correct = TRUE resolves conflicting IDs: conflicting field -> '0', consistent field kept", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("X", "Y", "A"),   # conflicting -> should become "0"
    female_parent = c("M", "M", "0")    # consistent  -> should stay "M"
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  corr_name  <- paste0(file_base, "_corrected")
  rep_name   <- paste0(file_base, "_report")
  on.exit({
    if (exists(corr_name, envir = .GlobalEnv)) rm(list = corr_name, envir = .GlobalEnv)
    if (exists(rep_name,  envir = .GlobalEnv)) rm(list = rep_name,  envir = .GlobalEnv)
  })

  check_ped(f, verbose = FALSE, correct = TRUE)
  corrected <- get(corr_name, envir = .GlobalEnv)
  a_row     <- corrected[corrected$id == "A", ]

  expect_equal(nrow(a_row), 1)
  expect_equal(a_row$male_parent,   "0")
  expect_equal(a_row$female_parent, "M")
})

# ===========================================================================
# 6. Inconsistent parent sex roles
# Tests that an individual ID appearing in both male_parent and female_parent
# columns is flagged. Also confirms that correct = TRUE does NOT auto-resolve
# this issue since the correct fix requires manual review.
# ===========================================================================
test_that("check_ped detects inconsistent parent sex roles", {
  ped <- data.frame(
    id            = c("child1", "child2", "P"),
    male_parent   = c("P",      "0",      "0"),
    female_parent = c("0",      "P",      "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE, correct = FALSE)

  expect_true("inconsistent_sex_roles" %in% names(out))
  expect_gt(nrow(out$inconsistent_sex_roles), 0)

  # Both rows that reference P should be flagged
  expect_true(any(out$inconsistent_sex_roles$male_parent   == "P" |
                    out$inconsistent_sex_roles$female_parent == "P"))
})

test_that("correct = TRUE does NOT auto-correct inconsistent sex roles", {
  ped <- data.frame(
    id            = c("child1", "child2", "P"),
    male_parent   = c("P",      "0",      "0"),
    female_parent = c("0",      "P",      "0")
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  corr_name  <- paste0(file_base, "_corrected")
  rep_name   <- paste0(file_base, "_report")
  on.exit({
    if (exists(corr_name, envir = .GlobalEnv)) rm(list = corr_name, envir = .GlobalEnv)
    if (exists(rep_name,  envir = .GlobalEnv)) rm(list = rep_name,  envir = .GlobalEnv)
  })

  check_ped(f, verbose = FALSE, correct = TRUE)
  corrected <- get(corr_name, envir = .GlobalEnv)

  # P should still appear with its original (inconsistent) parent assignments
  expect_true(any(corrected$male_parent   == "P" |
                    corrected$female_parent == "P"))
})

# ===========================================================================
# 7. Missing parent detection and founder row injection
# Tests that parent IDs referenced in male_parent/female_parent but absent
# from the id column are flagged. Confirms that correct = TRUE appends a
# founder row (both parents = "0") for each missing parent.
# ===========================================================================
test_that("check_ped detects missing parents", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "X"),
    female_parent = c("0", "Y")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE, correct = FALSE)

  expect_equal(nrow(out$missing_parents), 1)
  expect_true(all(out$missing_parents$id == "B"))
})

test_that("correct = TRUE adds missing parents as founder rows", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "X"),
    female_parent = c("0", "Y")
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  corr_name  <- paste0(file_base, "_corrected")
  rep_name   <- paste0(file_base, "_report")
  on.exit({
    if (exists(corr_name, envir = .GlobalEnv)) rm(list = corr_name, envir = .GlobalEnv)
    if (exists(rep_name,  envir = .GlobalEnv)) rm(list = rep_name,  envir = .GlobalEnv)
  })

  check_ped(f, verbose = FALSE, correct = TRUE)
  corrected <- get(corr_name, envir = .GlobalEnv)

  # X and Y should now appear as founder rows
  expect_true("X" %in% corrected$id)
  expect_true("Y" %in% corrected$id)
  x_row <- corrected[corrected$id == "X", ]
  y_row <- corrected[corrected$id == "Y", ]
  expect_equal(x_row$male_parent,   "0")
  expect_equal(x_row$female_parent, "0")
  expect_equal(y_row$male_parent,   "0")
  expect_equal(y_row$female_parent, "0")
})

# ===========================================================================
# 8. Pedigree cycle / dependency detection
# Tests that circular ancestry chains (e.g. A is parent of B and B is parent
# of A) are detected and reported. Also confirms that correct = TRUE does NOT
# silently remove cycle-involved individuals, as cycles require manual review.
# ===========================================================================
test_that("check_ped detects a direct two-node cycle", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("B", "A"),
    female_parent = c("0", "0")
  )
  out <- check_ped(write_ped(ped), verbose = FALSE, correct = FALSE)

  expect_gt(nrow(out$dependencies), 0)
  expect_true(all(c("A", "B") %in% out$dependencies$id))
})

test_that("check_ped does NOT auto-correct cycles", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("B", "A"),
    female_parent = c("0", "0")
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  corr_name  <- paste0(file_base, "_corrected")
  rep_name   <- paste0(file_base, "_report")
  on.exit({
    if (exists(corr_name, envir = .GlobalEnv)) rm(list = corr_name, envir = .GlobalEnv)
    if (exists(rep_name,  envir = .GlobalEnv)) rm(list = rep_name,  envir = .GlobalEnv)
  })

  check_ped(f, verbose = FALSE, correct = TRUE)
  corrected <- get(corr_name, envir = .GlobalEnv)

  # Cycle-involved IDs should still be present (no silent removal)
  expect_true("A" %in% corrected$id)
  expect_true("B" %in% corrected$id)
})

# ===========================================================================
# 9. Global environment assignments
# Tests that the report list is always written to the global environment,
# that the corrected pedigree is written only when correct = TRUE, and that
# the internal row_number column does not leak into the corrected output.
# ===========================================================================
test_that("report is always assigned to global environment", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "A"),
    female_parent = c("0", "0")
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  rep_name   <- paste0(file_base, "_report")
  on.exit(if (exists(rep_name, envir = .GlobalEnv))
    rm(list = rep_name, envir = .GlobalEnv))

  check_ped(f, verbose = FALSE, correct = FALSE)

  expect_true(exists(rep_name, envir = .GlobalEnv))
  expect_type(get(rep_name, envir = .GlobalEnv), "list")
})

test_that("corrected pedigree is assigned to global env when correct = TRUE", {
  ped <- data.frame(
    id            = c("A", "A", "B"),
    male_parent   = c("X", "Y", "A"),
    female_parent = c("M", "M", "0")
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  corr_name  <- paste0(file_base, "_corrected")
  rep_name   <- paste0(file_base, "_report")
  on.exit({
    if (exists(corr_name, envir = .GlobalEnv)) rm(list = corr_name, envir = .GlobalEnv)
    if (exists(rep_name,  envir = .GlobalEnv)) rm(list = rep_name,  envir = .GlobalEnv)
  })

  check_ped(f, verbose = FALSE, correct = TRUE)

  expect_true(exists(corr_name, envir = .GlobalEnv))
  corrected <- get(corr_name, envir = .GlobalEnv)
  expect_true(is.data.frame(corrected))

  # row_number must not leak into corrected output
  expect_false("row_number" %in% names(corrected))
  expect_true(all(c("id", "male_parent", "female_parent") %in% names(corrected)))
})

test_that("corrected pedigree is NOT assigned to global env when correct = FALSE", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "A"),
    female_parent = c("0", "0")
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  corr_name  <- paste0(file_base, "_corrected")
  rep_name   <- paste0(file_base, "_report")
  on.exit({
    if (exists(corr_name, envir = .GlobalEnv)) rm(list = corr_name, envir = .GlobalEnv)
    if (exists(rep_name,  envir = .GlobalEnv)) rm(list = rep_name,  envir = .GlobalEnv)
  })

  check_ped(f, verbose = FALSE, correct = FALSE)

  expect_false(exists(corr_name, envir = .GlobalEnv))
  expect_true(exists(rep_name,   envir = .GlobalEnv))
})

# ===========================================================================
# 10. Missing required columns → hard stop
# Tests that a file lacking the canonical column names (after alias remapping)
# causes an immediate error rather than silently producing empty results.
# Column names like animal_id/parent1/parent2 are not in the alias list and
# therefore cannot be remapped, triggering the validation stop.
# ===========================================================================
test_that("check_ped errors when required columns are missing", {
  bad_df <- data.frame(
    animal_id  = c("a", "b"),
    parent1    = c("0", "a"),
    parent2    = c("0", "0")
  )
  expect_error(
    check_ped(write_ped(bad_df), verbose = FALSE),
    "Input file is missing required column"
  )
})

# ===========================================================================
# 11. Integration test using the bundled check_ped_test.txt fixture
# End-to-end test against a real fixture file shipped with the package.
# Validates the overall row counts and that known problematic individuals
# (grandfather2, grandfather3) are flagged for inconsistent sex roles.
# ===========================================================================
test_that("integration test with bundled fixture file", {
  ped_file <- system.file("check_ped_test.txt", package = "BIGr")
  skip_if(ped_file == "", "Bundled fixture file not found; skipping integration test.")

  out <- check_ped(ped_file, seed = 101919, verbose = FALSE)

  expect_length(out, 5)

  # Inconsistent sex roles should flag grandfather2 & grandfather3
  expect_true(all(
    c("grandfather2", "grandfather3") %in% out$inconsistent_sex_roles$male_parent |
      c("grandfather2", "grandfather3") %in% out$inconsistent_sex_roles$female_parent
  ))
  expect_equal(nrow(out$missing_parents), 8)
})

# ===========================================================================
# 12. Zip export (save_zip = TRUE / save_corrected_zip = TRUE)
# Tests that save_zip = TRUE produces a .zip archive in the working directory,
# and that save_corrected_zip = TRUE includes the corrected pedigree file
# inside that archive when correct = TRUE.
# ===========================================================================
test_that("save_zip = TRUE creates a zip archive in the working directory", {
  ped <- data.frame(
    id            = c("A", "B", "C"),
    male_parent   = c("0", "A", "A"),
    female_parent = c("0", "0", "0")
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  zip_name   <- paste0(file_base, "_report.zip")
  rep_name   <- paste0(file_base, "_report")
  zip_path   <- file.path(getwd(), zip_name)

  on.exit({
    if (file.exists(zip_path))                file.remove(zip_path)
    if (exists(rep_name, envir = .GlobalEnv)) rm(list = rep_name, envir = .GlobalEnv)
  })

  check_ped(f, verbose = FALSE, correct = FALSE, save_zip = TRUE)

  expect_true(file.exists(zip_path))
})

test_that("save_zip + save_corrected_zip + correct = TRUE includes corrected file in zip", {
  ped <- data.frame(
    id            = c("A", "B"),
    male_parent   = c("0", "X"),
    female_parent = c("0", "Y")
  )
  f          <- write_ped(ped)
  file_base  <- tools::file_path_sans_ext(basename(f))
  zip_name   <- paste0(file_base, "_report.zip")
  corr_name  <- paste0(file_base, "_corrected")
  rep_name   <- paste0(file_base, "_report")
  zip_path   <- file.path(getwd(), zip_name)

  on.exit({
    if (file.exists(zip_path))                 file.remove(zip_path)
    if (exists(corr_name, envir = .GlobalEnv)) rm(list = corr_name, envir = .GlobalEnv)
    if (exists(rep_name,  envir = .GlobalEnv)) rm(list = rep_name,  envir = .GlobalEnv)
  })

  check_ped(f, verbose = FALSE, correct = TRUE,
            save_zip = TRUE, save_corrected_zip = TRUE)

  expect_true(file.exists(zip_path))
  zip_contents <- utils::unzip(zip_path, list = TRUE)$Name
  expect_true(any(grepl("_corrected", zip_contents)))
})
