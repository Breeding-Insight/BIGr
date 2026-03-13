library(testthat)
library(dplyr)
library(stringr)

# ── Helpers: create temporary test fixtures ──────────────────────────────────

#' Build a raw (unprocessed) MADC CSV with 7 filler rows of "*"
create_raw_madc <- function(path, clone_ids, allele_ids, sample_names, counts_matrix,
                            include_summary_cols = TRUE) {
  n_rows <- length(clone_ids)
  n_samples <- length(sample_names)

  # Summary columns matching the full set the function auto-removes
  summary_col_names <- c("ClusterConsensusSequence", "CallRate", "OneRatioRef",
                         "OneRatioSnp", "FreqHomRef", "FreqHomSnp", "FreqHets",
                         "PICRef", "PICSnp", "AvgPIC", "AvgCountRef",
                         "AvgCountSnp", "RatioAvgCountRefAvgCountSnp")
  extra_cols <- if (include_summary_cols) length(summary_col_names) else 0
  total_cols <- 3 + extra_cols + n_samples  # AlleleID, CloneID, AlleleSequence + summary + samples

  # 7 filler rows: col 1 = "*", other cols have placeholder metadata (like real DArT files)
  filler <- matrix("*", nrow = 7, ncol = total_cols)

  # Header row — real MADC: AlleleID (col1), CloneID (col2), AlleleSequence (col3)
  summary_names <- if (include_summary_cols) summary_col_names else character(0)
  header <- c("AlleleID", "CloneID", "AlleleSequence", summary_names, sample_names)

  # Data rows
  allele_seqs <- rep("ACTCTCAGGTGGAT", n_rows)
  data_rows <- cbind(
    AlleleID = allele_ids,
    CloneID = clone_ids,
    AlleleSequence = allele_seqs,
    if (include_summary_cols) {
      data.frame(
        ClusterConsensusSequence = rep("ACTCTCAGGTGGAT", n_rows),
        CallRate = 1.0, OneRatioRef = 1.0, OneRatioSnp = 1.0,
        FreqHomRef = 0.5, FreqHomSnp = 0.3, FreqHets = 0.2,
        PICRef = 0.4, PICSnp = 0.3, AvgPIC = 0.35,
        AvgCountRef = 100, AvgCountSnp = 80,
        RatioAvgCountRefAvgCountSnp = 1.25,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    },
    counts_matrix
  )

  # Write filler rows first (no header), then header + data
  con <- file(path, open = "wt")
  for (i in seq_len(nrow(filler))) {
    writeLines(paste(filler[i, ], collapse = ","), con)
  }
  writeLines(paste(header, collapse = ","), con)
  for (i in seq_len(nrow(data_rows))) {
    writeLines(paste(data_rows[i, ], collapse = ","), con)
  }
  close(con)
  invisible(path)
}

#' Build a preprocessed (no filler rows) MADC CSV
create_preprocessed_madc <- function(path, clone_ids, allele_ids, sample_names, counts_matrix,
                                     include_summary_cols = TRUE) {
  n_rows <- length(clone_ids)

  df <- data.frame(
    AlleleID = allele_ids,
    CloneID = clone_ids,
    AlleleSequence = rep("ACTCTCAGGTGGAT", n_rows),
    stringsAsFactors = FALSE
  )
  if (include_summary_cols) {
    df$ClusterConsensusSequence <- rep("ACTCTCAGGTGGAT", n_rows)
    df$CallRate <- 1.0
    df$OneRatioRef <- 1.0
    df$OneRatioSnp <- 1.0
    df$FreqHomRef <- 0.5
    df$FreqHomSnp <- 0.3
    df$FreqHets <- 0.2
    df$PICRef <- 0.4
    df$PICSnp <- 0.3
    df$AvgPIC <- 0.35
    df$AvgCountRef <- 100
    df$AvgCountSnp <- 80
    df$RatioAvgCountRefAvgCountSnp <- 1.25
  }
  df <- cbind(df, counts_matrix)
  colnames(df)[(ncol(df) - ncol(counts_matrix) + 1):ncol(df)] <- sample_names

  write.csv(df, path, row.names = FALSE)
  invisible(path)
}

#' Build a 3-column marker CSV
create_marker_file <- function(path, clone_ids, chromosomes, positions) {
  df <- data.frame(CloneID = clone_ids, Chr = chromosomes, Pos = positions,
                   stringsAsFactors = FALSE)
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}


# ── Shared fixtures ──────────────────────────────────────────────────────────

sample_names  <- c("SampleA", "SampleB", "SampleC")
clone_ids     <- c("marker1", "marker2", "marker3", "marker4")
allele_ids    <- paste0(clone_ids, "|Ref")
counts_matrix <- data.frame(
  SampleA = c(10, 20, 30, 40),
  SampleB = c(15, 25, 35, 45),
  SampleC = c(12, 22, 32, 42)
)























# ── 1. Basic happy-path with raw MADC (7 filler rows) ───────────────────────

test_that("fixMADC correctly processes a raw MADC file and replaces IDs", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix,
                  include_summary_cols = TRUE)
  create_marker_file(marker_path, clone_ids,
                     c("Chr01", "Chr01", "Chr02", "Chr02"),
                     c(100, 200, 300, 400))

  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path)

  # Should return a data.frame

  expect_s3_class(result, "data.frame")


  # CloneIDs should now be Chr_Pos format

  expect_equal(unname(result$CloneID), c("Chr01_000000100", "Chr01_000000200", "Chr02_000000300", "Chr02_000000400"))

  # AlleleID should have the new prefix with |Ref_0001 suffix
  expect_true(all(grepl("\\|Ref_0001$", result$AlleleID)))
  expect_equal(unname(result$AlleleID), c("Chr01_000000100|Ref_0001", "Chr01_000000200|Ref_0001", "Chr02_000000300|Ref_0001", "Chr02_000000400|Ref_0001"))

  # Summary columns should be removed
  summary_cols <- c("ClusterConsensusSequence", "CallRate", "OneRatioRef",
                    "OneRatioSnp", "FreqHomRef", "FreqHomSnp", "FreqHets",
                    "PICRef", "PICSnp", "AvgPIC", "AvgCountRef",
                    "AvgCountSnp", "RatioAvgCountRefAvgCountSnp")
  expect_false(any(summary_cols %in% colnames(result)))

  # Sample columns should still be present
  expect_true(all(sample_names %in% colnames(result)))
})


# ── 1b. Marker file in different order than MADC ─────────────────────────────

test_that("fixMADC correctly maps IDs when marker file order differs from MADC", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix,
                  include_summary_cols = TRUE)

  # Marker file in reversed order relative to MADC
  create_marker_file(marker_path,
                     rev(clone_ids),
                     c("Chr02", "Chr02", "Chr01", "Chr01"),
                     c(400, 300, 200, 100))

  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path)

  # marker1 -> Chr01_100, marker2 -> Chr01_200, marker3 -> Chr02_300, marker4 -> Chr02_400
  # (same mapping regardless of marker file row order)
  expect_equal(unname(result$CloneID), c("Chr01_000000100", "Chr01_000000200", "Chr02_000000300", "Chr02_000000400"))
  expect_equal(unname(result$AlleleID), c("Chr01_000000100|Ref_0001", "Chr01_000000200|Ref_0001", "Chr02_000000300|Ref_0001", "Chr02_000000400|Ref_0001"))
})


# ── 1c. Duplicate CloneIDs with different AlleleID suffixes ──────────────────

test_that("fixMADC correctly handles multiple rows per CloneID (Ref/Alt pairs)", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  # Two markers, each with a Ref and Alt row = 4 rows total
  dup_clone_ids  <- c("marker1", "marker1", "marker2", "marker2")
  dup_allele_ids <- c("marker1|Ref", "marker1|Alt", "marker2|Ref", "marker2|Alt")
  dup_counts     <- data.frame(
    SampleA = c(10, 5, 20, 8),
    SampleB = c(15, 3, 25, 6),
    SampleC = c(12, 4, 22, 7)
  )

  create_raw_madc(madc_path, dup_clone_ids, dup_allele_ids, sample_names, dup_counts,
                  include_summary_cols = TRUE)
  create_marker_file(marker_path,
                     c("marker1", "marker2"),
                     c("Chr01", "Chr02"),
                     c(100, 200))

  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path)

  # All 4 rows should be present

  expect_equal(nrow(result), 4)

  # CloneIDs should both be updated (duplicates are expected)
  expect_equal(unname(result$CloneID), c("Chr01_000000100", "Chr01_000000100", "Chr02_000000200", "Chr02_000000200"))

  # AlleleIDs should have updated suffixes: Ref -> Ref_0001, Alt -> Alt_0002
  expect_equal(unname(result$AlleleID),
               c("Chr01_000000100|Ref_0001", "Chr01_000000100|Alt_0002", "Chr02_000000200|Ref_0001", "Chr02_000000200|Alt_0002"))

  # Sample data should be intact
  expect_equal(result$SampleA, c(10, 5, 20, 8))
})


# ── 2. Preprocessed MADC (no filler rows) triggers error ────────────────────

test_that("fixMADC errors when MADC appears preprocessed (no filler rows)", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_preprocessed_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix)
  create_marker_file(marker_path, clone_ids,
                     c("Chr01", "Chr01", "Chr02", "Chr02"),
                     c(100, 200, 300, 400))

  expect_error(
    fixMADC(madc.file = madc_path, marker.file = marker_path),
    "already use fixed allele IDs"
  )

})


# ── 3. Manual n.summary.columns removal ─────────────────────────────────────

test_that("n.summary.columns removes the correct number of extra columns", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix,
                  include_summary_cols = TRUE)
  create_marker_file(marker_path, clone_ids,
                     c("Chr01", "Chr01", "Chr02", "Chr02"),
                     c(100, 200, 300, 400))

  # 13 summary columns start at col 4, so columns 4:16
  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path,
                    n.summary.columns = 13)

  expect_s3_class(result, "data.frame")
  expect_false("CallRate" %in% colnames(result))
  expect_false("ClusterConsensusSequence" %in% colnames(result))
  expect_true(all(sample_names %in% colnames(result)))
})


# ── 4. Marker file validation — duplicate CloneIDs ──────────────────────────

test_that("fixMADC errors on duplicate marker IDs in column 1", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix)
  # Duplicate first marker
  create_marker_file(marker_path,
                     c("marker1", "marker1", "marker2", "marker3"),
                     c("Chr01", "Chr01", "Chr02", "Chr02"),
                     c(100, 200, 300, 400))

  expect_error(
    fixMADC(madc.file = madc_path, marker.file = marker_path),
    "duplicate marker IDs"
  )
})


# ── 5. Marker file validation — duplicate Chr_Pos ───────────────────────────

test_that("fixMADC errors on duplicate Chr_Pos combinations", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix)
  create_marker_file(marker_path, clone_ids,
                     c("Chr01", "Chr01", "Chr02", "Chr02"),
                     c(100, 100, 300, 400))  # Chr01_100 appears twice

  expect_error(
    fixMADC(madc.file = madc_path, marker.file = marker_path),
    "duplicate Chr and Pos"
  )
})


# ── 6. Marker file validation — special chars in chromosome ─────────────────

test_that("fixMADC errors on special characters in chromosome column", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids[1:2], allele_ids[1:2], sample_names,
                  counts_matrix[1:2, ])
  create_marker_file(marker_path,
                     clone_ids[1:2],
                     c("Chr_01", "Chr02"),  # underscore in chr name
                     c(100, 200))

  expect_error(
    fixMADC(madc.file = madc_path, marker.file = marker_path),
    "Special characters.*chromosome"
  )
})


# ── 7. Marker file validation — special chars in position ───────────────────

test_that("fixMADC errors on special characters in position column", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids[1:2], allele_ids[1:2], sample_names,
                  counts_matrix[1:2, ])
  create_marker_file(marker_path,
                     clone_ids[1:2],
                     c("Chr01", "Chr02"),
                     c("1.5", "200"))  # decimal in position

  expect_error(
    fixMADC(madc.file = madc_path, marker.file = marker_path),
    "Special characters.*position|position.*must be numeric"
  )
})


# ── 8. Marker file validation — non-numeric position ────────────────────────

test_that("fixMADC errors when position column is not numeric", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids[1:2], allele_ids[1:2], sample_names,
                  counts_matrix[1:2, ])
  create_marker_file(marker_path,
                     clone_ids[1:2],
                     c("Chr01", "Chr02"),
                     c("abc", "def"))

  expect_error(
    fixMADC(madc.file = madc_path, marker.file = marker_path),
    "must be numeric"
  )
})


# ── 9. Missing CloneIDs in marker file triggers warning + removal ───────────

test_that("fixMADC warns and removes CloneIDs not in marker file", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix)
  # Marker file only has 3 of the 4 markers
  create_marker_file(marker_path,
                     clone_ids[1:3],
                     c("Chr01", "Chr01", "Chr02"),
                     c(100, 200, 300))

  expect_warning(
    result <- fixMADC(madc.file = madc_path, marker.file = marker_path),
    "not found in the marker file"
  )

  # Only 3 rows should remain
  expect_equal(nrow(result), 3)
  expect_false("marker4" %in% result$CloneID)
})


# ── 10. Output file is written when output.file is provided ─────────────────

test_that("fixMADC writes CSV when output.file is specified", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")
  output_path <- tempfile()  # function appends ".csv"

  create_raw_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix)
  create_marker_file(marker_path, clone_ids,
                     c("Chr01", "Chr01", "Chr02", "Chr02"),
                     c(100, 200, 300, 400))

  fixMADC(madc.file = madc_path,
          marker.file = marker_path,
          output.file = output_path)

  written_file <- paste0(output_path, "_fixedID.csv")
  expect_true(file.exists(written_file))

  saved <- read.csv(written_file)
  expect_equal(saved$CloneID, c("Chr01_000000100", "Chr01_000000200", "Chr02_000000300", "Chr02_000000400"))
})


# ── 11. AlleleID with Alt suffix is handled ─────────────────────────────────

test_that("fixMADC correctly handles |Alt suffixes in AlleleID", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  alt_allele_ids <- paste0(clone_ids, "|Alt")

  create_raw_madc(madc_path, clone_ids, alt_allele_ids, sample_names, counts_matrix)
  create_marker_file(marker_path, clone_ids,
                     c("Chr01", "Chr01", "Chr02", "Chr02"),
                     c(100, 200, 300, 400))

  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path)

  expect_true(all(grepl("\\|Alt_0002$", result$AlleleID)))
  expect_equal(unname(result$AlleleID[1]), "Chr01_000000100|Alt_0002")
})


# ── 11b. AlleleID with |RefMatch and |AltMatch suffixes ─────────────────────

test_that("fixMADC correctly handles |RefMatch and |AltMatch suffixes", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  match_allele_ids <- c("marker1|RefMatch", "marker2|AltMatch",
                        "marker3|RefMatch", "marker4|AltMatch")

  create_raw_madc(madc_path, clone_ids, match_allele_ids, sample_names, counts_matrix)
  create_marker_file(marker_path, clone_ids,
                     c("Chr01", "Chr01", "Chr02", "Chr02"),
                     c(100, 200, 300, 400))

  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path)

  expect_equal(unname(result$AlleleID),
               c("Chr01_000000100|RefMatch_0001", "Chr01_000000200|AltMatch_0001",
                 "Chr02_000000300|RefMatch_0001", "Chr02_000000400|AltMatch_0001"))
})


# ── 11c. Duplicate RefMatch/AltMatch suffixes are numbered ──────────────────

test_that("fixMADC appends _0001, _0002 etc. to duplicate RefMatch/AltMatch within a CloneID", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  # One marker with Ref, Alt, 2x RefMatch, 3x AltMatch = 7 rows
  dup_clone_ids  <- rep("marker1", 7)
  dup_allele_ids <- c("marker1|Ref", "marker1|Alt",
                      "marker1|RefMatch", "marker1|RefMatch",
                      "marker1|AltMatch", "marker1|AltMatch", "marker1|AltMatch")
  dup_counts <- data.frame(
    SampleA = c(10, 5, 8, 12, 3, 6, 2),
    SampleB = c(15, 3, 9, 11, 4, 7, 1),
    SampleC = c(12, 4, 7, 13, 2, 5, 3)
  )

  create_raw_madc(madc_path, dup_clone_ids, dup_allele_ids, sample_names, dup_counts)
  create_marker_file(marker_path,
                     c("marker1"),
                     c("Chr01"),
                     c(100))

  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path)

  expect_equal(nrow(result), 7)

  # Ref always gets _0001, Alt always gets _0002
  expect_equal(unname(result$AlleleID[1]), "Chr01_000000100|Ref_0001")
  expect_equal(unname(result$AlleleID[2]), "Chr01_000000100|Alt_0002")

  # Duplicate RefMatch rows should be numbered _0001, _0002
  expect_equal(unname(result$AlleleID[3]), "Chr01_000000100|RefMatch_0001")
  expect_equal(unname(result$AlleleID[4]), "Chr01_000000100|RefMatch_0002")

  # Duplicate AltMatch rows should be numbered _0001, _0002, _0003
  expect_equal(unname(result$AlleleID[5]), "Chr01_000000100|AltMatch_0001")
  expect_equal(unname(result$AlleleID[6]), "Chr01_000000100|AltMatch_0002")
  expect_equal(unname(result$AlleleID[7]), "Chr01_000000100|AltMatch_0003")
})


# ── 11d. All suffix types are numbered even when unique per CloneID ───────────

test_that("fixMADC numbers all suffix types including single RefMatch/AltMatch per CloneID", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  # Two markers, each with Ref, Alt, one RefMatch, one AltMatch
  mixed_clone_ids  <- c("marker1", "marker1", "marker1", "marker1",
                        "marker2", "marker2", "marker2", "marker2")
  mixed_allele_ids <- c("marker1|Ref", "marker1|Alt", "marker1|RefMatch", "marker1|AltMatch",
                        "marker2|Ref", "marker2|Alt", "marker2|RefMatch", "marker2|AltMatch")
  mixed_counts <- data.frame(
    SampleA = c(10, 5, 8, 3, 20, 8, 12, 6),
    SampleB = c(15, 3, 9, 4, 25, 6, 11, 7),
    SampleC = c(12, 4, 7, 2, 22, 7, 13, 5)
  )

  create_raw_madc(madc_path, mixed_clone_ids, mixed_allele_ids, sample_names, mixed_counts)
  create_marker_file(marker_path,
                     c("marker1", "marker2"),
                     c("Chr01", "Chr02"),
                     c(100, 200))

  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path)

  # All suffixes are always numbered: Ref_0001, Alt_0002, RefMatch_0001, AltMatch_0001
  expect_equal(unname(result$AlleleID),
               c("Chr01_000000100|Ref_0001", "Chr01_000000100|Alt_0002", "Chr01_000000100|RefMatch_0001", "Chr01_000000100|AltMatch_0001",
                 "Chr02_000000200|Ref_0001", "Chr02_000000200|Alt_0002", "Chr02_000000200|RefMatch_0001", "Chr02_000000200|AltMatch_0001"))
})


# ── 12. Whitespace in marker file is trimmed ────────────────────────────────

test_that("fixMADC trims whitespace in marker file columns", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids[1:2], allele_ids[1:2], sample_names,
                  counts_matrix[1:2, ])

  # Write marker file with extra whitespace
  df <- data.frame(
    CloneID = c("  marker1  ", " marker2"),
    Chr = c("Chr01 ", " Chr02"),
    Pos = c(" 100", "200 "),
    stringsAsFactors = FALSE
  )
  write.csv(df, marker_path, row.names = FALSE)

  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path)

  expect_equal(unname(result$CloneID), c("Chr01_000000100", "Chr02_000000200"))
})


# ── 13. Returns NULL (invisibly) when writing to file ────────────────────────

test_that("fixMADC returns NULL when output.file is provided", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")
  output_path <- tempfile()

  create_raw_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix)
  create_marker_file(marker_path, clone_ids,
                     c("Chr01", "Chr01", "Chr02", "Chr02"),
                     c(100, 200, 300, 400))

  result <- fixMADC(madc.file = madc_path,
                    marker.file = marker_path,
                    output.file = output_path)

  expect_null(result)
})


# ── 14. Empty marker file (no matching IDs) removes all rows ────────────────

test_that("fixMADC warns and returns empty df when no CloneIDs match", {
  madc_path   <- tempfile(fileext = ".csv")
  marker_path <- tempfile(fileext = ".csv")

  create_raw_madc(madc_path, clone_ids, allele_ids, sample_names, counts_matrix)
  create_marker_file(marker_path,
                     c("noMatch1", "noMatch2"),
                     c("Chr01", "Chr02"),
                     c(100, 200))

  expect_warning(
    result <- fixMADC(madc.file = madc_path, marker.file = marker_path),
    "not found in the marker file"
  )

  expect_equal(nrow(result), 0)
})
