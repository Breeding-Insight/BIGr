context("MADC summary")

test_that("madc_summary returns the expected tables", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
  res <- suppressWarnings(suppressMessages(madc_summary(madc_file, min.depth = 10)))

  expect_type(res, "list")
  for (nm in c("mhap_freq", "per_sample", "per_marker", "missingness",
               "depth_sweep", "marker_qc", "overall"))
    expect_s3_class(res[[nm]], "data.frame")

  # mhap_freq: rows 2..10 + ">10" = 10 rows; example distribution 15/3/3 at 2/3/4
  expect_equal(nrow(res$mhap_freq), 10)
  expect_equal(res$mhap_freq$n_markers[res$mhap_freq$n_mhaps == "2"], 15)
  expect_equal(res$mhap_freq$n_markers[res$mhap_freq$n_mhaps == "3"], 3)
  expect_equal(res$mhap_freq$n_markers[res$mhap_freq$n_mhaps == "4"], 3)
  expect_equal(sum(res$mhap_freq$n_markers), 21)

  expect_equal(nrow(res$per_sample), 10)   # 10 samples
  expect_equal(nrow(res$per_marker), 21)   # 21 CloneIDs
  expect_equal(nrow(res$missingness), 12)  # 12 thresholds
  expect_equal(nrow(res$overall), 1)

  # per_sample columns and value ranges (robust stats + read-weighted off-target)
  expect_true(all(c("sample", "depth", "median_depth", "depth_mad", "total_reads",
                    "missing_rate", "call_rate", "n_mhaps_observed",
                    "offtarget_fraction", "offtarget_fraction_weighted")
                  %in% names(res$per_sample)))
  expect_true(all(res$per_sample$missing_rate >= 0 & res$per_sample$missing_rate <= 1))
  expect_equal(res$per_sample$call_rate, 1 - res$per_sample$missing_rate)

  # per_marker splits defined vs observed mHap counts, defined >= observed
  expect_true(all(c("n_mhaps_defined", "n_mhaps_observed") %in% names(res$per_marker)))
  expect_true(all(res$per_marker$n_mhaps_defined >= res$per_marker$n_mhaps_observed))

  # depth_sweep: one row per threshold, monotone non-increasing cells_passing
  expect_equal(nrow(res$depth_sweep), 6)
  expect_true(all(diff(res$depth_sweep$cells_passing) <= 0))

  # marker_qc: one row per marker, valid status vocabulary
  expect_equal(nrow(res$marker_qc), 21)
  expect_true(all(res$marker_qc$qc_status %in% c("PASS", "FLAG", "FAIL")))
  # overall carries both the mean-locus and read-weighted off-target
  expect_true(all(c("mean_offtarget_fraction", "read_weighted_offtarget_fraction",
                    "median_depth", "overall_call_rate") %in% names(res$overall)))
})

test_that("marker_qc criteria flag/fail as configured", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
  # disable the depth FAIL, enable an mHap FLAG that the fixture triggers
  res <- suppressWarnings(suppressMessages(
    madc_summary(madc_file, min.locus.depth = NULL, max.mhaps.per.loci = 3)))
  flagged <- res$marker_qc$qc_status == "FLAG"
  expect_true(any(flagged))
  expect_true(all(grepl("paralog", res$marker_qc$qc_reason[flagged])))
  # with all criteria off, every marker passes
  res0 <- suppressWarnings(suppressMessages(
    madc_summary(madc_file, min.locus.depth = NULL)))
  expect_true(all(res0$marker_qc$qc_status == "PASS"))
})

test_that("target.only counts only Ref/Alt", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
  res <- suppressWarnings(suppressMessages(madc_summary(madc_file, target.only = TRUE)))
  # 20 markers keep both Ref and Alt (one marker has only a Ref target row)
  expect_equal(res$mhap_freq$n_markers[res$mhap_freq$n_mhaps == "2"], 20)
  expect_equal(res$overall$mean_offtarget_fraction, 0)
})

test_that("metadata adds a by_category table", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
  meta <- data.frame(sample  = paste0("Sample_", 1:10),
                     species = rep(c("A", "B"), 5),
                     stringsAsFactors = FALSE)
  res <- suppressWarnings(suppressMessages(
    madc_summary(madc_file, metadata = meta, group.col = "species")))
  expect_s3_class(res$by_category, "data.frame")
  expect_equal(nrow(res$by_category), 2)
  expect_equal(sort(res$by_category$group), c("A", "B"))
  expect_equal(sum(res$by_category$n_samples), 10)
})

test_that("output.file writes CSVs", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
  out <- tempfile()
  suppressWarnings(suppressMessages(madc_summary(madc_file, output.file = out)))
  expect_true(file.exists(paste0(out, "_mhap_freq.csv")))
  expect_true(file.exists(paste0(out, "_per_sample.csv")))
  expect_true(file.exists(paste0(out, "_missingness.csv")))
  expect_true(file.exists(paste0(out, "_depth_sweep.csv")))
  expect_true(file.exists(paste0(out, "_marker_qc.csv")))
})

test_that("output.format = 'xlsx' writes one workbook with a sheet per table", {
  skip_if_not_installed("writexl")
  skip_if_not_installed("readxl")
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
  out <- tempfile()
  suppressWarnings(suppressMessages(
    madc_summary(madc_file, output.file = out, output.format = "xlsx")))
  xl <- paste0(out, ".xlsx")
  expect_true(file.exists(xl))
  expect_true(all(c("mhap_freq", "per_sample", "per_marker", "missingness",
                    "depth_sweep", "marker_qc", "overall") %in%
                  readxl::excel_sheets(xl)))
})

test_that("raw MADC is rejected", {
  raw <- system.file("iris_DArT_MADC.csv", package = "BIGr")
  expect_error(suppressMessages(madc_summary(raw)), "HapApp")
})
