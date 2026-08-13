context("MADC summary")

test_that("madc_summary returns the expected tables", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
  res <- suppressWarnings(suppressMessages(madc_summary(madc_file, min.depth = 10)))

  expect_type(res, "list")
  for (nm in c("mhap_freq", "per_sample", "per_marker", "missingness", "overall"))
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

  # per_sample columns and value ranges
  expect_true(all(c("sample", "depth", "missing_rate", "n_mhaps",
                    "offtarget_fraction") %in% names(res$per_sample)))
  expect_true(all(res$per_sample$missing_rate >= 0 & res$per_sample$missing_rate <= 1))
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
})

test_that("raw MADC is rejected", {
  raw <- system.file("iris_DArT_MADC.csv", package = "BIGr")
  expect_error(suppressMessages(madc_summary(raw)), "HapApp")
})
