context("MADC plot")

madc_file <- function() system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")

test_that("single plot.type returns a ggplot", {
  p <- suppressWarnings(suppressMessages(madc_plot(madc_file(), plot.type = "missing")))
  expect_s3_class(p, "ggplot")
})

test_that("pca, heatmap, marker build without error", {
  expect_no_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "pca"))))
  expect_no_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "heatmap", fill = "depth"))))
  expect_no_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "heatmap", fill = "mhaps"))))
  # example CloneIDs are Chr_Pos -> marker plot works from CloneID directly
  expect_no_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "marker"))))
})

test_that("balance and depth plots build", {
  expect_s3_class(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "balance"))), "ggplot")
  # ploidy adds expected-ratio guides
  expect_s3_class(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "balance", ploidy = 2))), "ggplot")
  expect_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "balance", ploidy = 0))), "positive")
  expect_s3_class(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "depth"))), "ggplot")
  # both embed in a multi-panel figure
  res <- suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = c("depth", "balance"))))
  expect_true(grid::is.grob(res$panel))
})

test_that("circos plot builds without error", {
  skip_if_not_installed("circlize")
  tmp <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "circos", output.file = tmp))))
  expect_true(file.exists(tmp))
})

test_that("circos must be requested on its own", {
  expect_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = c("circos", "pca")))), "on its own")
})

test_that("circos honors a custom density.window", {
  skip_if_not_installed("circlize")
  tmp <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "circos", density.window = 5e5, output.file = tmp))))
})

test_that("circos honors depth.qc.maxmiss for the combined Depth ring", {
  skip_if_not_installed("circlize")
  tmp <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "circos", depth.qc.maxmiss = 0.3, output.file = tmp))))
  expect_true(file.exists(tmp))
})

test_that("circos accepts custom label.gap, paralog.flag, and track colors", {
  skip_if_not_installed("circlize")
  tmp <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "circos", label.gap = 20, paralog.flag = 5,
              mhap.col = "#1f78b4", depth.col = c("navy", "grey90", "orange"),
              output.file = tmp))))
  expect_true(file.exists(tmp))
})

test_that("heatmap groups by category and facets by chromosome", {
  meta <- data.frame(sample = paste0("Sample_", 1:10),
                     species = rep(c("A", "B"), 5), stringsAsFactors = FALSE)
  p <- suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "heatmap",
              metadata = meta, group.col = "species", facet.chrom = TRUE)))
  expect_s3_class(p, "ggplot")
})

test_that("marker plot works with a markers_info lookup", {
  cids <- unique(read.csv(madc_file(), check.names = FALSE)$CloneID)
  mi <- data.frame(CloneID = cids, Chr = "chr1.1",
                   Pos = as.numeric(sub(".*_", "", cids)),
                   stringsAsFactors = FALSE)
  expect_no_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "marker", markers_info = mi))))
})

test_that("multiple plot.types return list + panel", {
  res <- suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = c("pca", "missing"))))
  expect_type(res, "list")
  expect_s3_class(res$plots$pca, "ggplot")
  expect_true(grid::is.grob(res$panel))
})

test_that("bare default call builds the 4-panel figure (no circos)", {
  res <- suppressWarnings(suppressMessages(madc_plot(madc_file())))
  expect_type(res, "list")
  expect_equal(names(res$plots), c("pca", "marker", "heatmap", "missing"))
  expect_true(grid::is.grob(res$panel))
})

test_that("metadata colors PCA and groups the missing boxplot", {
  meta <- data.frame(sample = paste0("Sample_", 1:10),
                     species = rep(c("A", "B"), 5), stringsAsFactors = FALSE)
  p1 <- suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "pca", metadata = meta, group.col = "species")))
  p2 <- suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "missing", metadata = meta, group.col = "species")))
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("PCA supports shape.col and a custom palette", {
  meta <- data.frame(sample = paste0("Sample_", 1:10),
                     species = rep(c("A", "B"), 5),
                     plate   = rep(c("P1", "P2"), each = 5),
                     stringsAsFactors = FALSE)
  p <- suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "pca", metadata = meta,
              group.col = "species", shape.col = "plate",
              palette = c("#440154", "#21918c", "#fde725"))))
  expect_s3_class(p, "ggplot")
})

test_that("missing boxplot supports miss.sort and horizontal", {
  meta <- data.frame(sample = paste0("Sample_", 1:10),
                     species = rep(c("A", "B"), 5), stringsAsFactors = FALSE)
  for (s in c("none", "asc", "desc")) for (h in c(TRUE, FALSE)) {
    p <- suppressWarnings(suppressMessages(
      madc_plot(madc_file(), plot.type = "missing", metadata = meta,
                group.col = "species", miss.sort = s, horizontal = h)))
    expect_s3_class(p, "ggplot")
  }
  # horizontal also works without a grouping category
  expect_s3_class(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "missing", horizontal = TRUE))), "ggplot")
  expect_error(suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "missing", miss.sort = "bogus"))))
})

test_that("output.file writes an image", {
  out <- tempfile(fileext = ".png")
  suppressWarnings(suppressMessages(
    madc_plot(madc_file(), plot.type = "missing", output.file = out)))
  expect_true(file.exists(out))
})

test_that("raw MADC is rejected", {
  raw <- system.file("iris_DArT_MADC.csv", package = "BIGr")
  expect_error(suppressMessages(madc_plot(raw, plot.type = "pca")), "HapApp")
})
