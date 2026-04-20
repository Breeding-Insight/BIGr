test_that("check madc",{

  github_path <- "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/long_seq/test_madcs/"
  names <- c("Columns", "FixAlleleIDs", "IUPACcodes", "LowerCase", "Indels", "ChromPos", "allNAcol", "allNArow", "RefAltSeqs", "OtherAlleles")

  # raw madc
  report <- read.csv(paste0(github_path,"/alfalfa_madc_raw.csv"))

  res <- check_madc_sanity(report)
  exp <- c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
  names(exp) <- names
  expect_equal(res$checks, exp)

  # test lower case
  report <- read.csv(paste0(github_path,"/alfalfa_lowercase.csv"))

  res <- check_madc_sanity(report)
  exp <- c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
  names(exp) <- names
  expect_equal(res$checks, exp)

  # test IUPAC
  report <- read.csv(paste0(github_path,"/alfalfa_IUPAC.csv"))

  res <- check_madc_sanity(report)
  exp <- c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,TRUE, FALSE)
  names(exp) <- names
  expect_equal(res$checks, exp)

  # clean alfalfa madc (fixed allele IDs, Chr_Pos CloneIDs, no issues)
  report <- read.csv(paste0(github_path,"/alfalfa_madc.csv"))

  res <- check_madc_sanity(report)
  exp <- c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,TRUE, FALSE)
  names(exp) <- names
  expect_equal(res$checks, exp)

  # potato indel madc (ChromPos FALSE because IDs are not Chr_Pos)
  report <- read.csv(paste0(github_path,"/potato_indel_madc.csv"))

  res <- check_madc_sanity(report)
  exp <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,TRUE, FALSE)
  names(exp) <- names
  expect_equal(res$checks, exp)

  # potato indel IUPAC (IUPAC codes + lowercase + indels)
  report <- read.csv(paste0(github_path,"/potato_indel_IUPAC.csv"))

  res <- check_madc_sanity(report)
  exp <- c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE,TRUE, FALSE)
  names(exp) <- names
  expect_equal(res$checks, exp)

  # potato indel lowercase (lowercase + indels)
  report <- read.csv(paste0(github_path,"/potato_indel_lowercase.csv"))

  res <- check_madc_sanity(report)
  exp <- c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE,TRUE, FALSE)
  names(exp) <- names
  expect_equal(res$checks, exp)

  # potato more indels madc ChromPosFALSE
  report <- read.csv(paste0(github_path,"/potato_more_indels_madc_ChromPosFALSE.csv"))

  res <- check_madc_sanity(report)
  exp <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,TRUE, FALSE)
  names(exp) <- names
  expect_equal(res$checks, exp)
})

test_that("check_botloci remaps using Marker_ID", {
  botloci <- data.frame(V1 = c("1_0001", "2_0002"))
  report <- data.frame(
    CloneID = c("ProbeA_0001", "ProbeB_0002"),
    AlleleID = c("ProbeA_0001|Ref_0001", "ProbeB_0002|Ref_0001"),
    AlleleSequence = c("A", "T"),
    check.names = FALSE
  )
  mi_df <- data.frame(
    Marker_ID = c("ProbeA_0001", "ProbeB_0002"),
    Chr = c("1", "2"),
    Pos = c(1, 2)
  )

  res <- check_botloci(botloci, report, ChromPos = FALSE, mi_df = mi_df, verbose = FALSE)

  expect_equal(res[[2]]$CloneID, botloci$V1)
  expect_equal(res[[3]]$CloneID, botloci$V1)
})

test_that("check_botloci resolves Marker_ID before padding report CloneIDs", {
  botloci <- data.frame(V1 = "1_000000123")
  report <- data.frame(
    CloneID = "1_123",
    AlleleID = "1_123|Ref_0001",
    AlleleSequence = "A",
    check.names = FALSE
  )
  mi_df <- data.frame(
    Marker_ID = "1_123",
    Chr = "1",
    Pos = 123
  )

  res <- check_botloci(botloci, report, ChromPos = TRUE, mi_df = mi_df, verbose = FALSE)

  expect_equal(res[[2]]$CloneID, botloci$V1)
  expect_equal(res[[3]]$CloneID, botloci$V1)
})

test_that("pick_markers_info_id_col scores distinct markers not allele rows", {
  mi_df <- data.frame(
    CloneID = c("m1", "m2"),
    Marker_ID = c("m1", "m3")
  )
  query_ids <- c("m1", "m1", "m1", "m2")

  expect_equal(pick_markers_info_id_col(mi_df, query_ids), "CloneID")
})

test_that("check_madc_sanity returns FALSE for malformed CloneID positions", {
  report <- data.frame(
    CloneID = c("Chr_abc", "Chr_abc"),
    AlleleID = c("Chr_abc|Ref_0001", "Chr_abc|Alt_0002"),
    AlleleSequence = c("A", "T"),
    check.names = FALSE
  )

  res <- check_madc_sanity(report)

  expect_false(is.na(res$checks["ChromPos"]))
  expect_false(res$checks["ChromPos"])
})

test_that("check_botloci errors if widening MADC padding still does not match", {
  botloci <- data.frame(V1 = "1_0002")
  report <- data.frame(
    CloneID = "1_1",
    AlleleID = "1_1|Ref_0001",
    AlleleSequence = "A",
    check.names = FALSE
  )

  expect_error(
    check_botloci(botloci, report, ChromPos = TRUE, verbose = FALSE),
    "After matching padding, botloci markers still not found in MADC file. Check marker IDs."
  )
})

test_that("check_botloci keeps AlleleID synchronized after CloneID remap", {
  botloci <- data.frame(V1 = "1_0001")
  report <- data.frame(
    CloneID = "ProbeA_0001",
    AlleleID = "ProbeA_0001|Ref_0001",
    AlleleSequence = "A",
    check.names = FALSE
  )
  mi_df <- data.frame(
    Marker_ID = "ProbeA_0001",
    Chr = "1",
    Pos = 1
  )

  res <- check_botloci(botloci, report, ChromPos = TRUE, mi_df = mi_df, verbose = FALSE)

  expect_equal(res[[2]]$CloneID, "1_0001")
  expect_equal(res[[2]]$AlleleID, "1_0001|Ref_0001")
})
