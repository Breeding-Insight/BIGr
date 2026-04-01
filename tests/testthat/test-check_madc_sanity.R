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


