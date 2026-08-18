context("Filter MADC")


test_that("test filter madc",{
  #Input variables
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")

  #Calculations
  temp <- tempfile()

  #Reference values from the raw input, used to check filter invariants
  raw_in <- read.csv(madc_file, check.names = FALSE)
  in_counts <- table(raw_in$CloneID)
  n_refalt_in <- sum(grepl("\\|(Ref|Alt)_", raw_in$AlleleID))

  # Filtering (target only)
  filtered_df <- suppressMessages(filterMADC(madc_file,
                         min.locus.depth = NULL,
                         max.locus.depth = NULL,
                         max.mhaps.per.loci = NULL,
                         min.reads.per.site = 1,
                         min.ind.with.reads = NULL,
                         target.only = TRUE,
                         output.file = NULL))


  #Test that a valid output was provided
  expect_equal(nrow(filtered_df), 41)
  #Check that it is a dataframe
  expect_true(is.data.frame(filtered_df))

  # Checking for no filtering
  filtered_df <- suppressMessages(filterMADC(madc_file,
                            min.locus.depth = NULL,
                            max.locus.depth = NULL,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            output.file = NULL))

  expect_equal(nrow(filtered_df), 51)
  expect_equal(sum(filtered_df[,-c(1:3)]), 53952)
  expect_true(all(names(filtered_df[1:3]) == c("AlleleID", "CloneID", "AlleleSequence")))

  #Checking for min.locus.depth filtering (drops whole low-depth loci)
  filtered_df <- suppressMessages(filterMADC(madc_file,
                            min.locus.depth = 50,
                            max.locus.depth = NULL,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            output.file = NULL))

  expect_equal(nrow(filtered_df), 37)
  expect_equal(ncol(filtered_df), 13)
  #Whole-locus removal: every retained locus keeps ALL of its original mhap rows (no un-pairing)
  out_counts <- table(filtered_df$CloneID)
  expect_true(all(out_counts == in_counts[names(out_counts)]))

  #Checking for max.locus.depth filtering (drops whole high-depth loci)
  filtered_df <- suppressMessages(filterMADC(madc_file,
                            min.locus.depth = NULL,
                            max.locus.depth = 50,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            output.file = NULL))

  expect_equal(nrow(filtered_df), 14)
  expect_equal(ncol(filtered_df), 13)
  expect_equal(sum(filtered_df[,-c(1:3)]), 1122)

  #Remove max mhaps
  filtered_df <- suppressMessages(filterMADC(madc_file,
                            min.locus.depth = NULL,
                            max.locus.depth = NULL,
                            max.mhaps.per.loci = 3,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            output.file = NULL))

  expect_equal(nrow(filtered_df), 44)
  expect_equal(ncol(filtered_df), 13)

  #Presence filter (min ind with reads) - target Ref/Alt are protected
  filtered_df <- suppressMessages(filterMADC(madc_file,
                            min.locus.depth = NULL,
                            max.locus.depth = NULL,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 10,
                            min.ind.with.reads = 10,
                            target.only = FALSE,
                            output.file = NULL))

  expect_equal(nrow(filtered_df), 41)
  expect_equal(ncol(filtered_df), 13)
  expect_equal(sum(filtered_df[,-c(1:3)]), 47985)
  #Every input target Ref/Alt row is retained (only non-target mhaps can be pruned)
  expect_equal(sum(grepl("\\|(Ref|Alt)_", filtered_df$AlleleID)), n_refalt_in)

  #Check that the output file is created
  suppressMessages(filterMADC(madc_file,
                            min.locus.depth = NULL,
                            max.locus.depth = NULL,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            output.file = temp))

  expect_true(file.exists(paste0(temp,".csv")))

  #Now checking that all paramaters can work together
  filtered_df <- suppressMessages(filterMADC(madc_file,
                            min.locus.depth = 10,
                            max.locus.depth = 300,
                            max.mhaps.per.loci = 3,
                            min.reads.per.site = 10,
                            min.ind.with.reads = 5,
                            target.only = FALSE,
                            output.file = NULL))

  expect_equal(nrow(filtered_df), 29)
  expect_equal(ncol(filtered_df), 13)
  expect_equal(sum(filtered_df[,-c(1:3)]), 15375)


  #Reject a raw (non fixed allele ID) MADC file
  raw_file <- system.file("iris_DArT_MADC.csv", package="BIGr")
  expect_error(suppressMessages(filterMADC(raw_file)), "HapApp")


  #Fixed allele ID fixture containing |Other alleles to exercise the new
  #Ref/Alt-only retention in max.mhaps.per.loci and target.only
  tmp_other <- tempfile(fileext = ".csv")
  other_df <- data.frame(
    AlleleID = c("Chr1_001|Ref_0001", "Chr1_001|Alt_0002",
                 "Chr1_001|RefMatch_0001", "Chr1_001|AltMatch_0001",
                 "Chr1_001|Other_0001"),
    CloneID = "Chr1_001",
    AlleleSequence = c("ACGT", "ACGA", "ACGT", "ACGA", "TTTT"),
    Sample_1 = c(10, 8, 5, 4, 3),
    Sample_2 = c(12, 9, 6, 5, 2),
    check.names = FALSE
  )
  write.csv(other_df, tmp_other, row.names = FALSE)

  # max.mhaps.per.loci: 5 mhaps at the locus > 3 -> keep only |Ref/|Alt (drop Match AND Other)
  out_mhap <- suppressMessages(filterMADC(tmp_other, max.mhaps.per.loci = 3))
  expect_equal(nrow(out_mhap), 2)
  expect_true(all(grepl("\\|(Ref|Alt)_", out_mhap$AlleleID)))

  # target.only: retains strictly |Ref/|Alt, dropping Match AND Other
  out_target <- suppressMessages(filterMADC(tmp_other, target.only = TRUE))
  expect_equal(nrow(out_target), 2)
  expect_false(any(grepl("\\|Other", out_target$AlleleID)))

})
