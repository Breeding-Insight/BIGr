context("Filter MADC")


test_that("test filter madc",{
  #Input variables
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")

  #Calculations
  temp <- tempfile()

  # Filtering (target only)
  filtered_df <- filterMADC(madc_file,
                         min.mean.reads = NULL,
                         max.mean.reads = NULL,
                         max.mhaps.per.loci = NULL,
                         min.reads.per.site = 1,
                         min.ind.with.reads = NULL,
                         target.only = TRUE,
                         n.summary.columns = NULL,
                         output.file = NULL)


  #Test that a valid output was provided
  expect_equal(nrow(filtered_df), 41)
  #Check that it is a dataframe
  expect_true(is.data.frame(filtered_df))

  # Checking for no filtering
  filtered_df <- filterMADC(madc_file,
                            min.mean.reads = NULL,
                            max.mean.reads = NULL,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            n.summary.columns = NULL,
                            output.file = NULL)

  expect_equal(nrow(filtered_df), 51)
  expect_equal(sum(filtered_df[,-c(1:3)]), 53952)
  expect_true(all(names(filtered_df[1:3]) == c("AlleleID", "CloneID", "AlleleSequence")))

  #Checking for min.mean.reads filtering
  filtered_df <- filterMADC(madc_file,
                            min.mean.reads = 10,
                            max.mean.reads = NULL,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            n.summary.columns = NULL,
                            output.file = NULL)

  expect_equal(nrow(filtered_df), 36)
  expect_equal(ncol(filtered_df), 13)

  #Checking for max.mean.reads filtering
  filtered_df <- filterMADC(madc_file,
                            min.mean.reads = NULL,
                            max.mean.reads = 10,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            n.summary.columns = NULL,
                            output.file = NULL)

  expect_equal(nrow(filtered_df), 15)
  expect_equal(ncol(filtered_df), 13)

  #Remove max mhaps
  filtered_df <- filterMADC(madc_file,
                            min.mean.reads = NULL,
                            max.mean.reads = NULL,
                            max.mhaps.per.loci = 3,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            n.summary.columns = NULL,
                            output.file = NULL)

  expect_equal(nrow(filtered_df), 44)
  expect_equal(ncol(filtered_df), 13)

  #Remove min ind with reads
  filtered_df <- filterMADC(madc_file,
                            min.mean.reads = NULL,
                            max.mean.reads = NULL,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 10,
                            min.ind.with.reads = 10,
                            target.only = FALSE,
                            n.summary.columns = NULL,
                            output.file = NULL)

  expect_equal(nrow(filtered_df), 9)
  expect_equal(ncol(filtered_df), 13)
  expect_equal(sum(filtered_df[,-c(1:3)]), 31642)

  #Check that the output file is created
  filterMADC(madc_file,
                            min.mean.reads = NULL,
                            max.mean.reads = NULL,
                            max.mhaps.per.loci = NULL,
                            min.reads.per.site = 1,
                            min.ind.with.reads = NULL,
                            target.only = FALSE,
                            n.summary.columns = NULL,
                            output.file = temp)

  expect_true(file.exists(paste0(temp,".csv")))

  #Check that the plots are created in the console
  #filtered_df <- filterMADC(madc_file,
  #                          min.mean.reads = NULL,
  #                          max.mean.reads = NULL,
  #                          max.mhaps.per.loci = 3,
  #                          min.reads.per.site = 1,
  #                          min.ind.with.reads = NULL,
  #                          target.only = FALSE,
  #                          n.summary.columns = NULL,
  #                          plot.summary = TRUE,
  #                          output.file = NULL)

  #expect_true(is.numeric(dev.cur()))
  #expect_true(dev.cur() > 1)

  #Now checking that all paramaters can work together
  filtered_df <- filterMADC(madc_file,
                            min.mean.reads =1,
                            max.mean.reads = 150,
                            max.mhaps.per.loci = 3,
                            min.reads.per.site = 10,
                            min.ind.with.reads = 10,
                            target.only = FALSE,
                            n.summary.columns = NULL,
                            output.file = NULL)

  expect_equal(nrow(filtered_df), 3)
  expect_equal(ncol(filtered_df), 13)
  expect_equal(sum(filtered_df[,-c(1:3)]), 3960)


})
