context("Obtain Read Counts")


test_that("test madc read count extraction",{
  #Input variables
  input <- diversity_items <- list()
  input$madc_file$datapath <- system.file("iris_DArT_MADC.csv", package = "BIGr")

  madc <- input$madc_file$datapath

  #Retrieve read counts from MADC file for target markers only
  counts_matrices <- get_countsMADC(madc)

  #Checks
  expect_true(all(dim(counts_matrices[[1]]) == dim(counts_matrices[[2]])))
  expect_true(all(names(counts_matrices) == c("ref_matrix", "size_matrix")))
  expect_true(sum(counts_matrices[[1]]) == 15568260)
  expect_true(sum(counts_matrices[[2]]) == 23618990)

})
