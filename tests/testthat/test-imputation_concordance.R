context("Imputation Concordance")


test_that("test imputation",{
  #Input variables
  ignore_file <- system.file("imputation_ignore.txt", package="BIGr")
  ref_file <- system.file("imputation_reference.txt", package="BIGr")
  test_file <- system.file("imputation_test.txt", package="BIGr")

  #import files
  snps = read.table(ignore_file, header = TRUE)
  ref = read.table(ref_file, header = TRUE)
  test = read.table(test_file, header = TRUE)

  #Calculations
  result <- imputation_concordance(ref, test,snps_2_exclude = snps, missing_code =5, output = NULL, verbose = FALSE)

  #Check
  result2 <- sum(as.numeric(gsub("%","",result$Concordance)))

  expect_equal(result2, 1910.51, tolerance = 0.01)
  expect_true(nrow(result) == nrow(ref))

})
