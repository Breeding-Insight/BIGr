context("Merge MADCs")


test_that("test merge madc",{
  #Input variables
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  madc2_file <- system.file("example_MADC_to_merge.csv", package="BIGr")

  #Calculations
  temp <- tempfile(fileext = ".csv")
  temp2 <- tempfile(fileext = ".csv")

  merge_MADCs(madc_list = list(madc_file,madc2_file),
              out_madc=temp, run_ids=NULL)

  merge_MADCs(madc_list = list(madc_file,madc_file), out_madc=temp2, run_ids=NULL)

  merged_madc <- data.frame(read_csv(temp))
  merged_madc2 <- data.frame(read_csv(temp2))

  #Check
  count_sum <- sum(as.matrix(merged_madc[,-c(1,2,3)]))
  df_dim <- dim(merged_madc)


  expect_true(all(df_dim == c("61","23")))
  expect_true(count_sum == 86845)
  expect_error(merge_MADCs(madc_list = NULL,out_madc=temp, run_ids=NULL))
  expect_error(merge_MADCs(madc_list = list(madc_file,madc2_file), out_madc=NULL, run_ids=NULL))
  expect_error(merge_MADCs(madc_list = list(madc_file,madc2_file), out_madc=temp, run_ids="one"))
  expect_true(all(merged_madc2[,4:13] == merged_madc2[,14:23]))

  rm(count_sum,merged_madc,merged_madc2,df_dim)

})
