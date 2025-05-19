context("Imputation Concordance")


test_that("test imputation",{
  #Input variables
  ped_file <- system.file("check_ped_test.txt", package="BIGr")

  #Calculations
  output.list <- check_ped(ped_file,
                           seed = 101919,
                           verbose = FALSE)

  #Check
  df_length <- length(output.list)
  messy_parents <- output.list$messy_parents
  missing_parents <- output.list$missing_parents

  expect_true(df_length == 4)
  expect_true(all(messy_parents$id == c("grandfather2","grandfather3")))
  expect_true(nrow(missing_parents) == 13)

})
