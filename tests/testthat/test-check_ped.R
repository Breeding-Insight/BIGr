context("Imputation Concordance")


test_that("test imputation",{
  #Input variables
  ped_file <- system.file("check_ped_test.txt", package="BIGr")
  temp_ped <- tempfile("check_ped_test_", fileext = ".txt")
  file.copy(ped_file, temp_ped)

  #Calculations
  output.list <- check_ped(temp_ped,
                           seed = 101919,
                           verbose = FALSE)

  #Check
  df_length <- length(output.list)
  messy_parents <- output.list$messy_parents
  missing_parents <- output.list$missing_parents
  corrected_pedigree <- output.list$corrected_pedigree
  file_base <- tools::file_path_sans_ext(basename(temp_ped))

  expect_true(df_length == 6)
  expect_true(all(messy_parents$id == c("grandfather2","grandfather3")))
  expect_true(nrow(missing_parents) == 13)
  expect_s3_class(corrected_pedigree, "data.frame")
  expect_true(all(missing_parents$id %in% corrected_pedigree$id))
  expect_false(exists(paste0(file_base, "_corrected"), envir = .GlobalEnv, inherits = FALSE))
  expect_false(exists(paste0(file_base, "_report"), envir = .GlobalEnv, inherits = FALSE))

})
