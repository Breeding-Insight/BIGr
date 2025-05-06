context("BreedTools")


test_that("test breedtools poly",{
  #Input variables
  ref_file <- system.file("test_ref.txt", package="BIGr")
  val_file <- system.file("test_test.txt", package="BIGr")
  ref_ids <- system.file("ref_ids.txt", package="BIGr")

  #import files
  reference = read.table(ref_file, header = T, row.names = 1, sep = "\t")
  validation = read.table(val_file, header = T, row.names = 1, sep = "\t")
  reference_ids = read.table(ref_ids, header = T, sep = "\t")

  #Calculations
  ref_ids = lapply(as.list(reference_ids),as.character)

  freq = allele_freq_poly(reference, ref_ids, ploidy = 4)

  prediction = as.data.frame(solve_composition_poly(validation,freq, ploidy = 4))

  #Check
  freq_mean <- round(mean(as.numeric(freq)),6)
  pred_mean <- round(mean(as.numeric(prediction$R2)),6)


  expect_equal(freq_mean, 0.888889, tolerance = 0.01)
  expect_equal(pred_mean, 0.841454, tolerance = 0.01)
  expect_true(nrow(prediction) == 175)

})
