context("MADC 2 Gmatrix")


test_that("test madc2gmat",{
  #Input variables
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")

  #Calculations
  temp <- tempfile()

  # Converting to additive relationship matrix
  gmat <- madc2gmat(madc_file,
                    seed = 123,
                    output.file = NULL)

  #When output a file
  madc2gmat(madc_file,
            seed = 123,
            output.file = temp)

  #Test that a valid output was provided
  expect_true(file.exists(paste0(temp, ".csv")))

  #Check
  expect_true(all(dim(gmat) == c("10","10")))
  expect_true(all(row.names(gmat) == row.names(gmat)))
  expect_equal(sum(gmat), -1.480586e-15, tolerance = 1e-16)
  expect_true(is.matrix(gmat), "Output should be a matrix")

  # Read the output file
  output_data <- read.csv(paste0(temp,".csv"), row.names = 1)

  # Test the content of the output file
  expect_true(is.matrix(as.matrix(output_data)), "Data in output file should be a matrix")
  expect_true(all(dim(output_data) == c("10","10")))
  expect_identical(row.names(output_data), colnames(output_data), "Row and column names in output file should be identical")
  expect_equal(sum(output_data), 6.214647e-15, tolerance = 1e-15)
})
