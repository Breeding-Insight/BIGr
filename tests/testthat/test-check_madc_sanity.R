test_that("check madc",{
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  report <- read.csv(madc_file, check.names = FALSE)

  check_madc_sanity(report)


})


