context("Updog to VCF")


test_that("test updog conversion",{
  #Input variables
  load(system.file("extdata", "iris-multidog.rdata", package = "BIGr"))

  temp_file <- tempfile()

  # Convert updog to VCF
  updog2vcf(
    multidog.object = mout,
    output.file = temp_file,
    updog_version = "0.0.0",
    compress = TRUE
  )

  vcf_result <- read.vcfR(paste0(temp_file,".vcf.gz"), verbose = FALSE)

  DP <- sum(as.numeric(extract.gt(vcf_result, "DP")))

  expect_equal(DP, 11645)

  MPP <- sum(as.numeric(extract.gt(vcf_result, "MPP")))

  expect_equal(MPP, 39.99, tolerance = 0.01)

})
