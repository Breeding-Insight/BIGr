context("Updog to VCF")


test_that("test updog conversion",{
  #Input variables
  updog_file <- system.file("iris_updog_multidog.RData", package="BIGr")

  temp_file <- tempfile()
  load(updog_file)

  # Convert updog to VCF
  updog2vcf(
    multidog.object = mout,
    output.file = temp_file,
    updog_version = packageVersion("updog"),
    compress = TRUE
  )

  vcf_result <- read.vcfR(paste0(temp_file,".vcf.gz"), verbose = FALSE)

  DP <- sum(as.numeric(extract.gt(vcf_result, "DP")))

  expect_equal(DP, 23618990)

  MPP <- sum(as.numeric(extract.gt(vcf_result, "MPP")))

  expect_equal(MPP, 74519.94, tolerance = 0.01)

})
