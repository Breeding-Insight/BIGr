context("MADC to VCF")

test_that("test madc conversion",{
  #Input variables
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  bot_file <- system.file("example_SNPs_DArTag-probe-design_f180bp.botloci", package="BIGr")

  #Calculations
  temp <- tempfile(fileext = ".vcf")

  #Convert the dart files to vcf
  madc2vcf_targets(madc_file = madc_file, output.file = temp, get_REF_ALT = FALSE)

  #Test validity of VCF
  vcf <- read.vcfR(temp, verbose = FALSE)

  DP_df <- extract.gt(vcf, element = "DP")
  DP_df_mean <- mean(as.numeric(DP_df))

  #Checks
  expect_s4_class(vcf, "vcfR")
  expect_true(nrow(vcf) == 20)
  expect_true(ncol(vcf@gt) == 11)
  expect_true(all(dim(DP_df) == c(20,10)))
  expect_true(round(DP_df_mean, 3) == 227.22)

  rm(vcf)
  rm(temp)

  # Test with REF_ALT
  temp <- tempfile(fileext = ".vcf")
  madc2vcf_targets(madc_file = madc_file, output.file = temp, get_REF_ALT = TRUE, botloci_file = bot_file)

  #Test validity of VCF
  vcf <- read.vcfR(temp, verbose = FALSE)

  #Checks
  expect_s4_class(vcf, "vcfR")
  expect_true(nrow(vcf) == 20)
  expect_true(ncol(vcf@gt) == 11)
  expect_true(all(dim(DP_df) == c(20,10)))
  expect_true(round(DP_df_mean, 3) == 227.22)
  expect_true(all(vcf@fix[,4] != "."))
  expect_true(all(vcf@fix[,5] != "."))
  expect_true(all(vcf@fix[,4] != vcf@fix[,5]))
  expect_true(all(vcf@fix[1:5,4] == c("A", "G", "G", "G", "C")))
  expect_true(all(vcf@fix[1:5,5] == c("T", "A", "A", "A", "A")))

  rm(vcf)
  rm(temp)
})
