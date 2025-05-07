context("MADC to VCF")


test_that("test madc conversion",{
  #Input variables
  input <- diversity_items <- list()
  input$diversity_ploidy <- 2
  input$madc_file$datapath <- system.file("iris_DArT_MADC.csv", package = "BIGr")

  ploidy <- as.numeric(input$diversity_ploidy)
  madc <- input$madc_file$datapath

  #Temp location
  output_file1 <- tempfile(fileext = ".vcf")
  #output_file2 <- tempfile(fileext = ".vcf")

  #Convert the dart files to vcf
  madc2vcf(madc_file = madc, output.file = output_file1, get_REF_ALT = FALSE)
  #madc2vcf(madc_file = madc, output.file = output_file2, get_REF_ALT = TRUE) #Update the example MADC to support this eval

  #Test validity of VCF
  vcf <- read.vcfR(output_file1, verbose = FALSE)

  DP_df <- extract.gt(vcf, element = "DP")
  DP_df_mean <- mean(as.numeric(DP_df))

  #Checks
  expect_s4_class(vcf, "vcfR")
  expect_true(nrow(vcf) == 500)
  expect_true(ncol(vcf@gt) == 151)
  expect_true(all(dim(DP_df) == c(500,150)))
  expect_true(round(DP_df_mean, 3) == 314.92)

  rm(vcf)
  rm(output_file1)

})
