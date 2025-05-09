context("Calculating MAF")


test_that("test diversity",{
  #Input variables (need to add support for VCF file)
  input <- diversity_items <- list()
  input$diversity_ploidy <- 2
  input$diversity_file$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGr")

  ploidy <- as.numeric(input$diversity_ploidy)
  geno <- input$diversity_file$datapath

  #Import genotype information if in VCF format
  vcf <- read.vcfR(geno, verbose = FALSE)

  #Extract GT and convert to numeric calls
  geno_mat <- extract.gt(vcf, element = "GT")
  geno_mat <- apply(geno_mat, 2, convert_to_dosage)

  #Calculate maf
  maf_df <- calculate_MAF(geno_mat, ploidy = ploidy)

  #Check
  expect_s3_class(maf_df, "data.frame")
  expect_true(all(maf_df$MAF >= 0 & maf_df$MAF <= 0.5))
  expect_true(round(mean(maf_df$MAF), 7) == 0.3134469)
  expect_true(all(round(data.frame(maf(vcf))$Frequency,5) == round(maf_df$MAF,5)))
  rm(vcf) #Remove VCF

})
