context("Calculating Heterozygosity")


test_that("test heterozygosity",{
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
  het_df <- calculate_Het(data.frame(geno_mat, check.names=FALSE), ploidy = ploidy)

  #Check
  expect_s3_class(het_df, "data.frame")
  expect_true(all(het_df$ObservedHeterozygosity >= 0 & het_df$ObservedHeterozygosity <= 1))
  expect_true(round(mean(het_df$ObservedHeterozygosity), 7) == 0.5871075)
  rm(vcf) #Remove VCF

})
