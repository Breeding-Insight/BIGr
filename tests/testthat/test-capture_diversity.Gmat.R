context("Capturing Diversity Gmat")


test_that("test capture diversity",{
  #Input variables
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

  #Estimate diversity to capture
  set.seed(123) #Set seed for consistency
  result_df <- capture_diversity.Gmat(df = geno_mat, ploidy = ploidy, r2_threshold=0.9, iterations = 10, sample_list = NULL, parallel=FALSE, save.result=FALSE)

  #Check
  expect_s3_class(result_df, "data.frame")
  expect_true(result_df$Individuals == 8.4)
  rm(vcf) #Remove VCF

})
