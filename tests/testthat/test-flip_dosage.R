context("Flip Dosage Values")


test_that("test dosage flip",{
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

  rm(vcf)

  #Flip the dosage values
  geno_mat_flipped <- flip_dosage(df= geno_mat, ploidy = ploidy, is.reference = TRUE)
  geno_mat_flipped_non <- flip_dosage(df= geno_mat, ploidy = ploidy, is.reference = FALSE)

  #Checks
  expect_true(is(geno_mat_flipped, "matrix"))
  expect_true(all(as.numeric(table(geno_mat)) == rev(as.numeric(table(geno_mat_flipped)))))
  expect_true(all(dim(geno_mat) == dim(geno_mat_flipped)))
  expect_true(all(geno_mat == geno_mat_flipped_non)) #Confirm if the dosages aren't flipped that the output will equal the input

})
