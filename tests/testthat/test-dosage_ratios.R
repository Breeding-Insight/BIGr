context("Dosage Ratios")


test_that("test ratios",{
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

  #Calculate percentages
  ratios_df <- dosage_ratios(geno_mat, ploidy = ploidy)

  #Get means for each dosage value
  dosage_means <- ratios_df %>%
    group_by(Dosage) %>%
    summarize(Mean_Percentage = mean(Percentage))

  #Checks
  expect_s3_class(ratios_df, "data.frame")
  expect_true(all(ratios_df$Percentage >= 0 & ratios_df$Percentage <= 100))
  expect_true(all(round(dosage_means$Mean_Percentage,5) == c(34.30327, 58.71075, 6.98597)))

})
