context("Dosage Report to VCF")


test_that("test dosage report conversion",{
  #Input variables
  input <- diversity_items <- list()
  input$diversity_ploidy <- 2
  input$report_file$datapath <- system.file("iris_DArT_Allele_Dose_Report.csv", package = "BIGr")
  input$counts_file$datapath <- system.file("iris_DArT_Counts.csv", package = "BIGr")

  ploidy <- as.numeric(input$diversity_ploidy)
  report <- input$report_file$datapath
  counts <- input$counts_file$datapath

  #Temp location
  output_file <- tempfile()

  #Convert the dart files to vcf
  dosage2vcf(dart.report = report, dart.counts = counts, ploidy = ploidy, output.file = output_file)

  #Test validity of VCF
  vcf <- read.vcfR(paste0(output_file,".vcf"), verbose = FALSE)

  #Checks
  expect_s4_class(vcf, "vcfR")
  expect_true(nrow(vcf) == 500)
  expect_true(ncol(vcf@gt) == 151)
  expect_true(round(mean(data.frame(maf(vcf))$Frequency), 5) == 0.44771)

  rm(vcf)
  rm(output_file)

})
