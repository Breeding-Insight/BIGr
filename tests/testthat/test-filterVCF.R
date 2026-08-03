context("Filtering")


test_that("Filtering with Updog metrics",{

  #Variables
  filter_ploidy <- 2
  filter_maf <- 0.05
  size_depth <- 10
  snp_miss <- 50
  sample_miss <- 50
  OD_filter <- 0.05
  Bias <- c(0.5, 2)
  Bias_min <- Bias[1]
  Bias_max <- Bias[2]
  Prop_mis <- 0.05
  maxpostprob_filter <- 0.5
  max_post <- maxpostprob_filter
  output_name <- "out"
  snp_miss <- snp_miss/100
  sample_miss <- sample_miss/100
  ploidy <- filter_ploidy
  maf_filter <- filter_maf

  input <- filtering_files <- list()
  input$updog_rdata$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGr")

  temp_file <- tempfile(fileext = ".vcf.gz")

  #Input file
  vcf <- read.vcfR(input$updog_rdata$datapath, verbose = FALSE)

  # Identify if have updog parameters
  format_fields <- unique(vcf@gt[,1])
  info_fields <- vcf@fix[1,8]

  updog_par <- grepl("MPP", format_fields) & grepl("PMC", info_fields) & grepl("BIAS", info_fields)

  #Starting SNPs
  starting_snps <- nrow(vcf)
  #export INFO dataframe
  filtering_files$raw_vcf_df <- data.frame(vcf@fix)

  #Filtering
  vcf <- filterVCF(vcf.file = vcf,
                   ploidy=ploidy,
                   output.file=NULL,
                   filter.OD = OD_filter,
                   filter.BIAS.min = Bias_min,
                   filter.BIAS.max = Bias_max,
                   filter.DP = as.numeric(size_depth),
                   filter.PMC = Prop_mis,
                   filter.SAMPLE.miss = as.numeric(sample_miss),
                   filter.SNP.miss = as.numeric(snp_miss),
                   filter.MAF = as.numeric(maf_filter),
                   filter.MPP = max_post)

  #Getting missing data information
  #Add support for genotype matrix filtering?
  gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  filtering_files$snp_miss_df <- rowMeans(is.na(gt_matrix)) #SNP missing values
  filtering_files$sample_miss_df <- as.numeric(colMeans(is.na(gt_matrix))) #Sample missing values

  #These counts increased when the INFO parser was fixed to read scientific
  #notation (see the final_snps note below) - more SNPs are retained, so each
  #sample has more genotype calls.
  expect_true(all(table(gt_matrix[,10]) == c(37,36,18)))

  rm(gt_matrix) #Remove gt matrix

  #Writing file
  write.vcf(vcf, file = temp_file)

  #Get final_snps
  #Was 43 while filterVCF()'s INFO parser could not read scientific notation:
  #110 of the 499 SNPs in the test file have a PMC written as e.g. 4.78e-07,
  #which parsed as 4.78 and so failed filter.PMC = 0.05 despite being well
  #under it. Those SNPs are now retained, hence 93 rather than 43.
  final_snps <- nrow(vcf)
  expect_equal(final_snps, 93)

})

test_that("INFO filters read scientific notation",{

  #updog2vcf() writes INFO values with paste0(), so R's default formatting emits
  #scientific notation for small values. filterVCF() must read those as written.
  vcf <- read.vcfR(system.file("iris_DArT_VCF.vcf.gz", package = "BIGr"), verbose = FALSE)
  vcf <- vcf[1:4, ]

  #Two SNPs written in fixed notation, two in scientific notation. All four have
  #an OD and a PMC that pass the filters below.
  vcf@fix[,"INFO"] <- c("DP=100;ADS=60,40;BIAS=1.0;OD=0.001;PMC=0.001",
                        "DP=100;ADS=60,40;BIAS=1.0;OD=1.2e-05;PMC=4.78e-07",
                        "DP=100;ADS=60,40;BIAS=1.0;OD=0.002;PMC=0.002",
                        "DP=100;ADS=60,40;BIAS=1.0;OD=3.4e-08;PMC=9.10e-09")

  expect_equal(nrow(filterVCF(vcf, filter.PMC = 0.05, ploidy = 2)), 4)
  expect_equal(nrow(filterVCF(vcf, filter.OD = 0.05, ploidy = 2)), 4)

  #A value that genuinely exceeds the threshold is still removed, so the parser
  #is reading the exponent rather than ignoring it.
  vcf@fix[2,"INFO"] <- "DP=100;ADS=60,40;BIAS=1.0;OD=1.2e+01;PMC=5.0e-01"
  expect_equal(nrow(filterVCF(vcf, filter.PMC = 0.05, ploidy = 2)), 3)
  expect_equal(nrow(filterVCF(vcf, filter.OD = 0.05, ploidy = 2)), 3)

  #The field name is matched against a whole INFO entry, not as a substring, so
  #a longer field ending in the same name is not picked up by mistake.
  vcf@fix[,"INFO"] <- paste0("XOD=99;", vcf@fix[,"INFO"])
  expect_equal(nrow(filterVCF(vcf, filter.OD = 0.05, ploidy = 2)), 3)

})
