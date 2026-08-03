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
  #a longer field ending in the same name is not picked up by mistake. The
  #look-alike is placed after the real OD on purpose: the pattern this replaced
  #searched greedily and so took the rightmost match, meaning a look-alike placed
  #before the real field would be read correctly by accident.
  vcf@fix[,"INFO"] <- sub("(OD=[^;]*)", "\\1;MOD=77", vcf@fix[,"INFO"])
  expect_equal(nrow(filterVCF(vcf, filter.OD = 0.05, ploidy = 2)), 3)

})

test_that("Variants with unreadable filter values are removed, not inserted as NA",{

  #An NA in a logical index does not drop a row, it inserts an all-NA row, so a
  #variant whose filter value cannot be read must be removed explicitly.
  vcf <- read.vcfR(system.file("iris_DArT_VCF.vcf.gz", package = "BIGr"), verbose = FALSE)
  vcf <- vcf[1:4, ]
  starting_snps <- nrow(vcf)

  #Record 3 carries no PMC at all.
  vcf@fix[,"INFO"] <- "DP=100;ADS=60,40;BIAS=1.0;OD=0.001;PMC=0.001"
  vcf@fix[3,"INFO"] <- "DP=100;ADS=60,40;BIAS=1.0;OD=0.001"

  expect_warning(filtered <- filterVCF(vcf, filter.PMC = 0.05, ploidy = 2),
                 "1 of 4 variants had a missing or unreadable PMC value")
  expect_equal(nrow(filtered), starting_snps - 1)
  expect_false(any(is.na(filtered@fix[,"CHROM"])))
  expect_false(any(is.na(filtered@fix[,"POS"])))

  #The same applies to MAF, which vcfR::maf() reports as NA for a variant left
  #with no called genotypes - something the DP and MPP masking can produce. Only
  #the variants without calls should go, not the whole file.
  vcf@fix[,"INFO"] <- "DP=100;ADS=60,40;BIAS=1.0;OD=0.001;PMC=0.001"
  vcf@gt[2:3, -1] <- "./.:.:.:.:.:."
  expect_warning(masked <- filterVCF(vcf, filter.MAF = 0.05, ploidy = 2),
                 "2 of 4 variants had a missing or unreadable MAF value")
  expect_lte(nrow(masked), 2)
  expect_false(any(is.na(masked@fix[,"CHROM"])))

  #A file with nothing unreadable must not warn (progress messages are expected).
  expect_no_warning(suppressMessages(filterVCF(vcf, filter.PMC = 0.05, ploidy = 2)))

})

test_that("A requested INFO filter runs whatever the first record contains",{

  #Which INFO fields exist used to be read from the first record alone, so a
  #field missing from that one record disabled its filter for the whole file,
  #silently. A field may legally be absent from any individual record.
  vcf <- read.vcfR(system.file("iris_DArT_VCF.vcf.gz", package = "BIGr"), verbose = FALSE)
  vcf <- vcf[1:4, ]
  vcf@fix[,"INFO"] <- "DP=100;ADS=60,40;BIAS=1.0;OD=0.99;PMC=0.001"
  vcf@fix[1,"INFO"] <- "DP=100;ADS=60,40;BIAS=1.0;PMC=0.001"

  #Every record either fails OD < 0.05 or has no OD at all, so none survive.
  expect_warning(filtered <- filterVCF(vcf, filter.OD = 0.05, ploidy = 2),
                 "1 of 4 variants had a missing or unreadable OD value")
  expect_equal(nrow(filtered), 0)

  #Reordering the records must not change the result.
  reordered <- vcf
  reordered@fix[,"INFO"] <- reordered@fix[c(2,3,4,1),"INFO"]
  expect_warning(also <- filterVCF(reordered, filter.OD = 0.05, ploidy = 2),
                 "1 of 4 variants had a missing or unreadable OD value")
  expect_equal(nrow(also), nrow(filtered))

  #A filter that was never requested stays inactive.
  expect_no_warning(suppressMessages(filterVCF(vcf, filter.MAF = 0.05, ploidy = 2)))

})

test_that("A filter whose values are all unreadable removes every variant",{

  #Every filter behaves the same way here. Previously OD, BIAS and PMC returned
  #the data unfiltered in this situation, which reported success while applying
  #no filter at all.
  vcf <- read.vcfR(system.file("iris_DArT_VCF.vcf.gz", package = "BIGr"), verbose = FALSE)
  vcf <- vcf[1:4, ]

  vcf@fix[,"INFO"] <- "DP=100;ADS=60,40;BIAS=1.0;OD=.;PMC=0.001"
  expect_warning(dropped <- filterVCF(vcf, filter.OD = 0.05, ploidy = 2),
                 "No readable OD values were found")
  expect_equal(nrow(dropped), 0)

  #The warning names the INFO column, so the cause is identifiable.
  expect_warning(filterVCF(vcf, filter.OD = 0.05, ploidy = 2),
                 "INFO column")

  #MAF reports the cause that actually applies to it, rather than an INFO field.
  vcf@fix[,"INFO"] <- "DP=100;ADS=60,40;BIAS=1.0;OD=0.001;PMC=0.001"
  expect_warning(masked <- filterVCF(vcf, filter.DP = 1e6, filter.MAF = 0.05, ploidy = 2),
                 "no variant has any called genotypes")
  expect_equal(nrow(masked), 0)

})
