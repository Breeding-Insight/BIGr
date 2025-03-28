context("Get OffTargets")


test_that("test madc offtargets",{
  #Input variables
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  bot_file <- system.file("example_SNPs_DArTag-probe-design_f180bp.botloci", package="BIGr")
  db_file <- system.file("example_allele_db.fa", package="BIGr")

  #Calculations
  temp <- tempfile(fileext = ".vcf")
  temp_multi <- tempfile(fileext = ".vcf")

  set.seed(123)
  get_OffTargets(madc = madc_file,
                 botloci = bot_file,
                 hap_seq = db_file,
                 n.cores = 2,
                 rm_multiallelic_SNP = FALSE,
                 multiallelic_SNP_dp_thr = 0,
                 multiallelic_SNP_sample_thr = 0,
                 out_vcf = temp,
                 verbose = FALSE)

  set.seed(456)
  get_OffTargets(madc = madc_file,
                 botloci = bot_file,
                 hap_seq = db_file,
                 n.cores = 2,
                 rm_multiallelic_SNP = TRUE,
                 multiallelic_SNP_dp_thr = 0,
                 multiallelic_SNP_sample_thr = 0,
                 out_vcf = temp_multi,
                 verbose = FALSE)

  vcf <- read.vcfR(temp)
  vcf_multi <- read.vcfR(temp_multi)

  #Check
  expect_true(all(dim(vcf@gt) == c("33","11")))
  expect_true(all(dim(vcf_multi@gt) == c("32","11")))

  rm(vcf)
  rm(vcf_multi)

})
