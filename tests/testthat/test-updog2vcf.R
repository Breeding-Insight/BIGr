context("Updog to VCF")


test_that("test updog conversion",{
  #Input variables
  load(system.file("extdata", "iris-multidog.rdata", package = "BIGr"))

  temp_file <- tempfile()

  # Convert updog to VCF
  updog2vcf(
    multidog.object = mout,
    output.file = temp_file,
    updog_version = "0.0.0",
    compress = TRUE,
    reference = "test"
  )

  vcf_result <- read.vcfR(paste0(temp_file,".vcf.gz"), verbose = FALSE)

  DP <- sum(as.numeric(extract.gt(vcf_result, "DP")))

  expect_equal(DP, 11645)

  MPP <- sum(as.numeric(extract.gt(vcf_result, "MPP")))

  expect_equal(MPP, 39.99, tolerance = 0.01)

})

test_that("reference sets the ##reference line of the header", {
  load(system.file("extdata", "iris-multidog.rdata", package = "BIGr"))

  temp_file <- tempfile()
  updog2vcf(
    multidog.object = mout,
    output.file = temp_file,
    updog_version = "0.0.0",
    compress = FALSE,
    reference = "GCF_test.1"
  )

  meta <- grep("^##", readLines(paste0(temp_file, ".vcf")), value = TRUE)
  expect_equal(sum(meta == "##reference=GCF_test.1"), 1)
  expect_false(any(grepl("^##contig=", meta)))
  expect_true(any(grepl("^##BIGrCommandLine.updog2vcf", meta)))
})

test_that("SNPs without a position in the SNP ID are kept with the ID as CHROM and POS 1", {
  load(system.file("extdata", "iris-multidog.rdata", package = "BIGr"))

  # One ID without an underscore and two with a position that is not a number
  new_ids <- c(Chr1_10 = "markerA", Chr1_100 = "Chr1_abc", Chr4_810 = "Chr1_xyz")
  mout$snpdf$snp <- ifelse(mout$snpdf$snp %in% names(new_ids), new_ids[mout$snpdf$snp], mout$snpdf$snp)
  mout$inddf$snp <- ifelse(mout$inddf$snp %in% names(new_ids), new_ids[mout$inddf$snp], mout$inddf$snp)

  temp_file <- tempfile()
  expect_warning(
    updog2vcf(multidog.object = mout, output.file = temp_file, updog_version = "0.0.0",
              compress = FALSE, reference = "GCF_test.1"),
    "3 SNPs have a position that could not be read from the SNP ID"
  )

  lines <- readLines(paste0(temp_file, ".vcf"))
  body <- do.call(rbind, strsplit(lines[!startsWith(lines, "#")], "\t"))
  expect_equal(nrow(body), 10)
  changed <- body[body[, 3] %in% new_ids, , drop = FALSE]
  expect_equal(changed[, 1], changed[, 3])
  expect_true(all(changed[, 2] == "1"))
  expect_false(any(body[, 2] %in% c("0", "NA")))
})
