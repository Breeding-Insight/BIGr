context("Dosage Report to VCF")

read_dosage2vcf_body <- function(output_file) {
  vcf_lines <- readLines(paste0(output_file, ".vcf"))
  body_lines <- vcf_lines[!grepl("^##", vcf_lines)]
  read.table(text = paste(body_lines, collapse = "\n"),
             header = TRUE,
             sep = "\t",
             quote = "",
             comment.char = "",
             check.names = FALSE,
             stringsAsFactors = FALSE,
             colClasses = "character")
}

write_dart_counts_fixture <- function(path, sample_order = c("S2", "S1", "S3")) {
  meta_cols <- c("MarkerName", "AlleleSequence", "Variant", "CallRate",
                 "OneRatioRef", "OneRatioSnp", "FreqHomRef", "FreqHomSnp",
                 "FreqHets", "PICRef", "PICSnp", "AvgPIC", "AvgCountRef",
                 "AvgCountSnp", "RatioAvgCountRefAvgCountSnp")
  skipped <- rep(paste(rep("*", length(meta_cols) + length(sample_order)), collapse = ","), 6)
  sample_header <- paste(c(rep("*", length(meta_cols)), sample_order), collapse = ",")
  target_header <- paste(c(meta_cols, paste0("T", seq_along(sample_order))), collapse = ",")
  fixed_values <- rep("1", length(meta_cols) - 3)
  rows <- c(
    paste(c("Chr01_000120735", "AAA", "", fixed_values, c("10", "20", "30")), collapse = ","),
    paste(c("Chr01_000120735", "ATA", "-:A>T", fixed_values, c("1", "2", "3")), collapse = ","),
    paste(c("Chr02_000000005", "GGG", "", fixed_values, c("40", "50", "60")), collapse = ","),
    paste(c("Chr02_000000005", "GCG", "-:G>C", fixed_values, c("4", "5", "6")), collapse = ",")
  )
  writeLines(c(skipped, sample_header, target_header, rows), path)
}

write_snp_1row_fixture <- function(path, sample_order = c("S1", "S2", "S3")) {
  meta_cols <- c("MarkerName", "AlleleSequenceRef", "AlleleSequenceAlt", "Variant",
                 "CallRate", "OneRatioRef", "OneRatioSnp", "FreqHomRef",
                 "FreqHomSnp", "FreqHets", "PICRef", "PICSnp", "AvgPIC",
                 "AvgCountRef", "AvgCountSnp", "RatioAvgCountRefAvgCountSnp")
  skipped <- rep(paste(rep("*", length(meta_cols) + length(sample_order)), collapse = ","), 6)
  header <- paste(c(meta_cols, sample_order), collapse = ",")
  fixed_values <- rep("1", length(meta_cols) - 4)
  rows <- c(
    paste(c("Chr01_000120735", "AAA", "ATA", "-:A>T", fixed_values, c("0", "1", "2")), collapse = ","),
    paste(c("Chr02_000000005", "GGG", "GCG", "-:G>C", fixed_values, c("-", "", "0")), collapse = ",")
  )
  writeLines(c(skipped, header, rows), path)
}

write_snp_2row_fixture <- function(path, sample_order = c("S1", "S2", "S3")) {
  meta_cols <- c("MarkerName", "AlleleSequence", "Variant", "CallRate",
                 "OneRatioRef", "OneRatioSnp", "FreqHomRef", "FreqHomSnp",
                 "FreqHets", "PICRef", "PICSnp", "AvgPIC", "AvgCountRef",
                 "AvgCountSnp", "RatioAvgCountRefAvgCountSnp")
  skipped <- rep(paste(rep("*", length(meta_cols) + length(sample_order)), collapse = ","), 6)
  header <- paste(c(meta_cols, sample_order), collapse = ",")
  fixed_values <- rep("1", length(meta_cols) - 3)
  rows <- c(
    paste(c("Chr01_000120735", "AAA", "", fixed_values, c("1", "0", "1")), collapse = ","),
    paste(c("Chr01_000120735", "ATA", "-:A>T", fixed_values, c("0", "1", "1")), collapse = ","),
    paste(c("Chr02_000000005", "GGG", "", fixed_values, c("-", "-", "1")), collapse = ","),
    paste(c("Chr02_000000005", "GCG", "-:G>C", fixed_values, c("-", "-", "0")), collapse = ",")
  )
  writeLines(c(skipped, header, rows), path)
}


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

test_that("SNP/INDEL 1 row reports use diploid genotype codes and align counts by sample name", {
  report <- tempfile(fileext = ".csv")
  counts <- tempfile(fileext = ".csv")
  output_file <- tempfile()
  write_snp_1row_fixture(report)
  write_dart_counts_fixture(counts)

  dosage2vcf(dart.report = report, dart.counts = counts, ploidy = 2, output.file = output_file)
  vcf <- read_dosage2vcf_body(output_file)

  expect_equal(colnames(vcf)[10:12], c("S1", "S2", "S3"))
  expect_equal(vcf$`#CHROM`, c("Chr01", "Chr02"))
  expect_equal(vcf$POS, c("000120735", "000000005"))
  expect_equal(vcf$ID, c("Chr01_000120735", "Chr02_000000005"))
  expect_equal(vcf$S1[1], "0/0:2:22:20")
  expect_equal(vcf$S2[1], "1/1:0:11:10")
  expect_equal(vcf$S3[1], "0/1:1:33:30")
  expect_equal(vcf$S1[2], "./.:.:55:50")
  expect_equal(vcf$S2[2], "./.:.:44:40")
  expect_equal(vcf$S3[2], "0/0:2:66:60")
})

test_that("SNP/INDEL 2 row reports collapse to the same genotype codes as 1 row reports", {
  snp_1row <- tempfile(fileext = ".csv")
  snp_2row <- tempfile(fileext = ".csv")
  counts <- tempfile(fileext = ".csv")
  output_1row <- tempfile()
  output_2row <- tempfile()
  write_snp_1row_fixture(snp_1row)
  write_snp_2row_fixture(snp_2row)
  write_dart_counts_fixture(counts)

  dosage2vcf(dart.report = snp_1row, dart.counts = counts, ploidy = 2, output.file = output_1row)
  dosage2vcf(dart.report = snp_2row, dart.counts = counts, ploidy = 2, output.file = output_2row)

  expect_equal(read_dosage2vcf_body(output_2row), read_dosage2vcf_body(output_1row))
})

test_that("SNP/INDEL reports require diploid ploidy", {
  report <- tempfile(fileext = ".csv")
  counts <- tempfile(fileext = ".csv")
  write_snp_1row_fixture(report)
  write_dart_counts_fixture(counts)

  expect_error(
    dosage2vcf(dart.report = report, dart.counts = counts, ploidy = 4, output.file = tempfile()),
    "diploid genotype reports"
  )
})
