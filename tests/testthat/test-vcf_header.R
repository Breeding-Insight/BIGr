context("VCF header reference genome lines")

write_header_fixture <- function(path, contigs = c("Chr01", "Chr02"), reference = TRUE) {
  writeLines(c("##fileformat=VCFv4.3",
               if (reference) "##reference=file:test_genome.fna",
               paste0("##contig=<ID=", contigs, ",length=1000000>"),
               '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
               '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">'),
             path)
}

test_that("vcf_header gives the reference and contig lines only", {
  header <- tempfile(fileext = ".txt")
  write_header_fixture(header)

  res <- vcf_reference_lines(vcf_header = header)
  expect_equal(res$lines, c("##reference=file:test_genome.fna",
                            "##contig=<ID=Chr01,length=1000000>",
                            "##contig=<ID=Chr02,length=1000000>"))
  expect_equal(res$contigs, c("Chr01", "Chr02"))

  # reference replaces the ##reference line of the file
  res <- vcf_reference_lines(vcf_header = header, reference = "GCF_020283565.1")
  expect_equal(res$lines[1], "##reference=GCF_020283565.1")
  expect_equal(length(res$lines), 3)
})

test_that("reference alone gives the ##reference line", {
  res <- vcf_reference_lines(reference = "GCF_020283565.1")
  expect_equal(res$lines, "##reference=GCF_020283565.1")
  expect_null(res$contigs)
})

test_that("no reference information gives no lines and a warning", {
  expect_warning(res <- vcf_reference_lines(), "No reference genome information")
  expect_equal(res$lines, character(0))
  expect_null(res$contigs)
  expect_no_warning(vcf_reference_lines(warn = FALSE))
})

test_that("invalid vcf_header and reference inputs give errors", {
  expect_error(vcf_reference_lines(vcf_header = tempfile()), "vcf_header file not found")

  no_ref <- tempfile(fileext = ".txt")
  writeLines(c("##fileformat=VCFv4.3", '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">'), no_ref)
  expect_error(vcf_reference_lines(vcf_header = no_ref), "no ##reference or ##contig lines")

  two_ref <- tempfile(fileext = ".txt")
  writeLines(c("##reference=a.fna", "##reference=b.fna"), two_ref)
  expect_error(vcf_reference_lines(vcf_header = two_ref), "more than one ##reference line")

  bad_contig <- tempfile(fileext = ".txt")
  writeLines(c("##reference=a.fna", "##contig=<length=100>"), bad_contig)
  expect_error(vcf_reference_lines(vcf_header = bad_contig), "without an ID")

  expect_error(vcf_reference_lines(reference = c("a", "b")), "single character string")
  expect_error(vcf_reference_lines(reference = ""), "single character string")
})

test_that("chromosomes missing from the contig lines give a warning", {
  expect_warning(check_vcf_contigs(c("Chr01", "Chr03", "Chr03"), c("Chr01", "Chr02")),
                 "1 chromosome\\(s\\) in the VCF are not listed.*Chr03")
  expect_no_warning(check_vcf_contigs(c("Chr01", "Chr02"), c("Chr01", "Chr02")))
  expect_no_warning(check_vcf_contigs(c("Chr01", "Chr03"), NULL))
})

test_that("the meta lines of a VCF file can be replaced", {
  vcf <- tempfile(fileext = ".vcf")
  body <- c("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "Chr01\t100\tChr01_000000100\tA\tG\t.\t.\tDP=10")
  writeLines(c("##fileformat=VCFv4.3", "##source=polyRADv2", "##contig=<ID=Chr01>", body), vcf)

  expect_equal(read_vcf_meta(vcf), c("##fileformat=VCFv4.3", "##source=polyRADv2", "##contig=<ID=Chr01>"))

  new_meta <- c("##fileformat=VCFv4.3", "##reference=file:test_genome.fna",
                "##contig=<ID=Chr01,length=1000000>", "##source=polyRADv2")
  replace_vcf_meta(vcf, new_meta)
  expect_equal(readLines(vcf), c(new_meta, body))
})

test_that("reference genome lines of an existing header are replaced", {
  meta <- c("##fileformat=VCFv4.3", "##source=polyRADv2", "##contig=<ID=Chr01>", "##contig=<ID=Chr02>")
  ref_lines <- list(lines = c("##reference=file:test_genome.fna", "##contig=<ID=Chr01,length=1000000>"),
                    contigs = "Chr01")
  expect_equal(update_vcf_reference_lines(meta, ref_lines),
               c("##fileformat=VCFv4.3", "##reference=file:test_genome.fna",
                 "##contig=<ID=Chr01,length=1000000>", "##source=polyRADv2"))

  # A reference ID alone keeps the existing contig lines
  ref_lines <- list(lines = "##reference=GCF_test.1", contigs = NULL)
  expect_equal(update_vcf_reference_lines(meta, ref_lines),
               c("##fileformat=VCFv4.3", "##reference=GCF_test.1", "##source=polyRADv2",
                 "##contig=<ID=Chr01>", "##contig=<ID=Chr02>"))
})

test_that("quotes and backslashes are escaped for quoted header values", {
  expect_equal(vcf_quote_escape('> f("a.csv")'), '> f(\\"a.csv\\")')
  expect_equal(vcf_quote_escape("C:\\data\\a.csv"), "C:\\\\data\\\\a.csv")
  expect_equal(vcf_quote_escape("no quotes"), "no quotes")
})
