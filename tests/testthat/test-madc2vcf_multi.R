context("MADC to VCF via polyRAD")

# =======================================================================
# Using Breeding-Insight/BIGapp-PanelHub test files
# =======================================================================

test_that("madc2vcf_multi — alfalfa (BIGapp-PanelHub)", {

  skip_if_not_installed("polyRAD")
  skip_if_not_installed("VariantAnnotation")

 github_path <- "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/long_seq/"

  # External alfalfa test files
  alfalfa_madc           <- paste0(github_path, "test_madcs/alfalfa_madc.csv")
  alfalfa_madc_wrongID   <- paste0(github_path, "test_madcs/alfalfa_madc_wrongID.csv")
  alfalfa_madc_raw       <- paste0(github_path, "test_madcs/alfalfa_madc_raw.csv")       # raw DArT format (7-row header)
  alfalfa_iupac          <- paste0(github_path, "test_madcs/alfalfa_IUPAC.csv")
  alfalfa_lowercase      <- paste0(github_path, "test_madcs/alfalfa_lowercase.csv")
  alfalfa_botloci        <- paste0(github_path, "alfalfa/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_f180bp.botloci")          # botloci for alfalfa
  alfalfa_markers_info   <- paste0(github_path, "alfalfa/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_snpID_lut.csv") # markers_info: CloneID/BI_markerID, Chr, Pos, Ref, Alt
  alfalfa_markers_info_ChromPos   <- paste0(github_path, "test_madcs/alfalfa_marker_info_ChromPos.csv") # markers_info: CloneID/BI_markerID, Chr, Pos
  alfalfa_microhapDB <- paste0(github_path, "alfalfa/alfalfa_allele_db_v001.fa")

  # External potato test files
  potato_indel_madc                 <- paste0(github_path, "test_madcs/potato_indel_madc.csv")
  potato_indel_iupac                <- paste0(github_path, "test_madcs/potato_indel_IUPAC.csv")
  potato_indel_lowercase            <- paste0(github_path, "test_madcs/potato_indel_lowercase.csv")
  potato_more_indels_chrompos_false <- paste0(github_path, "test_madcs/potato_more_indels_madc_ChromPosFALSE.csv")
  potato_botloci                    <- paste0(github_path, "potato/potato_dartag_v2_3915markers_rm7dupTags_6traitMarkers_f150bp_ref_alt.botloci")
  potato_markers_info               <- paste0(github_path, "potato/potato_dartag_v2_3915markers_rm7dupTags_6traitMarkers_rm1dup_snpID_lut.csv") # CloneID/BI_markerID, Chr, Pos, Ref, Alt
  potato_markers_info_ChromPos   <- paste0(github_path, "test_madcs/potato_marker_info_chrompos.csv") # markers_info: CloneID/BI_markerID, Chr, Pos
  potato_microhapDB <- paste0(github_path, "potato/potato_allele_db_v001.fa")

  skip_if_offline("raw.githubusercontent.com")

  out <- tempfile(fileext = ".vcf")
  # Fixed allele ID format
  expect_no_error(
    madc2vcf_multi(
      madc_file    = alfalfa_madc,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      ploidy       = 4L,
      verbose      = TRUE,
      reference    = "test"
    )
  )

  vcf <- read.vcfR(out, verbose = FALSE)
  expect_s4_class(vcf, "vcfR")
  expect_equal(sum(grepl(",", vcf@fix[,5])), 281)
  GT <- extract.gt(vcf)
  expect_equal(GT[3,5],"0/0/0/3")

  # Bottom strand markers are reverse complemented by polyRAD (botloci must be
  # passed to polyRAD as IDs, not as a file path) and DArTstrand follows botloci
  bot_ids <- readLines(alfalfa_botloci)
  strand <- sub(".*DArTstrand=(.).*", "\\1", vcf@fix[, "INFO"])
  expect_equal(unname(strand == "-"), unname(vcf@fix[, "ID"] %in% bot_ids))
  madc <- read.csv(alfalfa_madc)
  bot_rec <- which(strand == "-")[1]
  madc_ref <- madc$AlleleSequence[madc$AlleleID == paste0(vcf@fix[bot_rec, "ID"], "|Ref_0001")]
  expect_equal(unname(vcf@fix[bot_rec, "REF"]),
               as.character(Biostrings::reverseComplement(Biostrings::DNAString(toupper(madc_ref)))))

  # Don't allow raw MADC
  out <- tempfile(fileext = ".vcf")
  expect_error(
    madc2vcf_multi(
      madc_file    = alfalfa_madc_raw,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      ploidy       = 4L,
      verbose      = FALSE,
      reference    = "test"
    ), regexp = "The MADC file does not have fixed AlleleIDs. Please process the MADC file through HapApp before using this function."
  )

  out <- tempfile(fileext = ".vcf")
  expect_no_error(
    madc2vcf_multi(
      madc_file    = alfalfa_madc,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      ploidy       = 4L,
      verbose      = TRUE,
      reference    = "test"
    )
  )

  # Wrong IDs
  out <- tempfile(fileext = ".vcf")
  expect_error(
    madc2vcf_multi(
      madc_file    = alfalfa_madc_wrongID,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      ploidy       = 4L,
      verbose      = TRUE,
      reference    = "test"
    ), regexp = "Check marker IDs in both MADC and botloci files. They should be the same."
  )

  out <- tempfile(fileext = ".vcf")
  expect_no_error(
    madc2vcf_multi(
      madc_file    = alfalfa_madc_wrongID,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      markers_info = alfalfa_markers_info_ChromPos,
      ploidy       = 4L,
      verbose      = TRUE,
      reference    = "test"
    )
  )

  vcf <- read.vcfR(out, verbose = FALSE)
  expect_s4_class(vcf, "vcfR")
  expect_equal(sum(grepl(",", vcf@fix[,5])), 281)
  GT <- extract.gt(vcf)
  expect_equal(GT[3,5],"0/0/0/3")

  ### Avoid IUPAC codes
  out <- tempfile(fileext = ".vcf")
  expect_error(
    madc2vcf_multi(
      madc_file    = alfalfa_iupac,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      ploidy       = 4L,
      verbose      = TRUE,
      reference    = "test"
    ), regexp = "MADC Allele Sequences contain different IUPAC \\(non-ATCG\\) codes in REF and ALT sequences. Please run HapApp to clean MADC file before using this function."
  )

  out <- tempfile(fileext = ".vcf")
  expect_error(
    madc2vcf_multi(
      madc_file    = alfalfa_lowercase,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      ploidy       = 4L,
      verbose      = TRUE,
      reference    = "test"
    ), regexp = "Not all Ref sequences have a corresponding Alt or vice versa. Please provide a complete MADC file before using this function."
  )

  out <- tempfile(fileext = ".vcf")
  expect_no_error(
    madc2vcf_multi(
      madc_file    = potato_indel_madc,
      botloci_file = potato_botloci,
      outfile      = out,
      markers_info = potato_markers_info_ChromPos,
      ploidy       = 4L,
      verbose      = TRUE,
      reference    = "test"
    )
  )

  vcf <- read.vcfR(out, verbose = FALSE)
  expect_s4_class(vcf, "vcfR")
  expect_equal(sum(grepl(",", vcf@fix[,5])), 277)
  GT <- extract.gt(vcf)
  expect_equal(GT[3,5],"0/1/1/6")

  # Indel alleles are aligned instead of cut by polyRAD: full lengths are kept
  # and every allele is matched to its MADC AlleleID
  indel <- vcf@fix[vcf@fix[, "ID"] == "chr05_004488021", ]
  expect_equal(unname(nchar(c(indel["REF"], strsplit(indel["ALT"], ",")[[1]]))), c(109, 116, 116))
  expect_false(any(grepl("hapID=\\.(,|;|$)|,\\.(,|;|$)", vcf@fix[, "INFO"])))

})


test_that("add_hap_info_records adds MADC IDs and strand", {
  report <- data.frame(
    AlleleID = c("Chr01_000100|Ref_0001", "Chr01_000100|Alt_0002",
                 "Chr02_000200|Ref_0001", "Chr02_000200|Alt_0002"),
    CloneID = c("Chr01_000100", "Chr01_000100", "Chr02_000200", "Chr02_000200"),
    AlleleSequence = c("AAACCCGGGT", "AAACCCTGGT", "AACCGGTTGA", "AACCGATTGA")
  )
  hap_info <- madc_hap_lookup(report, "Chr02_000200")

  # Chr02 is on the bottom strand; its VCF sequences are matched through the
  # reverse complement lookup
  lines <- c(
    "Chr01\t95\t.\tAAACCCGGGT\tAAACCCTGGT\t.\t.\tNS=2\tGT\t0/1",
    "Chr02\t190\t.\tTCAACCGGTT\tTCAATCGGTT\t.\t.\tNS=2\tGT\t0/0",
    "Chr03\t5\t.\tGGGGGGGGGG\t.\t.\t.\tNS=1\tGT\t0/0"
  )

  res <- add_hap_info_records(lines, hap_info)
  f <- do.call(rbind, strsplit(res$lines, "\t"))
  orig <- do.call(rbind, strsplit(lines, "\t"))
  expect_equal(res$n_missing, 1)
  expect_equal(f[, c(1, 2, 4, 5)], orig[, c(1, 2, 4, 5)])  # CHROM, POS, REF, ALT unchanged
  expect_equal(f[, 3], c("Chr01_000100", "Chr02_000200", "."))
  expect_equal(f[1, 8], "NS=2;targetSNP=Chr01_000100;DArTstrand=+;hapID=Chr01_000100|Ref_0001,Chr01_000100|Alt_0002")
  expect_equal(f[2, 8], "NS=2;targetSNP=Chr02_000200;DArTstrand=-;hapID=Chr02_000200|Ref_0001,Chr02_000200|Alt_0002")
  expect_equal(f[3, 8], "NS=1;targetSNP=.;DArTstrand=+;hapID=.")
  expect_equal(f[, 10], c("0/1", "0/0", "0/0"))

  # Two markers sharing a sequence (Chr04 Ref is the reverse complement of the
  # bottom strand Chr02 Ref): the marker must come from LU, not the sequence
  report2 <- rbind(report, data.frame(
    AlleleID = c("Chr04_000300|Ref_0001", "Chr04_000300|Alt_0002"),
    CloneID = c("Chr04_000300", "Chr04_000300"),
    AlleleSequence = c("TCAACCGGTT", "TCAACCGGTA")))
  hap_info <- madc_hap_lookup(report2, "Chr02_000200",
                              loci = c("Chr01_000100", "Chr02_000200", "Chr04_000300"))
  lines <- c("Chr02\t190\t.\tTCAACCGGTT\tTCAATCGGTT\t.\t.\tNS=2;LU=2\tGT\t0/0",
             "Chr04\t290\t.\tTCAACCGGTT\tTCAACCGGTA\t.\t.\tNS=2;LU=3\tGT\t0/0")
  res <- add_hap_info_records(lines, hap_info)
  f <- do.call(rbind, strsplit(res$lines, "\t"))
  expect_equal(res$n_missing, 0)
  expect_equal(f[, 3], c("Chr02_000200", "Chr04_000300"))
  expect_equal(f[1, 8], "NS=2;LU=2;targetSNP=Chr02_000200;DArTstrand=-;hapID=Chr02_000200|Ref_0001,Chr02_000200|Alt_0002")
  expect_equal(f[2, 8], "NS=2;LU=3;targetSNP=Chr04_000300;DArTstrand=+;hapID=Chr04_000300|Ref_0001,Chr04_000300|Alt_0002")

  meta <- add_hap_info_header(c("##fileformat=VCFv4.3", "##INFO=<ID=NS,Number=1>",
                                "##SAMPLE=<ID=S1,Ploidy=2>", "##source=polyRAD"))
  expect_false(any(startsWith(meta, "##SAMPLE")))
  expect_equal(sub("^##INFO=<ID=([^,]+),.*", "\\1", grep("^##INFO", meta, value = TRUE)),
               c("NS", "targetSNP", "DArTstrand", "hapID"))
})

test_that("align_madc_indels aligns indel alleles to Ref with '-' and '.'", {
  clone <- c(rep("M1", 3), rep("M2", 2))
  ids <- c("M1|Ref_0001", "M1|Alt_0002", "M1|AltMatch_0001", "M2|Ref_0001", "M2|Alt_0002")
  seqs <- c("ACGTACGTTTGCAGT",        # Ref
            "ACGTACGTTTCCCGCAGT",     # 3 bp insertion (CCC) after TT
            "ACGTACTTTGCAGT",         # 1 bp deletion (G)
            "AAAACCCC", "AAAACCCG")   # same length: not changed
  res <- align_madc_indels(clone, ids, seqs)
  expect_equal(res$n_markers, 1L)
  # Insertion columns (no Ref base) are "."; deletions are "-"
  expect_equal(res$seqs[1:3], c("ACGTACGTTT...GCAGT",
                                "ACGTACGTTTCCCGCAGT",
                                "ACGTAC-TTT...GCAGT"))
  expect_equal(res$seqs[4:5], seqs[4:5])
  expect_equal(gsub("[-.]", "", res$seqs), seqs)
  expect_equal(length(unique(nchar(res$seqs[1:3]))), 1L)

  # Target SNP index counts only Ref bases, not the gaps
  rad <- list(locTable = data.frame(Ref = c("AC..GTA", "ACGT"), row.names = c("L1", "L2")),
              alleleDepth = matrix(0L, 1, 4, dimnames = list(NULL, c("L1|Ref_0001", "L1|Alt_0002", "L2|Ref_0001", "L2|Alt_0002"))),
              alleles2loc = c(1L, 1L, 2L, 2L),
              alleleNucleotides = c("AC..GTA", "AC..GCA", "ACGT", "TCGT"))
  expect_equal(target_index_in_ref(rad), c(4L, 1L))
})

test_that("target_index_from_madc places target SNPs cut away by polyRAD", {
  clone <- c("T1", "T1", "T1", "B1", "B1")
  ids <- c("T1|Ref_0001", "T1|Alt_0002", "T1|RefMatch_0001", "B1|Ref_0001", "B1|Alt_0002")
  seqs <- c("ACGTACGTAC", "ACGTACGGAC", "ACGTA", "ACGTACGTAC", "ACGTACGGAC")
  # Top strand: index in the MADC sequence (8), even if beyond the cut (5 bp)
  # Bottom strand: counted from the end of the cut, reverse complemented Ref
  expect_equal(target_index_from_madc(clone, ids, seqs, loci = c("T1", "B1"),
                                      bottom = c(FALSE, TRUE), ref_len = c(5L, 10L)),
               c(8L, 3L))
})

test_that("SNP export helpers mask terminal gaps and anchor leading indels", {
  nuc <- c("ACGTA", "--GTA", "AC.TA", "-CGTA", "ACGTA")
  attr(nuc, "Variable_sites_only") <- FALSE
  rad <- list(alleleNucleotides = nuc, alleles2loc = c(1L, 1L, 1L, 2L, 2L),
              locTable = data.frame(Ref = c("ACGTA", "ACGTA"), Pos = c(100L, 200L)))
  m <- mask_terminal_gaps(rad)
  expect_equal(as.character(m$alleleNucleotides), c("ACGTA", "NNGTA", "AC.TA", "NCGTA", "ACGTA"))
  expect_false(attr(m$alleleNucleotides, "Variable_sites_only"))

  # An insertion column at the start ("." in Ref) needs a base before it
  rad2 <- list(alleleNucleotides = structure(c(".ACG", "TACG"), Variable_sites_only = FALSE),
               alleles2loc = c(1L, 1L), locTable = data.frame(Ref = ".ACG", Pos = 100L))
  a <- anchor_leading_gaps(rad2)
  expect_equal(as.character(a$alleleNucleotides), c("N.ACG", "NTACG"))
  expect_equal(a$locTable$Ref, "N.ACG")
  expect_equal(a$locTable$Pos, 99L)
  expect_identical(anchor_leading_gaps(m), m)  # no gap left at the start
})

test_that("add_hap_info_records names SNP records by position for asSNPs", {
  report <- data.frame(AlleleID = c("Chr01_000100|Ref_0001", "Chr01_000100|Alt_0002"),
                       CloneID = "Chr01_000100", AlleleSequence = c("AAACCCGGGT", "AAACCCTGGT"))
  hap_info <- madc_hap_lookup(report, character(0), loci = "Chr01_000100")
  hap_info$target_pos <- c(Chr01_000100 = 100L)
  lines <- c("Chr01\t95\t.\tA\tG\t.\t.\tNS=2;LU=1\tGT\t0/1",
             "Chr01\t100\t.\tG\tT\t.\t.\tNS=2;LU=1\tGT\t0/1",
             "Chr01\t101\t.\tGG\tG\t.\t.\tNS=2;LU=1\tGT\t0/1")
  f <- do.call(rbind, strsplit(add_hap_info_records(lines, hap_info)$lines, "\t"))
  expect_equal(f[, 3], c("Chr01_000095", "Chr01_000100", "Chr01_000101"))
  expect_equal(f[1, 8], "NS=2;LU=1;targetSNP=Chr01_000100;DArTstrand=+")

  # An indel written from the base before the target gets the target ID
  lines2 <- c("Chr01\t99\t.\tAG\tA\t.\t.\tNS=2;LU=1\tGT\t0/1",
              "Chr01\t105\t.\tC\tT\t.\t.\tNS=2;LU=1\tGT\t0/1")
  f2 <- do.call(rbind, strsplit(add_hap_info_records(lines2, hap_info)$lines, "\t"))
  expect_equal(f2[, 3], c("Chr01_000100", "Chr01_000105"))

  meta <- add_hap_info_header(c("##fileformat=VCFv4.3", "##INFO=<ID=NS,Number=1>"), hapID = FALSE)
  expect_equal(sub("^##INFO=<ID=([^,]+),.*", "\\1", grep("^##INFO", meta, value = TRUE)),
               c("NS", "targetSNP", "DArTstrand"))
})

test_that("madc2vcf_multi checks align_indels and asSNPs", {
  expect_error(madc2vcf_multi(test_path("nonexistent.csv"), "x", tempfile(), align_indels = NA))
  expect_error(madc2vcf_multi(test_path("nonexistent.csv"), "x", tempfile(), asSNPs = "yes"))
})
