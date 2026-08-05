context("rhAmpSeq to MADC conversion")

# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

make_test_files <- function() {
  fasta_file <- tempfile(fileext = ".fasta")
  geno_file  <- tempfile(fileext = ".txt")

  # Two loci with simple, hand-verifiable sequences (20 bp each).
  #
  # LocA: ref AAAAA…(20), allele #2 has A→C at pos 10, allele #3 has A→C at pos 15
  #   → polymorphic positions: {10, 15}
  #   → largest gap [1–9] → mid_pos = 5, A→C fake Alt
  #
  # LocB: ref TTTTT…(20), allele #2 has T→G at pos 11
  #   → polymorphic positions: {11}
  #   → largest gap [1–10] → mid_pos = 6, T→G fake Alt
  writeLines(c(
    ">LocA#1", "AAAAAAAAAAAAAAAAAAAA",
    ">LocA#2", "AAAAAAAAACAAAAAAAAAA",
    ">LocA#3", "AAAAAAAAAAAAAACAAAAA",
    ">LocB#1", "TTTTTTTTTTTTTTTTTTTT",
    ">LocB#2", "TTTTTTTTTTGTTTTTTTTT"
  ), fasta_file)

  writeLines(c(
    paste(c("Locus", "Haplotypes", "Sample1", "Sample2", "Sample3"),
          collapse = "\t"),
    paste(c("LocA", "1(0.4);2(0.35);3(0.25);",
            "1/2:50,30", "2/3:20,40", "./.:0,0"),
          collapse = "\t"),
    paste(c("LocB", "1(0.4);2(0.6);",
            "1:100",   "2:60",       "1/2:80,40"),
          collapse = "\t")
  ), geno_file)

  list(fasta = fasta_file, geno = geno_file)
}

# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------

test_that("rhampseq2madc returns a list with madc, target_positions, new_fasta", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  expect_equal(sum(out$madc$Sample1), 180)
  expect_equal(sum(out$madc[3,4:6]), 50)

  expect_type(out, "list")
  expect_named(out, c("madc", "target_positions", "new_fasta"))
  expect_s3_class(out$madc,             "data.frame")
  expect_s3_class(out$target_positions, "data.frame")
  expect_s4_class(out$new_fasta,        "DNAStringSet")
})

# ---------------------------------------------------------------------------
# MADC data.frame structure
# ---------------------------------------------------------------------------

test_that("madc has correct dimensions and column names", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  # LocA: Ref_0001 + Alt_0002 + RefMatch_0003 + RefMatch_0004 = 4 rows
  # LocB: Ref_0001 + Alt_0002 + RefMatch_0003                 = 3 rows
  expect_equal(nrow(out$madc), 7L)
  expect_true(all(c("AlleleID", "CloneID", "AlleleSequence",
                    "Sample1", "Sample2", "Sample3") %in% colnames(out$madc)))
})

test_that("AlleleID values follow Ref/Alt/RefMatch convention", {
  files    <- make_test_files()
  out      <- rhampseq2madc(files$geno, files$fasta)
  loca_ids <- out$madc$AlleleID[out$madc$CloneID == "LocA"]

  expect_equal(loca_ids[1], "LocA|Ref_0001")
  expect_equal(loca_ids[2], "LocA|Alt_0002")
  expect_equal(loca_ids[3], "LocA|RefMatch_0003")
  expect_equal(loca_ids[4], "LocA|RefMatch_0004")
})

# ---------------------------------------------------------------------------
# Read depth values
# ---------------------------------------------------------------------------

test_that("Ref allele depths are correctly extracted from genotype data", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  # LocA allele #1: appears in Sample1 (depth 50), not in Sample2 or Sample3
  ref_loca <- out$madc[out$madc$AlleleID == "LocA|Ref_0001", ]
  expect_equal(ref_loca$Sample1, 50L)
  expect_equal(ref_loca$Sample2, 0L)
  expect_equal(ref_loca$Sample3, 0L)

  # LocB allele #1: appears in Sample1 (100) and Sample3 (80)
  ref_locb <- out$madc[out$madc$AlleleID == "LocB|Ref_0001", ]
  expect_equal(ref_locb$Sample1, 100L)
  expect_equal(ref_locb$Sample2, 0L)
  expect_equal(ref_locb$Sample3, 80L)
})

test_that("RefMatch allele depths are correctly assigned", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  # LocA allele #2 → RefMatch_0003: Sample1=30, Sample2=20, Sample3=0
  rm0003 <- out$madc[out$madc$AlleleID == "LocA|RefMatch_0003", ]
  expect_equal(rm0003$Sample1, 30L)
  expect_equal(rm0003$Sample2, 20L)
  expect_equal(rm0003$Sample3, 0L)

  # LocA allele #3 → RefMatch_0004: Sample1=0, Sample2=40, Sample3=0
  rm0004 <- out$madc[out$madc$AlleleID == "LocA|RefMatch_0004", ]
  expect_equal(rm0004$Sample1, 0L)
  expect_equal(rm0004$Sample2, 40L)
  expect_equal(rm0004$Sample3, 0L)
})

test_that("Alt allele rows always have zero read depths (synthetic allele)", {
  files    <- make_test_files()
  out      <- rhampseq2madc(files$geno, files$fasta)
  alt_rows <- out$madc[grepl("\\|Alt_0002$", out$madc$AlleleID), ]
  depth_cols <- setdiff(colnames(alt_rows), c("AlleleID", "CloneID", "AlleleSequence"))

  expect_true(all(alt_rows[, depth_cols] == 0L))
})

test_that("missing genotypes (./.) contribute zero depth", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  # Sample3 had ./. for LocA — all LocA alleles should have 0 depth for Sample3
  loca_rows <- out$madc[out$madc$CloneID == "LocA", ]
  expect_true(all(loca_rows$Sample3 == 0L))
})

# ---------------------------------------------------------------------------
# Target positions
# ---------------------------------------------------------------------------

test_that("target_positions has one row per locus with correct columns", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  expect_equal(nrow(out$target_positions), 2L)
  expect_named(out$target_positions,
               c("AlleleID_rhAmpSeq", "AlleleID", "target_position", "target_base"))
})

test_that("target SNP is placed at the center of the largest polymorphism-free gap", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  # LocA: poly at {10,15}, largest gap [1-9] → mid = 5, A→C
  tp_a <- out$target_positions[out$target_positions$AlleleID_rhAmpSeq == "LocA#1", ]
  expect_equal(tp_a$target_position, 5L)
  expect_equal(tp_a$target_base,     "A/C")

  # LocB: poly at {11}, largest gap [1-10] → mid = round(5.5) = 6, T→G
  tp_b <- out$target_positions[out$target_positions$AlleleID_rhAmpSeq == "LocB#1", ]
  expect_equal(tp_b$target_position, 6L)
  expect_equal(tp_b$target_base,     "T/G")
})

# ---------------------------------------------------------------------------
# Alt allele sequence
# ---------------------------------------------------------------------------

test_that("Alt sequence differs from Ref at exactly the target position", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  for (locus in c("LocA", "LocB")) {
    ref_seq    <- as.character(out$new_fasta[paste0(locus, "|Ref_0001")])
    alt_seq    <- as.character(out$new_fasta[paste0(locus, "|Alt_0002")])
    target_pos <- out$target_positions$target_position[
      out$target_positions$AlleleID_rhAmpSeq == paste0(locus, "#1")
    ]
    diffs <- which(strsplit(ref_seq, "")[[1]] != strsplit(alt_seq, "")[[1]])
    expect_equal(length(diffs), 1L,
                 label = paste(locus, "Alt differs at exactly one position"))
    expect_equal(diffs, target_pos,
                 label = paste(locus, "Alt differs at the target position"))
  }
})

# ---------------------------------------------------------------------------
# new_fasta
# ---------------------------------------------------------------------------

test_that("new_fasta names match AlleleID column", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  expect_equal(length(out$new_fasta), nrow(out$madc))
  expect_equal(names(out$new_fasta),  out$madc$AlleleID)
})

test_that("new_fasta AlleleSequence matches madc AlleleSequence column", {
  files <- make_test_files()
  out   <- rhampseq2madc(files$geno, files$fasta)

  expect_equal(unname(as.character(out$new_fasta)), out$madc$AlleleSequence)
})

# ---------------------------------------------------------------------------
# prefix file writing
# ---------------------------------------------------------------------------

test_that("prefix writes _madc.csv, _target_positions.csv, and .fasta", {
  files <- make_test_files()
  pfx   <- file.path(tempdir(), paste0("rhampseq_test_", as.integer(Sys.time())))

  rhampseq2madc(files$geno, files$fasta, prefix = pfx)

  expect_true(file.exists(paste0(pfx, "_madc.csv")))
  expect_true(file.exists(paste0(pfx, "_target_positions.csv")))
  expect_true(file.exists(paste0(pfx, ".fasta")))
})

test_that("no files are written when prefix is NULL (default)", {
  files    <- make_test_files()
  before   <- list.files(tempdir(), full.names = TRUE)
  rhampseq2madc(files$geno, files$fasta)
  after    <- list.files(tempdir(), full.names = TRUE)
  new_files <- setdiff(after, before)

  expect_false(any(grepl("_madc\\.csv|_target_positions|new_fasta\\.fasta", new_files)))
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("error if hap_genotype_file has fewer than 3 columns", {
  geno_bad <- tempfile(fileext = ".txt")
  writeLines(c("Locus\tHaplotypes",
               "LocA\t1(1.0);"), geno_bad)

  files <- make_test_files()
  expect_error(rhampseq2madc(geno_bad, files$fasta), "at least 3 columns")
})

test_that("error if sample genotype cell has unexpected format", {
  geno_bad <- tempfile(fileext = ".txt")
  writeLines(c(
    paste(c("Locus", "Haplotypes", "Sample1"), collapse = "\t"),
    paste(c("LocA", "1(0.5);2(0.5);", "BAD"),  collapse = "\t")
  ), geno_bad)

  files <- make_test_files()
  expect_error(rhampseq2madc(geno_bad, files$fasta), "Unexpected genotype format")
})

test_that("error if FASTA sequence names do not follow LocusID#N convention", {
  fasta_bad <- tempfile(fileext = ".fasta")
  writeLines(c(">LocA_1", "AAAAAAAAAAAAAAAAAAAA"), fasta_bad)

  files <- make_test_files()
  expect_error(rhampseq2madc(files$geno, fasta_bad), "sequence names must end with '#N'")
})
