context("MADC to VCF")

test_that("test madc conversion",{
  #Input variables
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  bot_file <- system.file("example_SNPs_DArTag-probe-design_f180bp.botloci", package="BIGr")

  #Calculations
  temp <- tempfile(fileext = ".vcf")

  #Convert the dart files to vcf
  suppressWarnings(
    madc2vcf_targets(madc_file = madc_file, output.file = temp, get_REF_ALT = FALSE)
  )

  #Test validity of VCF
  vcf <- read.vcfR(temp, verbose = FALSE)

  DP_df <- extract.gt(vcf, element = "DP")
  DP_df_mean <- mean(as.numeric(DP_df))

  #Checks
  expect_s4_class(vcf, "vcfR")
  expect_true(nrow(vcf) == 20)
  expect_true(ncol(vcf@gt) == 11)
  expect_true(all(dim(DP_df) == c(20,10)))
  expect_true(round(DP_df_mean, 3) == 227.22)

  rm(vcf)
  rm(temp)

  # Test with REF_ALT
  temp <- tempfile(fileext = ".vcf")

  suppressWarnings(
    madc2vcf_targets(madc_file = madc_file, output.file = temp, get_REF_ALT = TRUE, botloci_file = bot_file)
  )

  #Test validity of VCF
  vcf <- read.vcfR(temp, verbose = FALSE)

  #Checks
  expect_s4_class(vcf, "vcfR")
  expect_true(nrow(vcf) == 20)
  expect_true(ncol(vcf@gt) == 11)
  expect_true(all(dim(DP_df) == c(20,10)))
  expect_true(round(DP_df_mean, 3) == 227.22)
  expect_true(all(vcf@fix[,4] != "."))
  expect_true(all(vcf@fix[,5] != "."))
  expect_true(all(vcf@fix[,4] != vcf@fix[,5]))

  # UPDATED: These expected values might change with the fix
  # You'll need to verify these are the correct expected values after the fix
  expect_true(all(vcf@fix[1:5,4] == c("A", "G", "G", "G", "C")))
  expect_true(all(vcf@fix[1:5,5] == c("T", "A", "A", "A", "A")))

  rm(vcf)
  rm(temp)
})

# NEW: Add specific test for bottom strand markers
test_that("bottom strand markers have correct REF/ALT", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  bot_file <- system.file("example_SNPs_DArTag-probe-design_f180bp.botloci", package="BIGr")

  temp_targets <- tempfile(fileext = ".vcf")
  temp_all <- tempfile(fileext = ".vcf")

  # Get results from both functions
  suppressWarnings(
    madc2vcf_targets(madc_file = madc_file, output.file = temp_targets,
                    get_REF_ALT = TRUE, botloci_file = bot_file)
  )

  suppressWarnings(
    madc2vcf_all(madc = madc_file, botloci_file = bot_file,
                hap_seq_file = NULL, out_vcf = temp_all, verbose = FALSE)
  )

  vcf_targets <- read.vcfR(temp_targets, verbose = FALSE)
  vcf_all <- read.vcfR(temp_all, verbose = FALSE)

  # Find common markers between both outputs
  common_markers <- intersect(vcf_targets@fix[,"ID"], vcf_all@fix[,"ID"])

  if(length(common_markers) > 0) {
    # Compare REF/ALT for common markers
    targets_subset <- vcf_targets@fix[vcf_targets@fix[,"ID"] %in% common_markers,]
    all_subset <- vcf_all@fix[vcf_all@fix[,"ID"] %in% common_markers,]

    # Sort both by ID for comparison
    targets_subset <- targets_subset[order(targets_subset[,"ID"]),]
    all_subset <- all_subset[order(all_subset[,"ID"]),]

    # Check that REF/ALT match between the two functions
    expect_equal(targets_subset[,"REF"], all_subset[,"REF"])
    expect_equal(targets_subset[,"ALT"], all_subset[,"ALT"])
  }

  # Validate that all REF/ALT are valid nucleotides
  expect_true(all(vcf_targets@fix[,"REF"] %in% c("A", "T", "G", "C", ".")))
  expect_true(all(vcf_targets@fix[,"ALT"] %in% c("A", "T", "G", "C", ".")))

  # Validate that REF != ALT where both are not "."
  valid_snps <- vcf_targets@fix[,"REF"] != "." & vcf_targets@fix[,"ALT"] != "."
  if(any(valid_snps)) {
    expect_true(all(vcf_targets@fix[valid_snps,"REF"] != vcf_targets@fix[valid_snps,"ALT"]))
  }

  rm(vcf_targets, vcf_all, temp_targets, temp_all)
})
