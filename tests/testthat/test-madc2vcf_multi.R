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
      verbose      = TRUE
    )
  )

  vcf <- read.vcfR(out, verbose = FALSE)
  expect_s4_class(vcf, "vcfR")
  expect_equal(sum(grepl(",", vcf@fix[,5])), 281)
  GT <- extract.gt(vcf)
  expect_equal(GT[3,5],"0/0/0/3")

  # Don't allow raw MADC
  out <- tempfile(fileext = ".vcf")
  expect_error(
    madc2vcf_multi(
      madc_file    = alfalfa_madc_raw,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      ploidy       = 4L,
      verbose      = FALSE
    ), regexp = "The MADC file does not have fixed AlleleIDs. Please process the MADC file through HapApp before using this function."
  )

  out <- tempfile(fileext = ".vcf")
  expect_no_error(
    madc2vcf_multi(
      madc_file    = alfalfa_madc,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      ploidy       = 4L,
      verbose      = TRUE
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
      verbose      = TRUE
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
      verbose      = TRUE
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
      verbose      = TRUE
    ), regexp = "MADC Allele Sequences contain IUPAC \\(non-ATCG\\) codes. Please run HapApp to clean MADC file before using this function."
  )

  out <- tempfile(fileext = ".vcf")
  expect_error(
    madc2vcf_multi(
      madc_file    = alfalfa_lowercase,
      botloci_file = alfalfa_botloci,
      outfile      = out,
      ploidy       = 4L,
      verbose      = TRUE
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
      verbose      = TRUE
    )
  )

  vcf <- read.vcfR(out, verbose = FALSE)
  expect_s4_class(vcf, "vcfR")
  expect_equal(sum(grepl(",", vcf@fix[,5])), 277)
  GT <- extract.gt(vcf)
  expect_equal(GT[3,5],"0/1/1/6")

})

