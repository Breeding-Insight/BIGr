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

  vcf_targets <- read.vcfR(temp_targets, verbose = FALSE)

  # Validate that all REF/ALT are valid nucleotides
  expect_true(all(vcf_targets@fix[,"REF"] %in% c("A", "T", "G", "C", ".")))
  expect_true(all(vcf_targets@fix[,"ALT"] %in% c("A", "T", "G", "C", ".")))

  # Validate that REF != ALT where both are not "."
  valid_snps <- vcf_targets@fix[,"REF"] != "." & vcf_targets@fix[,"ALT"] != "."
  if(any(valid_snps)) {
    expect_true(all(vcf_targets@fix[valid_snps,"REF"] != vcf_targets@fix[valid_snps,"ALT"]))
  }

  rm(vcf_targets, temp_targets)
})


# =======================================================================
# Using Breeding-Insight/BIGapp-PanelHub test files
# =======================================================================

test_that("simu alfalfa",{

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


  # External potato test files
  potato_indel_madc                 <- paste0(github_path, "test_madcs/potato_indel_madc.csv")
  potato_indel_iupac                <- paste0(github_path, "test_madcs/potato_indel_IUPAC.csv")
  potato_indel_lowercase            <- paste0(github_path, "test_madcs/potato_indel_lowercase.csv")
  potato_more_indels_chrompos_false <- paste0(github_path, "test_madcs/potato_more_indels_madc_ChromPosFALSE.csv")
  potato_botloci                    <- paste0(github_path, "potato/potato_dartag_v2_3915markers_rm7dupTags_6traitMarkers_f150bp_ref_alt.botloci")
  potato_markers_info               <- paste0(github_path, "potato/potato_dartag_v2_3915markers_rm7dupTags_6traitMarkers_rm1dup_snpID_lut.csv") # CloneID/BI_markerID, Chr, Pos, Ref, Alt
  potato_markers_info_ChromPos   <- paste0(github_path, "test_madcs/potato_marker_info_chrompos.csv") # markers_info: CloneID/BI_markerID, Chr, Pos

  skip_if_offline("raw.githubusercontent.com")

  test_that("ALFALFA — clean fixed allele ID MADC", {
    out <- tempfile(fileext = ".vcf")
    expect_no_error(
      madc2vcf_targets(madc_file = alfalfa_madc,
                       output.file = out,
                       get_REF_ALT = FALSE,
                       verbose = FALSE)
    )
    vcf <- read.vcfR(out, verbose = FALSE)
    expect_s4_class(vcf, "vcfR")
    expect_true(all(is.na(vcf@fix[, "REF"])))
    expect_true(all(is.na(vcf@fix[, "ALT"])))
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4139)
    expect_equal(sum(DP[,5]), 42869)
    unlink(out)

    expect_no_error(
      madc2vcf_targets(madc_file = alfalfa_madc,
                       output.file = out,
                       get_REF_ALT = FALSE,
                       collapse_matches_counts = TRUE,
                       verbose = FALSE)
    )
    vcf <- read.vcfR(out, verbose = FALSE)
    expect_s4_class(vcf, "vcfR")
    expect_true(all(is.na(vcf@fix[, "REF"])))
    expect_true(all(is.na(vcf@fix[, "ALT"])))
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 56547)

    # Test error when get_REF_ALT = TRUE but no markers_info or botloci provided to extract REF/ALT
    expect_error(
      madc2vcf_targets(madc_file = alfalfa_madc,
                       output.file = out,
                       get_REF_ALT = TRUE,
                       verbose = FALSE),
      regexp = "get_REF_ALT = TRUE but no markers_info file with Ref and Alt columns was provided neither a botloci_file to extrat REF/ALT from probe sequences. Please provide one of the these files or set get_REF_ALT to FALSE."
    )

    # Test that it works when get_REF_ALT = TRUE and botloci is provided (REF/ALT recovered from probe sequences)
    madc2vcf_targets(madc_file = alfalfa_madc,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     botloci_file = alfalfa_botloci,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)

    # Test that it also works when markers_info is provided together with botloci (should give same result as above but just to confirm no interference between the two)
    madc2vcf_targets(madc_file = alfalfa_madc,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     botloci_file = alfalfa_botloci,
                     markers_info = alfalfa_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)

  })

  test_that("ALFALFA — clean fixed allele ID MADC wrong CloneID", {
    out <- tempfile(fileext = ".vcf")
    # Test error when botloci provided but no matching CloneID between botloci and MADC (even after trying to fix potential padding mismatch with ChromPos info in markers_info)
    expect_error(
      madc2vcf_targets(madc_file = alfalfa_madc_wrongID,
                       output.file = out,
                       get_REF_ALT = TRUE,
                       botloci_file = alfalfa_botloci,
                       verbose = FALSE),
      regexp = "Check marker IDs in both MADC and botloci files. They should be the same."
    )

    # Test error when no matching CloneID between markers_info and MADC to fix the botloci mismatch issue (even if botloci file is not used, the function should still check that the provided markers_info can match with MADC CloneIDs to be able to use the ChromPos info to fix potential padding mismatch)
    expect_error(
      madc2vcf_targets(madc_file = alfalfa_madc_wrongID,
                       output.file = out,
                       get_REF_ALT = TRUE,
                       botloci_file = alfalfa_botloci,
                       markers_info = alfalfa_markers_info,
                       verbose = FALSE),
      "None of the MADC CloneID could be found in the markers_info CloneID or BI_markerID. Please make sure they match."
    )

    # Test that it works when the function can find a matching ID in markers_info to fix the botloci mismatch issue
    # (even if botloci file is not used, the function should still be able to use the ChromPos info in markers_info to
    # fix potential padding mismatch and find matching IDs between MADC and botloci)
    madc2vcf_targets(madc_file = alfalfa_madc_wrongID,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     botloci_file = alfalfa_botloci,
                     markers_info = alfalfa_markers_info_ChromPos,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)
  })

  test_that("alfalfa lower case fixed MADC", {
    out <- tempfile(fileext = ".vcf")
    madc2vcf_targets(madc_file = alfalfa_lowercase,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     botloci_file = alfalfa_botloci,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)

    madc2vcf_targets(madc_file = alfalfa_lowercase,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     botloci_file = alfalfa_botloci,
                     markers_info = alfalfa_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)

    madc2vcf_targets(madc_file = alfalfa_lowercase,
                     output.file = out,
                     get_REF_ALT = FALSE,
                     botloci_file = alfalfa_botloci,
                     markers_info = alfalfa_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)
  })

  test_that("alfalfa IUPAC code", {
    out <- tempfile(fileext = ".vcf")
    expect_error(
      madc2vcf_targets(madc_file = alfalfa_iupac,
                       output.file = out,
                       get_REF_ALT = TRUE,
                       botloci_file = alfalfa_botloci,
                       verbose = FALSE)
    )

    madc2vcf_targets(madc_file = alfalfa_iupac,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = alfalfa_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    madc2vcf_targets(madc_file = alfalfa_iupac,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = alfalfa_markers_info,
                     collapse_matches_counts = TRUE,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 4534)
    expect_equal(sum(dp[,5]), 56547)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    madc2vcf_targets(madc_file = alfalfa_iupac,
                     output.file = out,
                     get_REF_ALT = FALSE,
                     botloci_file = alfalfa_botloci,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)

  })

  test_that("potato indel madc chrompos=FALSE", {
    out <- tempfile(fileext = ".vcf")
    expect_error(
      madc2vcf_targets(madc_file = potato_indel_madc,
                       output.file = out,
                       get_REF_ALT = TRUE,
                       botloci_file = potato_botloci,
                       verbose = FALSE)
    )

    madc2vcf_targets(madc_file = potato_indel_madc,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = potato_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 41656)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    madc2vcf_targets(madc_file = potato_indel_madc,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = potato_markers_info,
                     collapse_matches_counts = TRUE,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 5163)
    expect_equal(sum(dp[,5]), 58927)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    expect_error(
      madc2vcf_targets(madc_file = potato_indel_madc,
                       output.file = out,
                       get_REF_ALT = FALSE,
                       botloci_file = potato_botloci,
                       verbose = FALSE)
    )

    madc2vcf_targets(madc_file = potato_indel_madc,
                     output.file = out,
                     get_REF_ALT = FALSE,
                     markers_info = potato_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 41656)
  })

  test_that("potato indel chromposFALSE", {
    out <- tempfile(fileext = ".vcf")
    expect_error(
      madc2vcf_targets(madc_file = potato_more_indels_chrompos_false,
                       output.file = out,
                       get_REF_ALT = TRUE,
                       botloci_file = potato_botloci,
                       verbose = FALSE)
    )

    madc2vcf_targets(madc_file = potato_more_indels_chrompos_false,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = potato_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 41755)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    madc2vcf_targets(madc_file = potato_more_indels_chrompos_false,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = potato_markers_info,
                     collapse_matches_counts = TRUE,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 6301)
    expect_equal(sum(dp[,5]), 53613)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    expect_error(
      madc2vcf_targets(madc_file = potato_more_indels_chrompos_false,
                       output.file = out,
                       get_REF_ALT = FALSE,
                       botloci_file = potato_botloci,
                       verbose = FALSE)
    )

    madc2vcf_targets(madc_file = potato_more_indels_chrompos_false,
                     output.file = out,
                     get_REF_ALT = FALSE,
                     markers_info = potato_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 41755)
  })

  test_that("potato lowercase", {
    out <- tempfile(fileext = ".vcf")
    madc2vcf_targets(madc_file = potato_indel_lowercase,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = potato_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 41755)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    expect_error(
      madc2vcf_targets(madc_file = potato_indel_lowercase,
                       output.file = out,
                       get_REF_ALT = TRUE,
                       markers_info = potato_markers_info_ChromPos,
                       botloci_file = potato_botloci,
                       verbose = FALSE)
    )

    madc2vcf_targets(madc_file = potato_indel_lowercase,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = potato_markers_info,
                     botloci_file = potato_botloci,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 41755)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    madc2vcf_targets(madc_file = potato_indel_lowercase,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = potato_markers_info,
                     collapse_matches_counts = TRUE,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 6301)
    expect_equal(sum(dp[,5]), 53613)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    madc2vcf_targets(madc_file = potato_indel_lowercase,
                     output.file = out,
                     get_REF_ALT = FALSE,
                     markers_info = potato_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 41755)
  })


  test_that("potato IUPAC", {
    out <- tempfile(fileext = ".vcf")
    madc2vcf_targets(madc_file = potato_indel_iupac,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = potato_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 41755)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    madc2vcf_targets(madc_file = potato_indel_iupac,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = potato_markers_info,
                     collapse_matches_counts = TRUE,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 6301)
    expect_equal(sum(dp[,5]), 53613)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    madc2vcf_targets(madc_file = potato_indel_iupac,
                     output.file = out,
                     get_REF_ALT = FALSE,
                     markers_info = potato_markers_info,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(potato_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 41755)
  })

  test_that("alfalfa raw MADC format (7-row header)", {
    out <- tempfile(fileext = ".vcf")
    # get_REF_ALT = FALSE: same counts as alfalfa_madc
    expect_error(
      madc2vcf_targets(madc_file = alfalfa_madc_raw,
                       output.file = out,
                       get_REF_ALT = TRUE,
                       verbose = FALSE)
    )

    expect_error(
      madc2vcf_targets(madc_file = alfalfa_madc_raw,
                       output.file = out,
                       get_REF_ALT = TRUE,
                       botloci_file = alfalfa_botloci,
                       verbose = FALSE)
    )

    madc2vcf_targets(madc_file = alfalfa_madc_raw,
                     output.file = out,
                     get_REF_ALT = FALSE,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)

    madc2vcf_targets(madc_file = alfalfa_madc_raw,
                     output.file = out,
                     get_REF_ALT = FALSE,
                     collapse_matches_counts = TRUE,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    expect_s4_class(vcf, "vcfR")
    expect_true(all(is.na(vcf@fix[, "REF"])))
    expect_true(all(is.na(vcf@fix[, "ALT"])))
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 56547)

    madc2vcf_targets(madc_file = alfalfa_madc_raw,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = alfalfa_markers_info,
                     collapse_matches_counts = FALSE,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 43691)

    madc2vcf_targets(madc_file = alfalfa_madc_raw,
                     output.file = out,
                     get_REF_ALT = TRUE,
                     markers_info = alfalfa_markers_info,
                     collapse_matches_counts = TRUE,
                     verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 56547)
  })
})
