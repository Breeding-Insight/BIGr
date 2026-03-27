context("Get OffTargets")


test_that("test madc offtargets",{
  #Input variables
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  bot_file <- system.file("example_SNPs_DArTag-probe-design_f180bp.botloci", package="BIGr")
  db_file <- system.file("example_allele_db.fa", package="BIGr")

  #Calculations
  temp <- tempfile(fileext = ".vcf")
  temp_multi <- tempfile(fileext = ".vcf")

  # With hap_seq provided
  set.seed(123)
  madc2vcf_all(madc = madc_file,
               botloci_file = bot_file,
               hap_seq_file = db_file,
               n.cores = 2,
               rm_multiallelic_SNP = FALSE,
               multiallelic_SNP_dp_thr = 0,
               multiallelic_SNP_sample_thr = 0,
               alignment_score_thr = 40,
               out_vcf = temp,
               verbose = TRUE)

  set.seed(456)
  madc2vcf_all(madc = madc_file,
               botloci_file = bot_file,
               hap_seq_file = db_file,
               n.cores = 2,
               rm_multiallelic_SNP = TRUE,
               multiallelic_SNP_dp_thr = 0,
               multiallelic_SNP_sample_thr = 0,
               alignment_score_thr = 40,
               out_vcf = temp_multi,
               verbose = TRUE)

  vcf <- read.vcfR(temp)
  vcf_multi <- read.vcfR(temp_multi)

  #Check
  expect_true(all(dim(vcf@gt) == c("33","11")))
  expect_true(all(dim(vcf_multi@gt) == c("32","11")))

  rm(vcf)
  rm(vcf_multi)

  # Without hap_seq provided
  set.seed(123)
  madc2vcf_all(madc = madc_file,
               botloci_file = bot_file,
               hap_seq_file = NULL,
               n.cores = 2,
               rm_multiallelic_SNP = FALSE,
               multiallelic_SNP_dp_thr = 0,
               multiallelic_SNP_sample_thr = 0,
               out_vcf = temp,
               verbose = TRUE)

  vcf <- read.vcfR(temp)

  #Check
  expect_true(all(dim(vcf@gt) == c("28","11")))
  rm(vcf)

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

  test_that("ALFALFA — clean fixed allele ID MADC", {
    out <- tempfile(fileext = ".vcf")
    expect_no_error(
      madc2vcf_all(madc = alfalfa_madc,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 2,
                   rm_multiallelic_SNP = TRUE,
                   multiallelic_SNP_sample_thr = 0,
                   multiallelic_SNP_dp_thr = 0,
                   alignment_score_thr = 40,
                   out_vcf = out,
                   verbose = TRUE)
    )
    vcf <- read.vcfR(out, verbose = FALSE)
    expect_s4_class(vcf, "vcfR")
    expect_true(all(!is.na(vcf@fix[, "REF"])))
    expect_true(all(!is.na(vcf@fix[, "ALT"])))
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 233482)
    unlink(out)

    expect_no_error(
      madc2vcf_all(madc = alfalfa_madc,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = NULL,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )
    vcf <- read.vcfR(out, verbose = FALSE)
    expect_s4_class(vcf, "vcfR")
    expect_true(all(is.na(vcf@fix[, "REF"])))
    expect_true(all(is.na(vcf@fix[, "ALT"])))
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 56547)

    # Test error when botloci_file is NULL
    expect_error(
      madc2vcf_all(madc = alfalfa_madc,
                   botloci_file = NULL,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )

    # Test that it works when hap_seq_file is provided (REF/ALT recovered from probe sequences)
    madc2vcf_all(madc = alfalfa_madc,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = alfalfa_microhapDB,
                 n.cores = 1,
                 out_vcf = out,
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

    # Test that it also works when markers_info is provided together with botloci
    madc2vcf_all(madc = alfalfa_madc,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = alfalfa_microhapDB,
                 n.cores = 1,
                 markers_info = alfalfa_markers_info,
                 out_vcf = out,
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
    # Test error when botloci provided but no matching CloneID between botloci and MADC
    expect_error(
      madc2vcf_all(madc = alfalfa_madc_wrongID,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE),
      regexp = "Check marker IDs in both MADC and botloci files. They should be the same."
    )

    # Test error when markers_info does not match MADC CloneIDs
    expect_error(
      madc2vcf_all(madc = alfalfa_madc_wrongID,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 1,
                   markers_info = alfalfa_markers_info,
                   out_vcf = out,
                   verbose = FALSE)
    )

    # Test error when markers_info_ChromPos is provided but IDs still don't match botloci
    expect_error(
      madc2vcf_all(madc = alfalfa_madc_wrongID,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 1,
                   markers_info = alfalfa_markers_info_ChromPos,
                   out_vcf = out,
                   verbose = FALSE)
    )
  })

  test_that("alfalfa lower case fixed MADC", {
    out <- tempfile(fileext = ".vcf")
    madc2vcf_all(madc = alfalfa_lowercase,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = alfalfa_microhapDB,
                 n.cores = 1,
                 out_vcf = out,
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

    madc2vcf_all(madc = alfalfa_lowercase,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = alfalfa_microhapDB,
                 n.cores = 1,
                 markers_info = alfalfa_markers_info,
                 out_vcf = out,
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

    madc2vcf_all(madc = alfalfa_lowercase,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = alfalfa_markers_info,
                 out_vcf = out,
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
    # IUPAC codes cause a stop in madc2vcf_all
    expect_error(
      madc2vcf_all(madc = alfalfa_iupac,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )

    madc2vcf_all(madc = alfalfa_iupac,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = alfalfa_microhapDB,
                 n.cores = 1,
                 markers_info = alfalfa_markers_info,
                 out_vcf = out,
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

    madc2vcf_all(madc = alfalfa_iupac,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = alfalfa_markers_info,
                 out_vcf = out,
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

    madc2vcf_all(madc = alfalfa_iupac,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 out_vcf = out,
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
    # Indels detected, no markers_info with Ref/Alt/Indel_pos -> error
    expect_error(
      madc2vcf_all(madc = potato_indel_madc,
                   botloci_file = potato_botloci,
                   hap_seq_file = potato_microhapDB,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )

    madc2vcf_all(madc = potato_indel_madc,
                 botloci_file = potato_botloci,
                 hap_seq_file = potato_microhapDB,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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

    madc2vcf_all(madc = potato_indel_madc,
                 botloci_file = potato_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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

    # ChromPos=FALSE and no markers_info -> error
    expect_error(
      madc2vcf_all(madc = potato_indel_madc,
                   botloci_file = potato_botloci,
                   hap_seq_file = NULL,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )

    madc2vcf_all(madc = potato_indel_madc,
                 botloci_file = potato_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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
    # Indels detected, no markers_info with Ref/Alt/Indel_pos -> error
    expect_error(
      madc2vcf_all(madc = potato_more_indels_chrompos_false,
                   botloci_file = potato_botloci,
                   hap_seq_file = potato_microhapDB,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )

    madc2vcf_all(madc = potato_more_indels_chrompos_false,
                 botloci_file = potato_botloci,
                 hap_seq_file = potato_microhapDB,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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

    madc2vcf_all(madc = potato_more_indels_chrompos_false,
                 botloci_file = potato_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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

    # ChromPos=FALSE and no markers_info -> error
    expect_error(
      madc2vcf_all(madc = potato_more_indels_chrompos_false,
                   botloci_file = potato_botloci,
                   hap_seq_file = NULL,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )

    madc2vcf_all(madc = potato_more_indels_chrompos_false,
                 botloci_file = potato_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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
    madc2vcf_all(madc = potato_indel_lowercase,
                 botloci_file = potato_botloci,
                 hap_seq_file = potato_microhapDB,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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

    # markers_info without Ref/Alt/Indel_pos while indels present -> error
    expect_error(
      madc2vcf_all(madc = potato_indel_lowercase,
                   botloci_file = potato_botloci,
                   hap_seq_file = potato_microhapDB,
                   n.cores = 1,
                   markers_info = potato_markers_info_ChromPos,
                   out_vcf = out,
                   verbose = FALSE)
    )

    madc2vcf_all(madc = potato_indel_lowercase,
                 botloci_file = potato_botloci,
                 hap_seq_file = potato_microhapDB,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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

    madc2vcf_all(madc = potato_indel_lowercase,
                 botloci_file = potato_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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

    madc2vcf_all(madc = potato_indel_lowercase,
                 botloci_file = potato_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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
    madc2vcf_all(madc = potato_indel_iupac,
                 botloci_file = potato_botloci,
                 hap_seq_file = potato_microhapDB,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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

    madc2vcf_all(madc = potato_indel_iupac,
                 botloci_file = potato_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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

    madc2vcf_all(madc = potato_indel_iupac,
                 botloci_file = potato_botloci,
                 hap_seq_file = NULL,
                 n.cores = 1,
                 markers_info = potato_markers_info,
                 out_vcf = out,
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
    # Raw format fails FixAlleleIDs check -> madc2vcf_all stops
    expect_error(
      madc2vcf_all(madc = alfalfa_madc_raw,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )

    expect_error(
      madc2vcf_all(madc = alfalfa_madc_raw,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = NULL,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )

    expect_error(
      madc2vcf_all(madc = alfalfa_madc_raw,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = NULL,
                   n.cores = 1,
                   markers_info = alfalfa_markers_info,
                   out_vcf = out,
                   verbose = FALSE)
    )

    expect_error(
      madc2vcf_all(madc = alfalfa_madc_raw,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 1,
                   markers_info = alfalfa_markers_info,
                   out_vcf = out,
                   verbose = FALSE)
    )
  })
})

