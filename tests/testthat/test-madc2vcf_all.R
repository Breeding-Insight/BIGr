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
               verbose = FALSE)

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
               verbose = FALSE)

  vcf <- read.vcfR(temp)
  vcf_multi <- read.vcfR(temp_multi)

  #Check
  expect_true(all(dim(vcf@gt) == c("33","11")))
  expect_true(all(dim(vcf_multi@gt) == c("32","11")))

  rm(vcf)
  rm(vcf_multi)

  # Without hap_seq provided
  set.seed(123)
  expect_warning(madc2vcf_all(madc = madc_file,
                              botloci_file = bot_file,
                              hap_seq_file = NULL,
                              n.cores = 2,
                              rm_multiallelic_SNP = FALSE,
                              multiallelic_SNP_dp_thr = 0,
                              multiallelic_SNP_sample_thr = 0,
                              out_vcf = temp,
                              verbose = FALSE))

  vcf <- read.vcfR(temp)

  #Check
  expect_true(all(dim(vcf@gt) == c("28","11")))
  rm(vcf)

})

test_that("madc2vcf_all preserves comma-separated AD for biallelic targets", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  bot_file <- system.file("example_SNPs_DArTag-probe-design_f180bp.botloci", package="BIGr")
  db_file <- system.file("example_allele_db.fa", package="BIGr")
  temp <- tempfile(fileext = ".vcf")

  madc2vcf_all(madc = madc_file,
               botloci_file = bot_file,
               hap_seq_file = db_file,
               n.cores = 1,
               out_vcf = temp,
               verbose = FALSE)

  vcf <- read.vcfR(temp, verbose = FALSE)
  ad <- extract.gt(vcf, "AD")
  biallelic <- !grepl(",", vcf@fix[, "ALT"])

  expect_true(all(grepl("^[0-9]+,[0-9]+$", ad[biallelic, 1])))
})

test_that("madc2vcf_all accepts BI_markerID matches when CloneID does not match", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  bot_file <- system.file("example_SNPs_DArTag-probe-design_f180bp.botloci", package="BIGr")
  db_file <- system.file("example_allele_db.fa", package="BIGr")
  temp <- tempfile(fileext = ".vcf")
  temp_markers <- tempfile(fileext = ".csv")

  report <- read.csv(madc_file, check.names = FALSE)
  marker_ids <- unique(report$CloneID)
  markers_info <- data.frame(
    CloneID = paste0("bogus_", seq_along(marker_ids)),
    BI_markerID = marker_ids
  )
  write.csv(markers_info, temp_markers, row.names = FALSE)

  expect_no_error(
    madc2vcf_all(madc = madc_file,
                 botloci_file = bot_file,
                 hap_seq_file = db_file,
                 markers_info = temp_markers,
                 n.cores = 1,
                 out_vcf = temp,
                 verbose = FALSE)
  )
})

test_that("madc2vcf_all surfaces missing-column validation error without crashing", {
  madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
  bot_file <- system.file("example_SNPs_DArTag-probe-design_f180bp.botloci", package="BIGr")
  temp_madc <- tempfile(fileext = ".csv")

  report <- read.csv(madc_file, check.names = FALSE)
  report$AlleleSequence <- NULL
  write.csv(report, temp_madc, row.names = FALSE)

  expect_error(
    madc2vcf_all(madc = temp_madc,
                 botloci_file = bot_file,
                 out_vcf = tempfile(fileext = ".vcf"),
                 verbose = FALSE),
    "One or more required columns missing"
  )
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

  skip_on_cran()
  skip_if_offline("raw.githubusercontent.com")

  test_that("ALFALFA — clean fixed allele ID MADC", {
    out <- tempfile(fileext = ".vcf")
    #out <- "test.vcf"
    # Default parameters
    expect_no_error(
      madc2vcf_all(madc = alfalfa_madc,
                   botloci_file = alfalfa_botloci,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 2,
                   rm_multiallelic_SNP = FALSE,
                   multiallelic_SNP_sample_thr = 0,
                   multiallelic_SNP_dp_thr = 0,
                   alignment_score_thr = 40,
                   out_vcf = out,
                   verbose = FALSE)
    )
    vcf <- read.vcfR(out, verbose = FALSE)
    expect_s4_class(vcf, "vcfR")
    expect_true(all(!is.na(vcf@fix[, "REF"])))
    expect_true(all(!is.na(vcf@fix[, "ALT"])))
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 235217)
    multi <- grepl(",", vcf@fix[,5])
    expect_true(any(multi)) # It has multiallelics
    expect_equal(sum(multi), 9)
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
    expect_true(all(!is.na(vcf@fix[, "REF"])))
    expect_true(all(!is.na(vcf@fix[, "ALT"])))
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 235217)
    multi <- grepl(",", vcf@fix[,5])
    expect_true(any(multi)) # It has multiallelics
    expect_equal(sum(multi), 9)

    # Test error when botloci_file is NULL
    expect_error(
      madc2vcf_all(madc = alfalfa_madc,
                   botloci_file = NULL,
                   hap_seq_file = alfalfa_microhapDB,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE)
    )

    # Test that it also works when markers_info is provided together with botloci
    madc2vcf_all(madc = alfalfa_madc,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = alfalfa_microhapDB,
                 multiallelic_SNP_dp_thr = 80,
                 multiallelic_SNP_sample_thr = 2,
                 n.cores = 1,
                 markers_info = alfalfa_markers_info,
                 out_vcf = out,
                 verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    expect_s4_class(vcf, "vcfR")
    expect_true(all(!is.na(vcf@fix[, "REF"])))
    expect_true(all(!is.na(vcf@fix[, "ALT"])))
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 234777)
    multi <- grepl(",", vcf@fix[,5])
    expect_true(any(multi)) # It has multiallelics
    expect_equal(sum(multi), 3)

    # Remove multiallelics
    madc2vcf_all(madc = alfalfa_madc,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = alfalfa_microhapDB,
                 rm_multiallelic_SNP = TRUE,
                 n.cores = 1,
                 markers_info = alfalfa_markers_info,
                 out_vcf = out,
                 verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    expect_s4_class(vcf, "vcfR")
    expect_true(all(!is.na(vcf@fix[, "REF"])))
    expect_true(all(!is.na(vcf@fix[, "ALT"])))
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 233482)
    multi <- grepl(",", vcf@fix[,5])
    expect_false(any(multi)) # It has multiallelics
    expect_equal(sum(multi), 0)
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
                   verbose = FALSE),
      regexp = "None of the markers_info CloneID, Marker_ID, or BI_markerID values match the MADC CloneID column. Please make sure they use the same marker IDs."
    )

    # Test error when markers_info_ChromPos is provided but IDs still don't match botloci
    madc2vcf_all(madc = alfalfa_madc_wrongID,
                 botloci_file = alfalfa_botloci,
                 hap_seq_file = alfalfa_microhapDB,
                 n.cores = 1,
                 markers_info = alfalfa_markers_info_ChromPos,
                 out_vcf = out,
                 verbose = FALSE)

    vcf <- read.vcfR(out, verbose = FALSE)
    expect_s4_class(vcf, "vcfR")
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    check <- check[-which(is.na(check$Pos)),]
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    DP <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(DP[1,]), 4534)
    expect_equal(sum(DP[,5]), 235217)
    multi <- grepl(",", vcf@fix[,5])

  })

  test_that("alfalfa lower case missing 3 ref and 1 alt fixed MADC", {
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
    check <- check[-which(is.na(check$Pos)),]
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 4534)
    expect_equal(sum(dp[,5]), 233719)

    expect_warning(madc2vcf_all(madc = alfalfa_lowercase,
                                botloci_file = alfalfa_botloci,
                                hap_seq_file = NULL,
                                n.cores = 1,
                                out_vcf = out,
                                verbose = FALSE))

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    check <- check[-which(is.na(check$Pos)),]
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 4534)
    expect_equal(sum(dp[,5]), 230415)

    expect_warning(madc2vcf_all(madc = alfalfa_lowercase,
                                botloci_file = alfalfa_botloci,
                                hap_seq_file = NULL,
                                n.cores = 1,
                                markers_info = alfalfa_markers_info,
                                out_vcf = out,
                                verbose = FALSE))

    vcf <- read.vcfR(out, verbose = FALSE)
    lut <- read.csv(alfalfa_markers_info)
    vcf_infos <- vcf@fix[,c(1:5)]
    lut_infos <- lut[match(vcf@fix[,3],lut$BI_markerID),c(2:6)]
    check <- cbind(vcf_infos,lut_infos)
    check <- check[-which(is.na(check$Pos)),]
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 4534)
    expect_equal(sum(dp[,5]), 230415)

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
                   verbose = FALSE),
      regexp = "IUPAC \\(non-ATCG\\) codes found in Ref_/Alt_ AlleleSequence that differ between alleles or are present in only one. These codes are not currently supported by BIGr/BIGapp"
    )
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
                   verbose = FALSE),
      regexp = "CloneID does not have the expected Chromosome_Position format. Please check your CloneIDs or provide a file with this information"

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
    check <- check[-which(is.na(check$Ref)),]
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 226838)
    expect_equal(sum(dp[3,]), 3996)
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
    check <- check[-which(is.na(check$Ref)),]
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 3937)
    expect_equal(sum(dp[,5]), 248571)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    # ChromPos=FALSE and no markers_info -> error
    expect_error(
      madc2vcf_all(madc = potato_indel_madc,
                   botloci_file = potato_botloci,
                   hap_seq_file = NULL,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE),
      regexp = "CloneID does not have the expected Chromosome_Position format. Please check your CloneIDs or provide a file with this information"
    )

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
                   verbose = FALSE),
      regexp = "CloneID does not have the expected Chromosome_Position format. Please check your CloneIDs or provide a file with this information"
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
    check <- check[-which(is.na(check$Ref)),]
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 5397)
    expect_equal(sum(dp[,5]), 215070)
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
    check <- check[-which(is.na(check$Ref)),]
    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[1,]), 5397)
    expect_equal(sum(dp[,5]), 215070)
    expect_equal(check$REF, check$Ref)
    expect_equal(check$ALT, check$Alt)

    # ChromPos=FALSE and no markers_info -> error
    expect_error(
      madc2vcf_all(madc = potato_more_indels_chrompos_false,
                   botloci_file = potato_botloci,
                   hap_seq_file = NULL,
                   n.cores = 1,
                   out_vcf = out,
                   verbose = FALSE),
      regexp = "CloneID does not have the expected Chromosome_Position format. Please check your CloneIDs or provide a file with this information"
    )
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
    check <- check[-which(is.na(check$Ref)),]

    expect_equal(as.numeric(check$POS), check$Pos)
    dp <- extract.gt(vcf, "DP", as.numeric = TRUE)
    expect_equal(sum(dp[,10]), 219742)
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
                   verbose = FALSE),
      regexp = "Indels detected in MADC file. The markers_info file must contain 'Ref', 'Alt', and 'Indel_pos' columns."
    )
  })

  test_that("potato IUPAC", {
    out <- tempfile(fileext = ".vcf")

    expect_error(
      madc2vcf_all(madc = potato_indel_iupac,
                   botloci_file = potato_botloci,
                   hap_seq_file = potato_microhapDB,
                   n.cores = 1,
                   markers_info = potato_markers_info,
                   out_vcf = out,
                   verbose = FALSE),
      regexp = "IUPAC \\(non-ATCG\\) codes found in Ref_/Alt_ AlleleSequence that differ between alleles or are present in only one. These codes are not currently supported by BIGr/BIGapp"
    )
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
                   verbose = FALSE), regexp = "MADC not processed by HapApp"
    )
  })
})

# =======================================================================
# Other allele filters and IUPAC codes (internal compare function)
# =======================================================================

test_that("Other alleles with too many close SNPs are discarded", {

  clone <- "Chr01_000000100"
  make_tag <- function(seqs, ids) {
    data.frame(AlleleID = paste0(clone, "|", ids), Chromosome = NA, SNP_position_in_Genome = NA,
               Ref = NA, Alt = NA, CloneID = clone, AlleleSequence = seqs,
               S1 = 10, S2 = 5, check.names = FALSE)
  }
  sub_bases <- function(seq, pos, bases) {
    x <- strsplit(seq, "")[[1]]
    x[pos] <- bases
    paste(x, collapse = "")
  }

  # Target SNP at position 30 of the sequence (position 100 in the genome)
  ref <- "ACGTTGCAAGCTTGACCTGAAGTCCATGGATCCGTACGATCGGATCCAAGTTCGACTGAT"
  alt <- sub_bases(ref, 30, "G")
  # Other_0003 has four SNPs 2bp apart, Other_0004 has three
  other4 <- sub_bases(ref, c(10, 12, 14, 16), c("C", "A", "A", "T"))
  other3 <- sub_bases(ref, c(10, 12, 14), c("C", "A", "A"))
  tag <- make_tag(c(ref, alt, other4, other3), c("Ref_0001", "Alt_0002", "Other_0003", "Other_0004"))
  botloci <- data.frame(V1 = "Chr02_000000500")

  # Defaults (others_min_dist = 5, others_max_close_snps = 3): only Other_0003 is discarded
  res <- compare(tag, botloci)
  expect_equal(res$n_rm_others_close_snps, 1L)
  expect_false(any(grepl("Other_0003", res$update_tag$AlleleID)))
  expect_equal(sum(grepl("Other_0004", res$update_tag$AlleleID)), 3)
  expect_equal(res$update_tag$SNP_position_in_Genome[grepl("Other_0004", res$update_tag$AlleleID)],
               c("80", "82", "84"))

  # The filter is not applied when others_max_close_snps is NULL
  res <- compare(tag, botloci, others_max_close_snps = NULL)
  expect_equal(res$n_rm_others_close_snps, 0L)
  expect_equal(sum(grepl("Other_0003", res$update_tag$AlleleID)), 4)

  # SNPs 2bp apart are not close when others_min_dist = 2
  res <- compare(tag, botloci, others_min_dist = 2)
  expect_equal(res$n_rm_others_close_snps, 0L)
  expect_equal(sum(grepl("Other_0003", res$update_tag$AlleleID)), 4)
})

test_that("identical IUPAC codes in Ref and Alt are not reported as SNPs", {

  clone <- "Chr01_000000100"
  ref <- "ACGTTGCAAGCTTGACCTGAAGTCCATGGATCCGTACGATCGGATNCAAGTTCGACTGAT"
  alt <- "ACGTTGCAAGCTTGACCTGAAGTCCATGGGTCCGTACGATCGGATNCAAGTTCGACTGAT"
  # RefMatch with one SNP at position 10 and the same N at position 46
  ref_match <- "ACGTTGCAACCTTGACCTGAAGTCCATGGATCCGTACGATCGGATNCAAGTTCGACTGAT"
  tag <- data.frame(AlleleID = paste0(clone, c("|Ref_0001", "|Alt_0002", "|RefMatch_0003")),
                    Chromosome = NA, SNP_position_in_Genome = NA, Ref = NA, Alt = NA,
                    CloneID = clone, AlleleSequence = c(ref, alt, ref_match),
                    S1 = 10, S2 = 5, check.names = FALSE)
  botloci <- data.frame(V1 = "Chr02_000000500")

  res <- compare(tag, botloci)
  # The tag is kept, with the target SNP and the RefMatch SNP only
  expect_null(res$rm_N)
  expect_null(res$rm_indels)
  expect_equal(nrow(res$update_tag), 3)
  expect_equal(res$update_tag$SNP_position_in_Genome, c("100", "100", "80"))
  expect_false(any(res$update_tag$Ref == "N" | res$update_tag$Alt == "N"))
})
