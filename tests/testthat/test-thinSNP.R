context("Thinning")


test_that("Thinning SNPs",{

  ##' # Create sample SNP data
  set.seed(123)
  n_snps <- 20
  snp_data <- data.frame(
      MarkerID = paste0("SNP", 1:n_snps),
      Chrom = sample(c("chr1", "chr2"), n_snps, replace = TRUE),
      ChromPosPhysical = c(
        sort(sample(1:1000, 5)), # SNPs on chr1
        sort(sample(1:1000, 5)) + 500, # More SNPs on chr1
        sort(sample(1:2000, 10))      # SNPs on chr2
      ),
      Allele = sample(c("A/T", "G/C"), n_snps, replace = TRUE)
      )
  # Ensure it's sorted by Chrom and ChromPosPhysical
  snp_data <- snp_data[order(snp_data$Chrom, snp_data$ChromPosPhysical), ]
  rownames(snp_data) <- NULL

  # Thin the SNPs, keeping a minimum distance of 100 units (e.g., bp)
  thinned_snps <- thinSNP(
                    df = snp_data,
                    chrom_col_name = "Chrom",
                    pos_col_name = "ChromPosPhysical",
                    min_distance = 100
                  )

  # Thin with a larger distance
  thinned_snps_large_dist <- thinSNP(
                                df = snp_data,
                                chrom_col_name = "Chrom",
                                pos_col_name = "ChromPosPhysical",
                                min_distance = 500
                              )

  # Check the results
  expect_true(all(dim(thinned_snps) == c(14, 4))) # Expect 14 SNPs
  expect_true(all(dim(thinned_snps_large_dist) == c(5, 4))) # Expect 5 SNPs
  expect_equal(ncol(thinned_snps), ncol(snp_data)) # Same number of columns


})
