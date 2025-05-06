test_that("Checking replicates",{
  example_vcf <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGr")

  check_tab <- check_replicates(path.vcf = example_vcf, select_samples = NULL)
  expect_equal(sum(check_tab$`%_matching_genotypes`), 799901)

  check_tab <- check_replicates(example_vcf, select_samples = paste0("Sample_",1:10))
  expect_equal(sum(check_tab$`%_matching_genotypes`), 3134.87, tolerance = 0.01)

})

test_that("Checking homozygous segregation by trios",{

  example_vcf <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGr")

  parents_candidates <- paste0("Sample_",1:10)
  progeny_candidates <- paste0("Sample_",11:20)

  check_tab <- check_homozygous_trios(path.vcf = example_vcf,
                                      ploidy = 2,
                                      parents_candidates = parents_candidates,
                                      progeny_candidates = progeny_candidates)

  expect_equal(sum(check_tab$homoRef_x_homoRef_match), 36562.35, tolerance = 0.01)

})
