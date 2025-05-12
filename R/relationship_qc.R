#' Check Homozygous Loci in Trios
#'
#' This function analyzes homozygous loci segregation in trios (parents and progeny) using genotype data from a VCF file.
#' It calculates the percentage of homozygous loci in the progeny that match the expected segregation patterns based on the tested parents.
#'
#' @param path.vcf A string specifying the path to the VCF file containing genotype data.
#' @param ploidy An integer specifying the ploidy level of the samples. Default is 4.
#' @param parents_candidates A character vector of parent sample names to be tested. Must be provided.
#' @param progeny_candidates A character vector of progeny sample names to be tested. Must be provided.
#' @param verbose A logical value indicating whether to print the number of combinations tested. Default is TRUE.
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item \code{parent1}: The name of the first parent in the pair.
#'     \item \code{parent2}: The name of the second parent in the pair.
#'     \item \code{progeny}: The name of the progeny sample.
#'     \item \code{homoRef_x_homoRef_n}: Number of loci where both parents are homozygous reference.
#'     \item \code{homoRef_x_homoRef_match}: Percentage of matching loci in the progeny for homozygous reference parents.
#'     \item \code{homoAlt_x_homoAlt_n}: Number of loci where both parents are homozygous alternate.
#'     \item \code{homoAlt_x_homoAlt_match}: Percentage of matching loci in the progeny for homozygous alternate parents.
#'     \item \code{homoRef_x_homoAlt_n}: Number of loci where one parent is homozygous reference and the other is homozygous alternate.
#'     \item \code{homoRef_x_homoAlt_match}: Percentage of matching loci in the progeny for mixed homozygous parents.
#'     \item \code{homoalt_x_homoRef_n}: Number of loci where one parent is homozygous alternate and the other is homozygous reference.
#'     \item \code{homoalt_x_homoRef_match}: Percentage of matching loci in the progeny for mixed homozygous parents (alternate-reference).
#'     \item \code{missing}: The number of loci with missing genotype data in the comparison.
#'   }
#'
#' @details This function is designed to validate the segregation of homozygous loci in trios, ensuring that the progeny genotypes align with the expected patterns based on the parental genotypes. It requires both parent and progeny candidates to be specified. The function validates the ploidy level and ensures that all specified samples are present in the VCF file. The results include detailed statistics for each combination of parents and progeny. Reciprocal comparisons (e.g., A vs. B and B vs. A) and self-comparisons (e.g., A vs. A) are removed to avoid redundancy. Missing genotype data is also accounted for and reported in the results.
#'
#' @importFrom vcfR read.vcfR extract.gt
#'
#' @export
check_homozygous_trios <- function(path.vcf, ploidy = 4, parents_candidates = NULL, progeny_candidates = NULL, verbose = TRUE) {

  # Check if parents and progeny are not NULL
  if (is.null(parents_candidates) || is.null(progeny_candidates)) {
    stop("Please provide both parents and progeny candidates.")
  }

  # Load the VCF file
  vcf <- read.vcfR(path.vcf, verbose = FALSE)

  # Extract the genotype data
  GT <- extract.gt(vcf, element = "GT")

  count_ones <- function(x) {
    sapply(strsplit(x, "/"), function(v) sum(v == "1"))
  }

  # Apply to entire matrix
  GT_counts <- apply(GT, c(1, 2), count_ones)
  if (max(GT_counts) != ploidy) stop("Ploidy level is not correct, check the VCF file")

  if (!all(parents_candidates %in% colnames(GT))) stop("Some parents are not in the VCF file")
  if (!all(progeny_candidates %in% colnames(GT))) stop("Some progeny are not in the VCF file")

  # Get all combinations of samples
  combinations <- expand.grid(parents_candidates, parents_candidates)
  filtered_combinations <- combinations[!duplicated(t(apply(combinations, 1, sort))), ] # remove reciprocal combinations
  filtered_combinations <- filtered_combinations[filtered_combinations$Var1 != filtered_combinations$Var2, ] # remove self-comparisons

  progeny_rep <- rep(progeny_candidates, each = nrow(filtered_combinations))
  rownames(filtered_combinations) <- NULL
  filtered_combinations <- cbind(filtered_combinations, progeny_rep)
  colnames(filtered_combinations) <- c("parent1", "parent2", "progeny")
  if(verbose) cat("Number of combinations tested: ", nrow(filtered_combinations), "\n")

  # Initialize a data frame to store results
  homo1 <- 0
  homo2 <- ploidy

  matches <- mapply(function(parent1, parent2, progeny, homo1, homo2) {
    # Get the genotypes for the two parents
    gt1 <- GT_counts[, parent1]
    gt2 <- GT_counts[, parent2]
    gt3 <- GT_counts[, progeny]

    # Check if the genotypes are compatible
    homoRef_x_homoRef <- which(gt1 == homo1 & gt2 == homo1)
    homoAlt_x_homoAlt <- which(gt1 == homo2 & gt2 == homo2)
    homoRef_x_homoAlt <- which(gt1 == homo1 & gt2 == homo2)
    homoalt_x_homoRef <- which(gt1 == homo2 & gt2 == homo1)
    homoRef_x_homoRef_n <- length(homoRef_x_homoRef)
    homoAlt_x_homoAlt_n <- length(homoAlt_x_homoAlt)
    homoRef_x_homoAlt_n <- length(homoRef_x_homoAlt)
    homoalt_x_homoRef_n <- length(homoalt_x_homoRef)
    homoRef_x_homoRef_match <- ifelse(homoRef_x_homoRef_n > 0, round((sum(gt3[homoRef_x_homoRef] == homo1, na.rm = TRUE) / homoRef_x_homoRef_n) * 100, 2), NA)
    homoAlt_x_homoAlt_match <- ifelse(homoAlt_x_homoAlt_n > 0, round((sum(gt3[homoAlt_x_homoAlt] == homo2, na.rm = TRUE) / homoAlt_x_homoAlt_n) * 100, 2), NA)
    homoRef_x_homoAlt_match <- ifelse(homoRef_x_homoAlt_n > 0, round((sum(gt3[homoRef_x_homoAlt] == ploidy / 2, na.rm = TRUE) / homoRef_x_homoAlt_n) * 100, 2), NA)
    homoalt_x_homoRef_match <- ifelse(homoalt_x_homoRef_n > 0, round((sum(gt3[homoalt_x_homoRef] == ploidy / 2, na.rm = TRUE) / homoalt_x_homoRef_n) * 100, 2), NA)

    miss <- sum(is.na(gt1) | is.na(gt2) | is.na(gt3))

    # Return the result
    result <- c(
      homoRef_x_homoRef_n, homoRef_x_homoRef_match, homoAlt_x_homoAlt_n,
      homoAlt_x_homoAlt_match, homoRef_x_homoAlt_n, homoRef_x_homoAlt_match,
      homoalt_x_homoRef_n, homoalt_x_homoRef_match, miss
    )
    return(result)
  }, filtered_combinations$parent1, filtered_combinations$parent2, filtered_combinations$progeny, homo1, homo2)

  all_comb <- cbind(filtered_combinations, t(matches))
  head(all_comb)
  colnames(all_comb) <- c(
    "parent1", "parent2", "progeny", "homoRef_x_homoRef_n", "homoRef_x_homoRef_match",
    "homoAlt_x_homoAlt_n", "homoAlt_x_homoAlt_match", "homoRef_x_homoAlt_n",
    "homoRef_x_homoAlt_match", "homoalt_x_homoRef_n", "homoalt_x_homoRef_match", "missing"
  )

  return(all_comb)
}

#' Compatibility Between Samples Genotypes
#'
#' This function checks the compatibility between sample genotypes in a VCF file by comparing all pairs of samples.
#'
#' @param path.vcf A string specifying the path to the VCF file containing genotype data.
#' @param select_samples An optional character vector of sample names to be selected for comparison. If NULL (default), all samples in the VCF file are used.
#' @param verbose A logical value indicating whether to print the number of combinations tested. Default is TRUE.
#'
#' @return A data frame with four columns:
#'   \itemize{
#'     \item \code{sample1}: The name of the first sample in the pair.
#'     \item \code{sample2}: The name of the second sample in the pair.
#'     \item \code{\%_matching_genotypes}: The percentage of compatible genotypes between the two samples.
#'     \item \code{\%_missing_genotypes}: The percentage of missing genotypes in the comparison.
#'   }
#'
#' @details The function removes reciprocal comparisons (e.g., A vs. B and B vs. A) and self-comparisons (e.g., A vs. A) to avoid redundancy. Compatibility is calculated as the percentage of matching genotypes between two samples, excluding missing values. The percentage of missing genotypes is also reported for each pair.
#'
#' @importFrom vcfR read.vcfR extract.gt
#'
#' @export
check_replicates <- function(path.vcf, select_samples = NULL, verbose = TRUE) {
  # Load the VCF file
  vcf <- read.vcfR(path.vcf, verbose = FALSE)

  # Extract the genotype data
  GT <- extract.gt(vcf, element = "GT", convertNA = TRUE)
  if(any(GT == "./.")) GT[which(GT == "./.")] <- NA

  # Select samples
  if (is.null(select_samples)) {
    samples <- colnames(GT)
  } else {
    if (!all(select_samples %in% colnames(GT))) stop("Some samples are not in the VCF file")
    samples <- select_samples
  }

  # Get all combinations of samples
  combinations <- expand.grid(samples, samples)
  filtered_combinations <- combinations[!duplicated(t(apply(combinations, 1, sort))), ] # remove reciprocal combinations
  filtered_combinations <- filtered_combinations[filtered_combinations$Var1 != filtered_combinations$Var2, ] # remove self-comparisons

  if(verbose) cat("Number of combinations tested: ", nrow(filtered_combinations), "\n")

  compatibility <- mapply(function(sample1, sample2) {

    # Get the genotypes for the two samples
    gt1 <- GT[, sample1]
    gt2 <- GT[, sample2]

    # Check if the genotypes are compatible
    compatible <- (sum(gt1 == gt2, na.rm = TRUE) / length(gt1)) * 100
    miss.perc <- (sum(is.na(gt1) | is.na(gt2))/ length(gt1)) * 100

    # Return the result
    return(c(compatible = compatible, miss.perc = miss.perc))
  }, filtered_combinations$Var1, filtered_combinations$Var2)


  result <- cbind(filtered_combinations, t(compatibility))
  colnames(result) <- c("sample1", "sample2", "%_matching_genotypes", "%_missing_genotypes")
  return(result)
}

