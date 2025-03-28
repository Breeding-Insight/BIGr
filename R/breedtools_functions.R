#' Computes allele frequencies for specified populations given SNP array data
#'
#' @param geno matrix of genotypes coded as the dosage of allele B \code{{0, 1, 2, ..., ploidy}}
#'  with individuals in rows (named) and SNPs in columns (named)
#' @param populations list of named populations. Each population has a vector of IDs
#'  that belong to the population. Allele frequencies will be derived from all animals
#' @param ploidy integer indicating the ploidy level (default is 2 for diploid)
#' @return data.frame consisting of allele_frequencies for populations (columns) for
#'  each SNP (rows)
#' @references Funkhouser SA, Bates RO, Ernst CW, Newcom D, Steibel JP. Estimation of genome-wide and locus-specific
#' breed composition in pigs. Transl Anim Sci. 2017 Feb 1;1(1):36-44.
#'
#' @examples
#' # Example inputs
#' geno_matrix <- matrix(
#' c(4, 1, 4, 0, # S1
#'   2, 2, 1, 3, # S2
#'   0, 4, 0, 4, # S3
#'   3, 3, 2, 2, # S4
#'   1, 4, 2, 3),# S5
#' nrow = 4, ncol = 5, byrow = FALSE, # individuals=rows, SNPs=cols
#' dimnames = list(paste0("Ind", 1:4), paste0("S", 1:5))
#' )
#'
#'pop_list <- list(
#' PopA = c("Ind1", "Ind2"),
#' PopB = c("Ind3", "Ind4")
#' )
#'
#' allele_freqs <- allele_freq_poly(geno = geno_matrix, populations = pop_list, ploidy = 4)
#' print(allele_freqs)
#'
#' @export
allele_freq_poly <- function(geno, populations, ploidy = 2) {

  # Initialize returned df
  X <- matrix(NA, nrow = ncol(geno), ncol = length(populations))

  # Subset geno into different populations
  for (i in 1:length(populations)) {

    # Get name of ith item in the list (population name)
    pop_name <- names(populations[i])

    # Subset geno to only include genotypes of IDs in pop
    pop_geno <- geno[rownames(geno) %in% populations[[i]], ]

    # Calculate allele frequencies
    al_freq <- colMeans(pop_geno, na.rm = TRUE) / ploidy

    # Add to X
    X[, i] <- al_freq
  }

  # Label X with populations and SNPs
  colnames(X) <- names(populations)
  rownames(X) <- colnames(geno)

  return(X)
}


#' Performs whole genome breed composition prediction.
#'
#' @param Y numeric vector of genotypes (with names as SNPs) from a single animal.
#'   coded as dosage of allele B \code{{0, 1, 2, ..., ploidy}}
#' @param X numeric matrix of allele frequencies from reference animals
#' @param p numeric indicating number of breeds represented in X
#' @param names character names of breeds
#' @return data.frame of breed composition estimates
#' @import quadprog
#' @importFrom stats cor
#' @references Funkhouser SA, Bates RO, Ernst CW, Newcom D, Steibel JP. Estimation of genome-wide and locus-specific
#' breed composition in pigs. Transl Anim Sci. 2017 Feb 1;1(1):36-44.
#'
#' @noRd
QPsolve <- function(Y, X) {

  # Remove NAs from Y and remove corresponding
  #   SNPs from X. Ensure Y is numeric
  Ymod <- Y[!is.na(Y)]
  Xmod <- X[names(Ymod), ]

  # Determine properties from X matrix - the number of parameters (breeds) p
  #   and the names of those parameters.
  p <- ncol(X)
  names <- colnames(X)

  # perfom steps needed to solve OLS by framing
  # as a QP problem
  # Rinv - matrix should be of dimensions px(p+1) where p is the number of variables in X
  Rinv <- solve(chol(t(Xmod) %*% Xmod))

  # C - the first column is a sum restriction (all equal to 1) and the rest of the columns an identity matrix
  C <- cbind(rep(1, p), diag(p))

  # b2 - This should be a vector of length p+1 the first element is the value of the sum (1)
  #   the other elements are the restriction of individual coefficients (>)
  #   so a value 0 produces positive coefficients
  b2 <- c(1, rep(0, p))

  # dd - this should be a matrix NOT a vector
  dd <- (t(Ymod) %*% Xmod)

  qp <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = dd, Amat = C, bvec = b2, meq = 1)
  beta <- qp$solution
  rr <- cor(Ymod, Xmod %*% beta)^2
  result <- c(beta, rr)
  names(result) <- c(names, "R2")
  return(result)
}


#' Compute genome-wide breed composition
#'
#' Computes genome-wide breed/ancestry composition using quadratic programming on a
#' batch of animals.
#'
#' @param Y numeric matrix of genotypes (columns) from all animals (rows) in population
#'  coded as dosage of allele B \code{{0, 1, 2, ..., ploidy}}
#' @param X numeric matrix of allele frequencies (rows) from each reference panel (columns). Frequencies are
#'  relative to allele B.
#' @param ped data.frame giving pedigree information. Must be formatted "ID", "Sire", "Dam"
#' @param groups list of IDs categorized by breed/population. If specified, output will be a list
#'  of results categorized by breed/population.
#' @param mia logical. Only applies if ped argument is supplied. If true, returns a data.frame
#'  containing the inferred maternally inherited allele for each locus for each animal instead
#'  of breed composition results.
#' @param sire logical. Only applies if ped argument is supplied. If true, returns a data.frame
#'  containing sire genotypes for each locus for each animal instead of breed composition results.
#' @param dam logical. Only applies if ped argument is supplied. If true, returns a data.frame
#'  containing dam genotypes for each locus for each animal instead of breed composition results.
#' @param ploidy integer. The ploidy level of the species (e.g., 2 for diploid, 3 for triploid, etc.).
#' @return A data.frame or list of data.frames (if groups is !NULL) with breed/ancestry composition
#'  results
#' @import quadprog
#' @references Funkhouser SA, Bates RO, Ernst CW, Newcom D, Steibel JP. Estimation of genome-wide and locus-specific
#' breed composition in pigs. Transl Anim Sci. 2017 Feb 1;1(1):36-44.
#'
#' @examples
#' # Example inputs for solve_composition_poly (ploidy = 4)
#'
#' # (This would typically be the output from allele_freq_poly)
#' allele_freqs_matrix <- matrix(
#'   c(0.625, 0.500,
#'     0.500, 0.500,
#'     0.500, 0.500,
#'     0.750, 0.500,
#'     0.625, 0.625),
#'   nrow = 5, ncol = 2, byrow = TRUE,
#'   dimnames = list(paste0("SNP", 1:5), c("VarA", "VarB"))
#' )
#'
#' # Validation Genotypes (individuals x SNPs)
#' val_geno_matrix <- matrix(
#'   c(2, 1, 2, 3, 4,  # Test1 dosages for SNP1-5
#'     3, 4, 2, 3, 0), # Test2 dosages for SNP1-5
#'   nrow = 2, ncol = 5, byrow = TRUE,
#'   dimnames = list(paste0("Test", 1:2), paste0("SNP", 1:5))
#' )
#'
#' # Calculate Breed Composition
#' composition <- solve_composition_poly(Y = val_geno_matrix,
#'                                       X = allele_freqs_matrix,
#'                                       ploidy = 4)
#' print(composition)
#'
#' @export
solve_composition_poly <- function(Y,
                                   X,
                                   ped = NULL,
                                   groups = NULL,
                                   mia = FALSE,
                                   sire = FALSE,
                                   dam = FALSE,
                                   ploidy = 2) {

  # Functions require Y to be animals x SNPs. Transpose
  Y <- t(Y)

  # SNPs in Y should only be the ones present in X
  Y <- Y[rownames(Y) %in% rownames(X), ]

  # If ped is supplied, use QPsolve_par to compute genomic composition using
  #   only animals who have genotyped parents (by incorporating Sire genotype).
  if (!is.null(ped)) {
    mat_results <- lapply(colnames(Y),
                          QPsolve_par,
                          Y,
                          X,
                          ped,
                          mia = mia,
                          sire = sire,
                          dam = dam)

    mat_results_tab <- do.call(rbind, mat_results)
    return (mat_results_tab)

    # Else if groups supplied - perform regular genomic computation
    #   and list results by groups
  } else if (!is.null(groups)) {

    # When using regular genomic computation - adjust dosage based on ploidy
    Y <- Y / ploidy #(default is 2)

    grouped_results <- lapply(groups, QPseparate, Y, X)
    return (grouped_results)

    # If neither using the ped or grouping option - just perform normal, unsegregated
    #   calculation
  } else {

    # Adjust dosage based on ploidy (default is 2)
    Y <- Y / ploidy

    results <- t(apply(Y, 2, QPsolve, X))
    return (results)
  }
}
