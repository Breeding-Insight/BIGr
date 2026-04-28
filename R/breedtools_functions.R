#' Compute allele frequencies for populations
#'
#' Computes allele frequencies for specified populations from SNP array data
#' coded as dosage of allele B.
#'
#' @param geno Numeric matrix of genotypes coded as dosage of allele B
#'   \code{{0, 1, 2, ..., ploidy}}, with individuals in rows (named) and
#'   SNPs in columns (named).
#' @param populations Named list of populations, each containing a character
#'   vector of individual IDs belonging to that population.
#' @param ploidy Integer. Ploidy level of the species. Default is `2`.
#'
#' @return A data frame of allele frequencies with SNPs as rows and populations
#'   as columns.
#'
#' @author Josué Chinchilla-Vargas
#'
#' @references Funkhouser SA, Bates RO, Ernst CW, Newcom D, Steibel JP.
#'   Estimation of genome-wide and locus-specific breed composition in pigs.
#'   _Transl Anim Sci._ 2017;1(1):36–44.
#'
#' @examples
#' geno_matrix <- matrix(
#'   c(4, 1, 4, 0,
#'     2, 2, 1, 3,
#'     0, 4, 0, 4,
#'     3, 3, 2, 2,
#'     1, 4, 2, 3),
#'   nrow = 4, ncol = 5, byrow = FALSE,
#'   dimnames = list(paste0("Ind", 1:4), paste0("S", 1:5))
#' )
#' pop_list <- list(
#'   PopA = c("Ind1", "Ind2"),
#'   PopB = c("Ind3", "Ind4")
#' )
#' allele_freqs <- allele_freq_poly(geno = geno_matrix, populations = pop_list, ploidy = 4)
#' print(allele_freqs)
#'
#' @export
allele_freq_poly <- function(geno, populations, ploidy = 2) {
  # Initialize result matrix
  X <- matrix(NA, nrow = ncol(geno), ncol = length(populations))
  for (i in 1:length(populations)) {
    pop_name <- names(populations[i])
    pop_geno <- geno[base::rownames(geno) %in% populations[[i]], ]
    al_freq  <- base::colMeans(pop_geno, na.rm = TRUE) / ploidy
    X[, i]   <- al_freq
  }
  base::colnames(X) <- base::names(populations)
  base::rownames(X) <- base::colnames(geno)
  return(X)
}


#' Solve breed composition for a single animal via quadratic programming
#'
#' Internal helper that solves the constrained OLS problem for one animal,
#' returning breed proportion estimates and the R² of the fit.
#'
#' @param Y Numeric vector of genotypes (named by SNP) for a single animal,
#'   coded as dosage of allele B \code{{0, 1, 2, ..., ploidy}}.
#' @param X Numeric matrix of allele frequencies with SNPs as rows and
#'   breeds/populations as columns.
#'
#' @return A named numeric vector of breed proportions plus `R2`.
#'
#' @references Funkhouser SA, Bates RO, Ernst CW, Newcom D, Steibel JP.
#'   Estimation of genome-wide and locus-specific breed composition in pigs.
#'   _Transl Anim Sci._ 2017;1(1):36–44.
#'
#' @import quadprog
#' @importFrom stats cor
#'
#' @noRd
QPsolve <- function(Y, X) {
  Ymod  <- Y[!base::is.na(Y)]
  Xmod  <- X[base::names(Ymod), ]
  p     <- base::ncol(X)
  nms   <- base::colnames(X)
  Rinv  <- base::solve(base::chol(base::t(Xmod) %*% Xmod))
  C     <- base::cbind(base::rep(1, p), base::diag(p))
  b2    <- c(1, base::rep(0, p))
  dd    <- base::t(Ymod) %*% Xmod
  qp    <- quadprog::solve.QP(Dmat = Rinv, factorized = TRUE, dvec = dd,
                              Amat = C, bvec = b2, meq = 1)
  beta  <- qp$solution
  rr    <- stats::cor(Ymod, Xmod %*% beta)^2
  result <- c(beta, rr)
  base::names(result) <- c(nms, "R2")
  return(result)
}


#' Compute genome-wide breed composition
#'
#' Estimates genome-wide breed/ancestry composition for a batch of animals
#' using quadratic programming, with optional pedigree-assisted and
#' grouped-output modes.
#'
#' @param Y Numeric matrix of genotypes with individuals as rows and SNPs as
#'   columns, coded as dosage of allele B \code{{0, 1, 2, ..., ploidy}}.
#' @param X Numeric matrix of allele frequencies with SNPs as rows and
#'   breeds/populations as columns.
#' @param ped Optional data frame with pedigree information formatted with
#'   columns `ID`, `Sire`, and `Dam`. When supplied, `QPsolve_par` is used
#'   and only animals with genotyped parents are included.
#' @param groups Optional named list of IDs grouped by breed/population.
#'   When supplied, results are returned as a named list partitioned by group.
#' @param mia Logical. Only applies when `ped` is supplied. If `TRUE`, returns
#'   the inferred maternally inherited allele per locus per animal. Default `FALSE`.
#' @param sire Logical. Only applies when `ped` is supplied. If `TRUE`, returns
#'   sire genotypes per locus per animal. Default `FALSE`.
#' @param dam Logical. Only applies when `ped` is supplied. If `TRUE`, returns
#'   dam genotypes per locus per animal. Default `FALSE`.
#' @param ploidy Integer. Ploidy level of the species. Default is `2`.
#'
#' @return A data frame, or a named list of data frames when `groups` is
#'   supplied, containing breed/ancestry composition estimates.
#'
#' @author Josué Chinchilla-Vargas
#'
#' @references Funkhouser SA, Bates RO, Ernst CW, Newcom D, Steibel JP.
#'   Estimation of genome-wide and locus-specific breed composition in pigs.
#'   _Transl Anim Sci._ 2017;1(1):36–44.
#'
#' @examples
#' allele_freqs_matrix <- matrix(
#'   c(0.625, 0.500,
#'     0.500, 0.500,
#'     0.500, 0.500,
#'     0.750, 0.500,
#'     0.625, 0.625),
#'   nrow = 5, ncol = 2, byrow = TRUE,
#'   dimnames = list(paste0("SNP", 1:5), c("VarA", "VarB"))
#' )
#' val_geno_matrix <- matrix(
#'   c(2, 1, 2, 3, 4,
#'     3, 4, 2, 3, 0),
#'   nrow = 2, ncol = 5, byrow = TRUE,
#'   dimnames = list(paste0("Test", 1:2), paste0("SNP", 1:5))
#' )
#' composition <- solve_composition_poly(Y = val_geno_matrix,
#'                                       X = allele_freqs_matrix,
#'                                       ploidy = 4)
#' print(composition)
#'
#' @import quadprog
#'
#' @export
solve_composition_poly <- function(Y,
                                   X,
                                   ped    = NULL,
                                   groups = NULL,
                                   mia    = FALSE,
                                   sire   = FALSE,
                                   dam    = FALSE,
                                   ploidy = 2) {
  # Transpose so Y is SNPs x animals, as required internally
  Y <- base::t(Y)
  # Retain only SNPs present in X
  Y <- Y[base::rownames(Y) %in% base::rownames(X), ]

  if (!base::is.null(ped)) {
    # Pedigree-assisted mode: use QPsolve_par for animals with genotyped parents
    mat_results     <- base::lapply(base::colnames(Y), QPsolve_par, Y, X, ped,
                                    mia = mia, sire = sire, dam = dam)
    mat_results_tab <- base::do.call(rbind, mat_results)
    return(mat_results_tab)

  } else if (!base::is.null(groups)) {
    # Grouped mode: standard computation, results partitioned by group
    Y               <- Y / ploidy
    grouped_results <- base::lapply(groups, QPseparate, Y, X)
    return(grouped_results)

  } else {
    # Standard mode: unsegregated computation across all animals
    Y       <- Y / ploidy
    results <- base::t(base::apply(Y, 2, QPsolve, X))
    return(results)
  }
}
