#' Estimate Minimum Number of Individuals to Sample to Capture Population Genomic Diversity (Genotype Matrix)
#'
#' This function can be used to estimate the number of individuals to sample from a population
#' in order to capture a desired percentage of the genomic diversity.
#' It assumes that the samples are the columns, and the genomic markers are in rows. Missing data should
#' be set as NA, which will then be ignored for the calculations. All samples must have the same ploidy.
#' This function was adapted from a previously developed Python method (Sandercock et al., 2023)
#' (https://github.com/alex-sandercock/Capturing_genomic_diversity/)
#'
#' @param df Genotype matrix or data.frame with the count of alternate alleles (0=homozygous reference, 1 = heterozygous, 2 = homozygous alternate)
#' @param ploidy The ploidy of the species being analyzed
#' @param r2_threshold The ratio of diversity to capture (default = 0.9)
#' @param iterations The number of iterations to perform to estimate the average result (default = 10)
#' @param sample_list The list of samples to subset from the dataset (optional)
#' @param parallel Run the analysis in parallel (True/False) (default = FALSE)
#' @param save.result Save the results to a .txt file? (default = TRUE)
#' @return A data.frame with minimum number of samples required to match or exceed the input ratio
#' @import foreach
#' @import doParallel
#' @import dplyr
#' @importFrom Rdpack reprompt
#' @importFrom stats lm qt sd
#' @importFrom utils write.table
#' @references
#' Sandercock, A. M., Westbrook, J. W., Zhang, Q., & Holliday, J. A. (2024). The road to restoration: Identifying and conserving the adaptive legacy of American chestnut. PNAS (in press).
#' @export
capture_diversity.Gmat <- function(df, ploidy, r2_threshold=0.9, iterations = 10, sample_list = NULL, parallel=FALSE, save.result=TRUE) {
##Need to make sure these two packages are loaded with BIGr (vcfR and dplyr,"foreach","doParallel"

  #Add batch
  batch=1

  # This  will subset it based on the user-supplied list
  if (!is.null(sample_list)) {
    df <- df[, colnames(df) %in% sample_list]
  }

  ############ 2. Perform bootstrap sampling ###########################

  bootstrap_batch_sample_regression <- function(df, batch=1, target_values, r2_threshold) {
    # Perform bootstrap sampling until threshold of diversity captured is met based on R2 value
    df_merged <- data.frame(matrix(ncol = 0, nrow = nrow(df)))
    rownames(df_merged) <- rownames(df)
    sampling_round <- 0
    diversity_captured <- FALSE

    while (!diversity_captured) {
      sampling_round <- sampling_round + 1
      sampled <- df[, sample(colnames(df), as.numeric(batch))]
      df_merged <- cbind(df_merged, sampled)
      af_df <- calculate_MAF(df_merged, ploidy)
      af_values <- af_df$AF
      lm_model <- stats::lm(target_values ~ af_values)
      r2 <- summary(lm_model)$r.squared

      if (r2 >= r2_threshold) {
        diversity_captured <- TRUE
      }
    }

    return(as.numeric(sampling_round))
  }

  ############### 3. Calculate statistics ################################

  calculate_statistics <- function(sampling_round_list, batch) {
    # Convert sampling rounds to number of individuals
    sampling_round_list <- lapply(sampling_round_list, function(x) x * as.numeric(batch))
    mean_ind <- mean(unlist(sampling_round_list))
    stand_dev <- stats::sd(unlist(sampling_round_list))
    CI <- stats::qt(0.975, df = length(sampling_round_list) - 1) * (stand_dev / sqrt(length(sampling_round_list)))
    CI_lower <- mean_ind - CI
    CI_upper <- mean_ind + CI
    #Add results to dataframe
    results_df <- data.frame(
      Individuals = mean_ind,
      CI_Lower = CI_lower,
      CI_Upper = CI_upper,
      Iterations = length(sampling_round_list)
    )
    return(results_df)
  }

  target_AF_values <- calculate_MAF(df, ploidy)$AF
  sample_round_list <- list()

  #Perform iterations depending on user parallel selection
  if (!parallel) {
    for (i in 1:iterations) {
      sample_round_list[i] <- bootstrap_batch_sample_regression(df, batch, target_AF_values, r2_threshold)
    }
  } else {
    ###Parallel script
    # Detect the number of available cores and subtract 1
    #Need to make sure num_cores is 1 if that is all that is available..
    num_cores <- parallel::detectCores() - 1
    if (num_cores == 0) {
      num_cores = 1
    }

    # Set up parallel backend
    cl <- parallel::makeCluster(num_cores)
    registerDoParallel(cl)
    # Perform boostrap sampling over user input or default iterations
    sample_round_list <- foreach(i = 1:iterations, .combine = c) %dopar% {
      bootstrap_batch_sample_regression(df, batch, target_AF_values, r2_threshold)
    }
    # Stop the parallel backend
    parallel::stopCluster(cl)
  }

  final_df <- calculate_statistics(sample_round_list, batch)

  #Save results to a .txt file
  if (save.result){
    utils::write.table(final_df, file= "capture_diversity_output.txt", row.names=FALSE)
  }else{
    cat("Number of individuals to sample =", final_df$Individuals, "\n95% Confidence Intervals =", final_df$CI_Lower, "-", final_df$CI_Upper, "\nIterations performed =", final_df$Iterations, "\n")
  }

  #Return the results dataframe
  return(final_df)
}


