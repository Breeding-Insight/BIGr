#' Filter a VCF file
#'
#' This function will filter a VCF file or vcfR object and export the updated version
#'
#' This function will input a VCF file or vcfR object and filter based on the user defined options.
#' The output file will be saved to the location and with the name that is specified.
#' The VCF format is v4.3
#'
#' @param vcf.file vcfR object or path to VCF file. Can be unzipped (.vcf) or gzipped (.vcf.gz).
#' @param filter.OD Updog filter
#' @param filter.BIAS.min Updog filter (requires a value for both BIAS.min and BIAS.max)
#' @param filter.BIAS.max Updog filter (requires a value for both BIAS.min and BIAS.max)
#' @param filter.DP Total read depth at each SNP filter
#' @param filter.MPP Updog filter
#' @param filter.PMC Updog filter
#' @param filter.MAF Minor allele frequency filter
#' @param filter.SAMPLE.miss Sample missing data filter
#' @param filter.SNP.miss SNP missing data filter
#' @param ploidy The ploidy of the species being analyzed
#' @param output.file output file name (optional). If no output.file name provided, then a vcfR object will be returned.
#' @return A gzipped vcf file
#' @importFrom vcfR read.vcfR
#' @importFrom vcfR write.vcf
#' @importFrom vcfR maf
#' @importFrom vcfR extract.gt
#' @examples
#' ## Use file paths for each file on the local system
#'
#' #Temp location (only for example)
#' output_file <- tempfile()
#'
#' filterVCF(vcf.file = system.file("iris_DArT_VCF.vcf.gz", package = "BIGr"),
#'            filter.OD = 0.5,
#'            filter.MAF = 0.05,
#'            ploidy = 2,
#'            output.file = output_file)
#'
#' # Removing the output for the example
#' rm(output_file)
#'
#' ##The function will output the filtered VCF to the current working directory
#'
#' @export
filterVCF <- function(vcf.file,
                       filter.OD = NULL,
                       filter.BIAS.min = NULL,
                       filter.BIAS.max = NULL,
                       filter.DP = NULL,
                       filter.MPP = NULL,
                       filter.PMC = NULL,
                       filter.MAF = NULL,
                       filter.SAMPLE.miss = NULL,
                       filter.SNP.miss = NULL,
                       ploidy,
                       output.file = NULL) {

  #Should allow for any INFO field to be entered to be filtered

  # Import VCF (can be .vcf or .vcf.gz)
  if (!inherits(vcf.file, "vcfR")) {
    vcf <- read.vcfR(vcf.file, verbose = FALSE)
  } else {
    vcf <- vcf.file
    #rm(vcf.file)
  }

  #Update header based on user filtering parameters
  param_list <- list(
    filter.OD = filter.OD,
    filter.BIAS.min = filter.BIAS.min,
    filter.BIAS.max = filter.BIAS.max,
    filter.DP = filter.DP,
    filter.MPP = filter.MPP,
    filter.PMC = filter.PMC,
    filter.MAF = filter.MAF,
    filter.SAMPLE.miss = filter.SAMPLE.miss,
    filter.SNP.miss = filter.SNP.miss,
    ploidy = ploidy
  )

  # Filter out NULL values
  param_list <- param_list[!sapply(param_list, is.null)]

  #Update header lines and append
  header_line <- paste0('##BIGr_filterVCFparameters, ', paste(names(param_list), unlist(param_list), sep="=", collapse=", "),"; ",Sys.time())
  vcf@meta <- c(vcf@meta, paste0('##BIGr_filterVCF=', packageVersion("BIGr")), header_line)

  #Getting starting number of SNPs and Samples
  starting_snps <- nrow(vcf)
  starting_samples <- ncol(vcf@gt)-1 #subtract 1 to not include the FORMAT column

  # Determine the number of items in the FORMAT field
  format_string <- vcf@gt[1, "FORMAT"]
  format_fields <- strsplit(format_string, ":")[[1]]
  num_fields <- length(format_fields)
  gt_pos <- which(format_fields == "GT")
  # Create the NA format string, replacing only GT with "./."
  missing_gt <- paste(rep(".", ploidy), collapse = "/")
  na_fields <- rep(".", num_fields)
  na_fields[gt_pos] <- missing_gt
  na_format <- paste(na_fields, collapse = ":")

  # Extract the DP values
  if ("DP" %in% format_fields && !is.null(filter.DP)) {
    message("Filtering by DP\n")
    dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
    # Identify cells to modify based on the DP threshold
    threshold <- as.numeric(filter.DP)
    to_modify <- dp < threshold
    # Replace cells in the vcf@gt matrix with NA format string where to_modify is TRUE
    vcf@gt[, -1][to_modify] <- na_format
    # Remove extra matrices
    rm(to_modify)
    rm(dp)
  }

  #Filter if the MPP field is present
  if ("MPP" %in% format_fields && !is.null(filter.MPP)) {
    message("Filtering by MPP\n")
    # Extract the MPP values
    mpp <- extract.gt(vcf, element = "MPP", as.numeric = TRUE)
    # Identify cells to modify based on the DP threshold
    threshold <- as.numeric(filter.MPP) #Need to make a variable for user to enter
    to_modify <- mpp < threshold
    # Replace cells in the vcf@gt matrix with NA format string where to_modify is TRUE
    vcf@gt[, -1][to_modify] <- na_format
    #remove extra matrices
    rm(to_modify)
    rm(mpp)
  }

  ## Filter based on INFO column (example: DP > 10)

  # Get INFO column
  info <- vcf@fix[, "INFO"] #Need to get after each filter..

  # Function to extract a specific INFO field value
  extract_info_value <- function(info, field) {
    pattern <- paste0(".*", field, "=([0-9.]+).*")
    values <- as.numeric(sub(pattern, "\\1", info))
    return(values)
  }

  # Function to extract INFO IDs from a single INFO string
  extract_info_ids <- function(info_string) {
    # Split the INFO string by ';'
    info_parts <- strsplit(info_string, ";")[[1]]
    # Extract the part before the '=' in each segment
    info_ids <- gsub("=.*", "", info_parts)
    return(info_ids)
  }

  # Apply the function to the first INFO string
  info_ids <- extract_info_ids(info[1])

  # Filtering by OD
  if ("OD" %in% info_ids && !is.null(filter.OD)) {
    info <- vcf@fix[, "INFO"] #Need to get after each filter..
    message("Filtering by OD\n")
    od_values <- extract_info_value(info, "OD")
    # Ensure no NA values before filtering
    if (!all(is.na(od_values))) {
      vcf <- vcf[od_values < as.numeric(filter.OD), ]
    } else {
      warning("No valid OD values found.\n")
    }
  }

  info <- vcf@fix[, "INFO"] #Need to get after each filter..

  # Filtering by BIAS
  if ("BIAS" %in% info_ids && !is.null(filter.BIAS.min) && !is.null(filter.BIAS.max)) {
    info <- vcf@fix[, "INFO"] #Need to get after each filter..
    message("Filtering by BIAS\n")
    bias_values <- extract_info_value(info, "BIAS")
    # Ensure no NA values before filtering
    if (!all(is.na(bias_values))) {
      vcf <- vcf[bias_values > as.numeric(filter.BIAS.min) & bias_values < as.numeric(filter.BIAS.max), ]
    } else {
      warning("No valid BIAS values found.\n")
    }
  }

  # Filtering by PMC
  if ("PMC" %in% info_ids && !is.null(filter.PMC)) {
    info <- vcf@fix[, "INFO"] #Need to get after each filter..
    message("Filtering by PMC\n")
    pmc_values <- extract_info_value(info, "PMC")
    # Ensure no NA values before filtering
    if (!all(is.na(pmc_values))) {
      vcf <- vcf[pmc_values < as.numeric(filter.PMC), ]
    } else {
      warning("No valid PMC values found.\n")
    }
  }

  # Example: Filter based on missing data for samples and SNPs
  if (!is.null(filter.SAMPLE.miss) || !is.null(filter.SNP.miss)){
    info <- vcf@fix[, "INFO"] #Need to get after each filter..
    gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)#as.matrix(vcfR2genlight(vcf))

    if (!is.null(filter.SNP.miss)) {
      message("Filtering by SNP missing data\n")
      snp_missing_data <- rowMeans(is.na(gt_matrix))
      vcf <- vcf[snp_missing_data < as.numeric(filter.SNP.miss), ]
      gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
    }

    if (!is.null(filter.SAMPLE.miss)) {
      message("Filtering by Sample missing data\n")
      # Calculate the proportion of missing data for each sample
      sample_missing_data <- colMeans(is.na(gt_matrix))
      # Identify samples to keep based on the missing data threshold
      samples_to_keep <- names(sample_missing_data)[sample_missing_data < as.numeric(filter.SAMPLE.miss)]
      # Include "FORMAT" column in the samples to keep
      samples_to_keep <- c("FORMAT", samples_to_keep)
      # Subset the VCF object to keep only the desired samples
      vcf <- vcf[, colnames(vcf@gt) %in% samples_to_keep]
    }

    # Remove matrices
    rm(gt_matrix)
  }

  ##MAF filter
  if (!is.null(filter.MAF)) {
    message("Filtering by MAF\n")
    maf_df <- data.frame(vcfR::maf(vcf, element = 2))
    vcf <- vcf[maf_df$Frequency > as.numeric(filter.MAF), ]
  }
  ### Export the modified VCF file (this exports as a .vcf.gz, so make sure to have the name end in .vcf.gz)
  message("Exporting VCF\n")
  if (!inherits(vcf.file, "vcfR")) {
    if (!is.null(output.file)) {
      output_name <- paste0(output.file, ".vcf.gz")
      vcfR::write.vcf(vcf, file = output_name)
    } else {
      return(vcf)
    }
  } else {
    if (!is.null(output.file)) {
      output_name <- paste0(output.file, "_filtered.vcf.gz")
      vcfR::write.vcf(vcf, file = output_name)
    } else {
      return(vcf)
    }
  }

  #Message
  samples_removed <- starting_samples - (ncol(vcf@gt)-1)
  SNPs_removed <- starting_snps - nrow(vcf)
  message("Samples removed due to filtering: ",samples_removed)
  message("SNPs removed due to filtering: ",SNPs_removed)
  message("Complete!")
}
