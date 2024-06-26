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
#' @param output.file output file name
#' @return A gzipped vcf file
#' @importFrom vcfR read.vcfR
#' @importFrom vcfR write.vcf
#' @importFrom vcfR maf
#' @importFrom vcfR extract.gt
#' @examples
#' ## Use file paths for each file on the local system
#'
#'
#' #filterVCF(vcf.file = "example_dart_Dosage_Report.csv",
#'  #          filter.OD = 0.5,
#'  #          ploidy = 2,
#'  #          output.file = "name_for_vcf")
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
                       output.file) {

  #Should allow for any INFO field to be entered to be filtered

  # Import VCF (can be .vcf or .vcf.gz)
  if (class(vcf.file) != "vcfR"){
    vcf <- read.vcfR(vcf.file)
  }

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
    cat("Filtering by DP\n")
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
    cat("Filtering by MPP\n")
    # Extract the MPP values
    mpp <- extract.gt(vcf, element = "MPP", as.numeric = TRUE)
    # Identify cells to modify based on the DP threshold
    threshold <- as.numeric(filter.MPP) #Need to make a variable for user to enter
    to_modify <- mpp > threshold
    # Replace cells in the vcf@gt matrix with NA format string where to_modify is TRUE
    vcf@gt[, -1][to_modify] <- na_format
    #remove extra matrices
    rm(to_modify)
    rm(mpp)
  }

  ## Filter based on INFO column (example: DP > 10)
  #Get INFO column
  info <- vcf@fix[, "INFO"]
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

  ##Filter by the INFO fields if present (maybe make small internal function)
  #DP
  #if ("DP" %in% info_ids) {
  #  dp_values <- as.numeric(sub(".*DP=([0-9]+);.*", "\\1", info))
  #  vcf <- vcf[dp_values > 1000, ]
  #}

  #OD
  if ("OD" %in% info_ids && !is.null(filter.OD)) {
    cat("Filtering by OD\n")
    od_values <- as.numeric(sub(".*OD=([0-9]+);.*", "\\1", info))
    vcf <- vcf[od_values < as.numeric(filter.OD), ]
    snps_removed <- starting_snps - nrow(vcf)
  }

  #BIAS (need to add user variables)
  if ("BIAS" %in% info_ids && !is.null(filter.BIAS.min) && !is.null(filter.BIAS.max)) {
    cat("Filtering by BIAS\n")
    bias_values <- as.numeric(sub(".*BIAS=([0-9]+);.*", "\\1", info))
    vcf <- vcf[bias_values > as.numeric(filter.BIAS.min) & bias_values < as.numeric(filter.BIAS.max), ]
  }

  #PMC (need to add user variables)
  if ("PMC" %in% info_ids && !is.null(filter.PMC)) {
    cat("Filtering by PMC\n")
    pmc_values <- as.numeric(sub(".*PMC=([0-9]+);.*", "\\1", info))
    vcf <- vcf[pmc_values < as.numeric(filter.PMC), ]
  }

  # Example: Filter based on missing data for samples and SNPs
  if (!is.null(filter.SAMPLE.miss) || !is.null(filter.SNP.miss)){
    gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)#as.matrix(vcfR2genlight(vcf))

    if (!is.null(filter.SNP.miss)) {
      cat("Filtering by SNP missing data\n")
      snp_missing_data <- rowMeans(is.na(gt_matrix))
      vcf <- vcf[snp_missing_data < as.numeric(filter.SNP.miss), ]
      gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
    }

    if (!is.null(filter.SAMPLE.miss)) {
      cat("Filtering by Sample missing data\n")
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

  ##Convert GT to dosage
  #gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)#as.matrix(vcfR2genlight(vcf))

  # Function to determine the ploidy level from a genotype string
  #determine_ploidy <- function(gt) {
  #  if (is.na(gt)) {
  #    return(NA)
  #  }
  #  return(length(strsplit(gt, "[|/]")[[1]]))
  #}

  # Function to find a non-NA example genotype to determine ploidy
  #find_example_gt <- function(matrix) {
  #  for (i in seq_len(nrow(matrix))) {
  #    for (j in seq_len(ncol(matrix))) {
  #      if (!is.na(matrix[i, j])) {
  #        return(matrix[i, j])
  #      }
  #    }
  #  }
  #  return(NA)  # Return NA if no non-NA genotype is found
  #}

  # Find a non-NA example genotype
  #example_gt <- find_example_gt(gt_matrix)

  # Determine the ploidy level
  #if (!is.na(example_gt)) {
  #  ploidy <- determine_ploidy(example_gt)
  #} else {
  #  stop("No non-NA genotype found to determine ploidy.")
  #}

  # Generate lookup table for genotypes to dosage conversion
  #generate_lookup_table <- function(ploidy) {
  #  possible_alleles <- 0:ploidy
  #  genotypes <- expand.grid(rep(list(possible_alleles), ploidy))
  #  genotypes <- apply(genotypes, 1, function(x) paste(x, collapse = "/"))
  #  dosage_values <- rowSums(expand.grid(rep(list(possible_alleles), ploidy)))
  #  lookup_table <- setNames(dosage_values, genotypes)
  #  return(lookup_table)
  #}

  # Generate the lookup table
  #lookup_table <- generate_lookup_table(ploidy)

  # Function to convert genotype to dosage using the lookup table
  #genotype_to_dosage <- function(gt, lookup_table) {
  #  if (is.na(gt)) {
  #    return(NA)
  #  }
  #  return(lookup_table[[gt]])
  #}

  # Function to convert genotype matrix to dosage matrix using vectorized operations
  #convert_genotypes_to_dosage <- function(gt_matrix, lookup_table) {
  #  unique_gts <- unique(gt_matrix)
  #  gt_to_dosage <- setNames(rep(NA, length(unique_gts)), unique_gts)
  #  valid_gts <- unique_gts[unique_gts %in% names(lookup_table)]
  #  gt_to_dosage[valid_gts] <- lookup_table[valid_gts]
  #  dosage_matrix <- gt_to_dosage[gt_matrix]
  #colnames(dosage_matrix) <- colnames(gt_matrix)
  #row.names(dosage_matrix) <- row.names(gt_matrix)
  #  return(matrix(as.numeric(dosage_matrix), nrow = nrow(gt_matrix), ncol = ncol(gt_matrix)))
  #}

  # Convert the genotype matrix to dosage matrix
  #dosage_matrix <- convert_genotypes_to_dosage(gt_matrix, lookup_table)

  ##MAF filter
  #Compare my lengthy process to estimate MAF with vcfR::maf() function
  #The BIGr::calculate_MAF(dosage_matrix, ploidy) is the exact same as the vcfR::maf() calculations
  #The step where I extract UD and calculate MAF is different...
  #if ("UD" %in% format_fields) {
  #  maf_df <- BIGr::calculate_MAF(extract.gt(vcf, element = "UD", as.numeric = TRUE), ploidy = ploidy)
  #} else {
  #convert genotypes to dosage and filter
  #  maf_df <- BIGr::calculate_MAF(dosage_matrix, ploidy)
  #}
  #Need to confirm that vcfR::maf will work with any ploidy...if not, use my code
  if (!is.null(filter.MAF)) {
    cat("Filtering by MAF\n")
    maf_df <- data.frame(vcfR::maf(vcf, element = 2))
    vcf <- vcf[maf_df$Frequency > as.numeric(filter.MAF), ]
  }
  ### Export the modified VCF file (this exports as a .vcf.gz, so make sure to have the name end in .vcf.gz)
  cat("Exporting VCF\n")
  if (!class(vcf.file) == "vcfR"){
    base_name <- basename(vcf.file)
    output_name <- paste0("filtered_",base_name,".vcf.gz")
    vcfR::write.vcf(vcf, file = output_name)
  }else{
    base_name <- output.file
    output_name <- paste0("filtered_",base_name,".vcf.gz")
    vcfR::write.vcf(vcf, file = output_name)
  }

  #Message that includes the output vcf stats
  print(vcf)

  #Message
  samples_removed <- starting_samples - (ncol(vcf@gt)-1)
  SNPs_removed <- starting_snps - nrow(vcf)
  message("Samples removed due to filtering: ",samples_removed)
  message("SNPs removed due to filtering: ",SNPs_removed)
  message("Complete!")
}
#This is not reliable, so no longer use this shortcut to get dosage matrix
#test2 <- vcfR2genlight(vcf)


#####Testing custom VCF reading function######
# Open the gzipped VCF file
#con <- gzfile("/Users/ams866/Desktop/output.vcf", "rt")

# Read in the entire file
#lines <- readLines(con)
#close(con)
# Read in the entire file
#lines <- readLines("/Users/ams866/Desktop/output.vcf")
# Filter out lines that start with ##
#filtered_lines <- lines[!grepl("^##", lines)]
# Create a temporary file to write the filtered lines
#temp_file <- tempfile()
#writeLines(filtered_lines, temp_file)
# Read in the filtered data using read.table or read.csv
#vcf_data <- read.table(temp_file, header = TRUE, sep = "\t", comment.char = "", check.names = FALSE)
# Clean up the temporary file
#unlink(temp_file)

##Extract INFO column and Filter SNPs by those values
#Update the filtering options by the items present in the INFO column?

# Load required library
#library(dplyr)

# Split INFO column into key-value pairs
#vcf_data_parsed <- vcf_data %>%
#  mutate(INFO_PARSED = strsplit(INFO, ";")) %>%
#  unnest(INFO_PARSED) %>%
#  separate(INFO_PARSED, into = c("KEY", "VALUE"), sep = "=") %>%
#  spread(KEY, VALUE)

#Filter by DP
#filtered_vcf_data <- vcf_data_parsed %>%
#  filter(as.numeric(DP) > 10)

# View the filtered dataframe
#print(filtered_vcf_data)

##Extracting and filtering by FORMAT column
# Identify the columns that are not sample columns
#non_sample_cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
# Identify the sample columns
#sample_cols <- setdiff(names(vcf_data), non_sample_cols)
# Extract FORMAT keys
#format_keys <- strsplit(as.character(vcf_data$FORMAT[1]), ":")[[1]]
# Split SAMPLE columns based on FORMAT
#vcf_data_samples <- vcf_data %>%
#  mutate(across(all_of(sample_cols), ~strsplit(as.character(.), ":"))) %>%
#  mutate(across(all_of(sample_cols), ~map(., ~setNames(as.list(.), format_keys)))) %>%
#  unnest_wider(all_of(sample_cols), names_sep = "_")

# View the parsed dataframe
#print(head(vcf_data_samples))

# Create separate dataframes for each FORMAT variable
#format_dfs <- lapply(format_keys, function(format_key) {
#  vcf_data_samples %>%
#    select(ID, ends_with(paste0("_", format_key))) %>%
#    column_to_rownames("ID")
#})

# Assign names to the list elements
#names(format_dfs) <- format_keys

# Access the separate dataframes
#gt_df <- format_dfs$GT  # Genotype dataframe
#ad_df <- format_dfs$AD  # Allelic depths dataframe

#*I think the above method is okay if you only need to filter at the INFO level,
#*But I think if you want to filter for FORMAT, that vcfR is probably best,
#*Will need to explore further if I can easily just filter for MPP by checking if it is above a
#*threshold, and then converting the GT and UD values to NA if so...
#*If that is efficient and works, then I will just use this custom VCF method...
