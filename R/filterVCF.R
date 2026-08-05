#' Filter a VCF file
#'
#' This function will filter a VCF file or vcfR object and export the updated version
#'
#' This function will input a VCF file or vcfR object and filter based on the user defined options.
#' The output file will be saved to the location and with the name that is specified.
#' The VCF format is v4.3
#'
#' Filters are applied in a fixed order, which matters because each one acts on
#' what the previous ones left behind. `filter.DP` and `filter.MPP` run first and
#' set individual genotype calls to missing. `filter.OD`, `filter.BIAS.min` and
#' `filter.BIAS.max`, and `filter.PMC` then remove variants using the INFO column.
#' `filter.SNP.miss` removes variants and `filter.SAMPLE.miss` removes samples,
#' both counting the calls that were set to missing above. `filter.MAF` runs last,
#' so it is calculated from the samples that remain.
#'
#' A variant is removed if the value a requested filter needs cannot be read from
#' it, for instance because the INFO field is absent from that record or because
#' every genotype call for the variant is missing. The number of variants removed
#' for this reason is reported in a warning.
#'
#' @param vcf.file vcfR object or path to VCF file. Can be unzipped (.vcf) or gzipped (.vcf.gz).
#' @param filter.OD Maximum overdispersion, read from the `OD` field of the INFO
#'   column as estimated by updog. Variants with an `OD` below this value are kept.
#' @param filter.BIAS.min Minimum allele bias, read from the `BIAS` field of the
#'   INFO column as estimated by updog. Variants with a `BIAS` above this value are
#'   kept. Has no effect unless `filter.BIAS.max` is also supplied.
#' @param filter.BIAS.max Maximum allele bias, read from the `BIAS` field of the
#'   INFO column. Variants with a `BIAS` below this value are kept. Has no effect
#'   unless `filter.BIAS.min` is also supplied.
#' @param filter.DP Minimum read depth for a genotype call. Calls whose FORMAT/DP
#'   is below this value are set to missing. This does not remove variants, though
#'   the calls it sets to missing then count towards `filter.SNP.miss`,
#'   `filter.SAMPLE.miss` and `filter.MAF`.
#' @param filter.MPP Minimum posterior probability for a genotype call, as
#'   reported by updog. Calls whose FORMAT/MPP is below this value are set to
#'   missing, on the same terms as `filter.DP`.
#' @param filter.PMC Maximum proportion of individuals misclassified, read from
#'   the `PMC` field of the INFO column as estimated by updog. Variants with a
#'   `PMC` below this value are kept.
#' @param filter.MAF Minimum minor allele frequency. Variants with a minor allele
#'   frequency above this value are kept. Calculated after any sample removal, so
#'   it reflects only the samples that survive `filter.SAMPLE.miss`.
#' @param filter.SAMPLE.miss Maximum proportion of missing genotype calls a sample
#'   may have, between 0 and 1. Samples missing a smaller proportion than this are
#'   kept. Note that this is a proportion and not a percentage.
#' @param filter.SNP.miss Maximum proportion of missing genotype calls a variant
#'   may have, between 0 and 1. Variants missing a smaller proportion than this are
#'   kept. Note that this is a proportion and not a percentage.
#' @param ploidy The ploidy of the species being analyzed. Required. Used to write
#'   the missing genotype, so a diploid is recorded as `./.` and a tetraploid as
#'   `./././.` when a call is filtered out.
#' @param output.file Output file name, without an extension (optional). When
#'   supplied, the filtered VCF is written as a gzipped file and nothing is
#'   returned. `.vcf.gz` is appended when `vcf.file` is a path, and
#'   `_filtered.vcf.gz` when `vcf.file` is a vcfR object. When not supplied, a
#'   vcfR object is returned instead and no file is written.
#' @return A vcfR object when `output.file` is not supplied. Otherwise the
#'   filtered VCF is written to disk and nothing is returned.
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

  # A threshold that is not a single readable number turns the comparison it
  # feeds into an NA index, which silently corrupts or skips the filter rather
  # than erroring, so all thresholds are validated before anything else happens.
  validate_filter_number <- function(x, name) {
    if (is.null(x)) return(NULL)
    if (!(is.numeric(x) || is.character(x)) || length(x) != 1L) {
      stop(name, " must be a single numeric value or NULL.", call. = FALSE)
    }
    value <- suppressWarnings(as.numeric(x))
    if (!is.finite(value)) {
      stop(name, " must be a single finite numeric value or NULL, not ",
           deparse(x), ".", call. = FALSE)
    }
    value
  }

  filter.OD          <- validate_filter_number(filter.OD, "filter.OD")
  filter.BIAS.min    <- validate_filter_number(filter.BIAS.min, "filter.BIAS.min")
  filter.BIAS.max    <- validate_filter_number(filter.BIAS.max, "filter.BIAS.max")
  filter.DP          <- validate_filter_number(filter.DP, "filter.DP")
  filter.MPP         <- validate_filter_number(filter.MPP, "filter.MPP")
  filter.PMC         <- validate_filter_number(filter.PMC, "filter.PMC")
  filter.MAF         <- validate_filter_number(filter.MAF, "filter.MAF")
  filter.SAMPLE.miss <- validate_filter_number(filter.SAMPLE.miss, "filter.SAMPLE.miss")
  filter.SNP.miss    <- validate_filter_number(filter.SNP.miss, "filter.SNP.miss")

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
  # The field is located by splitting on ';' so that the name is matched against a
  # whole INFO entry, and the value is passed to as.numeric() unaltered so that any
  # valid numeric representation is accepted. This includes scientific notation
  # (e.g. PMC=4.78e-07), which updog2vcf() writes whenever a value is small enough
  # for R's default formatting to use it.
  extract_info_value <- function(info, field) {
    prefix <- paste0(field, "=")
    values <- vapply(strsplit(info, ";", fixed = TRUE), function(parts) {
      hit <- parts[startsWith(parts, prefix)]
      if (length(hit) == 0L) {
        NA_real_
      } else {
        suppressWarnings(as.numeric(substring(hit[1L], nchar(prefix) + 1L)))
      }
    }, numeric(1))
    return(values)
  }

  # Build a row selection from a filter comparison. An NA in a logical index does
  # not drop a variant, it inserts an all-NA row into the VCF, so variants whose
  # value could not be read are removed explicitly here. The count is reported so
  # that they are never discarded silently.
  select_variants <- function(keep, values, label, hint) {
    unusable <- is.na(values)
    if (length(values) > 0 && all(unusable)) {
      # Losing every variant is usually a problem with the input rather than a
      # genuine result, so this case says what to look at instead of only
      # reporting the count.
      warning("No readable ", label, " values were found, so all ", length(values),
              " variants were removed. ", hint, call. = FALSE)
    } else if (any(unusable)) {
      warning(sum(unusable), " of ", length(values), " variants had a missing or ",
              "unreadable ", label, " value and were removed.", call. = FALSE)
    }
    return(!is.na(keep) & keep & !unusable)
  }

  # A requested filter is applied to every record. Which INFO fields are present
  # is not decided from the first record, because a field may legally be absent
  # from any individual record; records without a usable value are removed by
  # select_variants() and reported there.

  # Filtering by OD
  if (!is.null(filter.OD)) {
    info <- vcf@fix[, "INFO"] #Need to get after each filter..
    message("Filtering by OD\n")
    od_values <- extract_info_value(info, "OD")
    vcf <- vcf[select_variants(od_values < as.numeric(filter.OD),
                               od_values, "OD",
                               "Check that the INFO column of the VCF records OD values."), ]
  }

  info <- vcf@fix[, "INFO"] #Need to get after each filter..

  # Filtering by BIAS
  if (!is.null(filter.BIAS.min) && !is.null(filter.BIAS.max)) {
    info <- vcf@fix[, "INFO"] #Need to get after each filter..
    message("Filtering by BIAS\n")
    bias_values <- extract_info_value(info, "BIAS")
    vcf <- vcf[select_variants(bias_values > as.numeric(filter.BIAS.min) &
                                 bias_values < as.numeric(filter.BIAS.max),
                               bias_values, "BIAS",
                               "Check that the INFO column of the VCF records BIAS values."), ]
  }

  # Filtering by PMC
  if (!is.null(filter.PMC)) {
    info <- vcf@fix[, "INFO"] #Need to get after each filter..
    message("Filtering by PMC\n")
    pmc_values <- extract_info_value(info, "PMC")
    vcf <- vcf[select_variants(pmc_values < as.numeric(filter.PMC),
                               pmc_values, "PMC",
                               "Check that the INFO column of the VCF records PMC values."), ]
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
    # maf() returns NA for a variant with no called genotypes, which the DP and
    # MPP masking above can produce.
    vcf <- vcf[select_variants(maf_df$Frequency > as.numeric(filter.MAF),
                               maf_df$Frequency, "MAF",
                               "This happens when no variant has any called genotypes, which the filter.DP and filter.MPP masking can cause."), ]
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
