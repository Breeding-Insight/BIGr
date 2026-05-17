#' Convert DArTag genotype reports and counts to VCF
#'
#' This function will convert DArT genotype report and Counts files to VCF format
#'
#' This function will convert Allele Dose Report or SNP/INDEL report files and Counts files from DArT into a VCF file.
#' These two files are received directly from DArT for a given sequencing project.
#' SNP/INDEL one-row and two-row reports are treated as diploid genotype reports
#' with 0 = reference homozygote, 1 = alternate homozygote, 2 = heterozygote,
#' and - = missing. Allele Dose reports are interpreted as reference allele
#' dosages using the supplied ploidy.
#' The output file will be saved to the location and with the name that is specified.
#' The VCF format is v4.3
#'
#' @param dart.report Path to the DArT Allele Dose Report or SNP/INDEL report .csv file.
#' @param dart.counts Path to the DArT counts .csv file. Typically contains "Counts" in the file name.
#' @param ploidy The ploidy of the species being analyzed
#' @param output.file output file name and path
#' @return A vcf file
#' @import dplyr
#' @importFrom readr read_csv
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom utils packageVersion read.csv tail write.csv write.table
#' @examples
#' ## Use file paths for each file on the local system
#'
#' #The files are directly from DArT for a given sequencing project.
#' #The are labeled with Dosage_Report or Counts in the file names.
#'
#' #Temp location (only for example)
#' output_file <- tempfile()
#'
#' dosage2vcf(dart.report = system.file("iris_DArT_Allele_Dose_Report_small.csv", package = "BIGr"),
#'            dart.counts = system.file("iris_DArT_Counts_small.csv", package = "BIGr"),
#'            ploidy = 2,
#'            output.file = output_file)
#'
#' # Removing the output for the example
#' rm(output_file)
#'
#' ##The function will output the converted VCF using information from the DArT files
#'
#' @export
dosage2vcf <- function(dart.report, dart.counts, ploidy, output.file) {

  dosage_report <- dart.report
  counts_file <- dart.counts
  ploidy <- as.numeric(ploidy)
  output.file <- paste0(output.file,".vcf")

  message("Reading input files\n")

  #Make a header separate from the dataframe
  vcf_header <- c(
    "##fileformat=VCFv4.3",
    paste0("##BIGr_Dosage2VCF=",packageVersion("BIGr")),
    "##reference=NA",
    "##contig=<ID=NA,length=NA>",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##INFO=<ID=ADS,Number=R,Type=Integer,Description="Depths for the ref and each alt allele in the order listed">',
    '##INFO=<ID=BIAS,Number=1,Type=Float,Description="The estimated allele bias of the SNP from updog">',
    '##INFO=<ID=OD,Number=1,Type=Float,Description="The estimated overdispersion parameter of the SNP from updog">',
    '##INFO=<ID=PMC,Number=1,Type=Float,Description="The estimated proportion of individuals misclassified in the SNP from updog">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, where 1 is the count of alternate alleles">',
    '##FORMAT=<ID=UD,Number=1,Type=Integer,Description="Dosage count of reference alleles from updog, where 0 = homozygous alternate">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=RA,Number=1,Type=Integer,Description="Reference allele read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    '##FORMAT=<ID=MPP,Number=1,Type=Float,Description="Maximum posterior probability for that dosage call from updog">'
  )

  ## Helper functions ---------------------------------------------------------
  is_missing_value <- function(x) {
    x_chr <- trimws(as.character(x))
    is.na(x) | x_chr %in% c("", "-", "NA")
  }

  missing_gt <- function(ploidy) {
    paste(rep(".", ploidy), collapse = "/")
  }

  as_clean_character_matrix <- function(x) {
    mat <- as.matrix(x)
    mat <- matrix(trimws(as.character(mat)), nrow = nrow(mat),
                  ncol = ncol(mat), dimnames = dimnames(mat))
    mat[is_missing_value(mat)] <- NA_character_
    mat
  }

  ensure_unique <- function(x, label) {
    duplicated_values <- unique(x[duplicated(x)])
    if (length(duplicated_values) > 0) {
      stop(label, " must be unique. Duplicated values include: ",
           paste(utils::head(duplicated_values, 5), collapse = ", "))
    }
  }

  check_matching_sets <- function(x, y, x_label, y_label) {
    missing_in_y <- setdiff(x, y)
    missing_in_x <- setdiff(y, x)
    if (length(missing_in_y) > 0 || length(missing_in_x) > 0) {
      msg <- c()
      if (length(missing_in_y) > 0) {
        msg <- c(msg, paste0("in ", x_label, " only: ",
                             paste(utils::head(missing_in_y, 5), collapse = ", ")))
      }
      if (length(missing_in_x) > 0) {
        msg <- c(msg, paste0("in ", y_label, " only: ",
                             paste(utils::head(missing_in_x, 5), collapse = ", ")))
      }
      stop("Mismatched ", x_label, " and ", y_label, " values (", paste(msg, collapse = "; "), ").")
    }
  }

  metadata_cols <- c("MarkerName", "MarkerID", "AlleleID", "CloneID",
                     "AlleleSequence", "AlleleSequenceRef", "AlleleSequenceAlt",
                     "Variant", "CallRate", "OneRatioRef", "OneRatioSnp",
                     "FreqHomRef", "FreqHomSnp", "FreqHets", "PICRef", "PICSnp",
                     "AvgPIC", "AvgCountRef", "AvgCountSnp",
                     "RatioAvgCountRefAvgCountSnp", "Chrom", "ChromPos")

  get_sample_cols <- function(df) {
    sample_cols <- names(df)[!(names(df) %in% metadata_cols)]
    ensure_unique(sample_cols, "Sample columns")
    sample_cols
  }

  parse_variant_bases <- function(variant) {
    variant <- trimws(as.character(variant))
    invalid <- is_missing_value(variant) | !grepl(">", variant, fixed = TRUE)
    if (any(invalid)) {
      stop("Counts Variant values must be present and contain '>' for REF/ALT parsing.")
    }
    clean_variant <- sub("^-:", "", variant)
    parts <- strsplit(clean_variant, ">", fixed = TRUE)
    valid_parts <- vapply(parts, length, integer(1)) == 2
    if (!all(valid_parts)) {
      stop("Counts Variant values must contain exactly one REF>ALT allele definition.")
    }
    data.frame(
      REF = vapply(parts, `[`, character(1), 1),
      ALT = vapply(parts, `[`, character(1), 2),
      stringsAsFactors = FALSE
    )
  }

  derive_coordinates <- function(marker_info) {
    if (all(c("Chrom", "ChromPos") %in% names(marker_info)) &&
        !any(is_missing_value(marker_info$Chrom)) &&
        !any(is_missing_value(marker_info$ChromPos))) {
      return(data.frame(CHROM = as.character(marker_info$Chrom),
                        POS = as.character(marker_info$ChromPos),
                        stringsAsFactors = FALSE))
    }

    marker_names <- as.character(marker_info$MarkerName)
    cannot_split <- !grepl("_", marker_names, fixed = TRUE)
    if (any(cannot_split)) {
      stop("Chrom/ChromPos columns are absent or incomplete, and MarkerName values cannot be split at '_'. Examples: ",
           paste(utils::head(marker_names[cannot_split], 5), collapse = ", "))
    }

    data.frame(
      CHROM = sub("_[^_]*$", "", marker_names),
      POS = sub("^.*_", "", marker_names),
      stringsAsFactors = FALSE
    )
  }

  precompute_genotype_strings <- function(ploidy) {
    genotype_strings <- character(ploidy + 1)
    for (dosage_value in 0:ploidy) {
      ref_count <- dosage_value
      alt_count <- ploidy - dosage_value
      genotype_strings[dosage_value + 1] <- paste(c(rep("0", ref_count), rep("1", alt_count)), collapse = "/")
    }
    genotype_strings
  }

  convert_dosage2gt <- function(dosage_matrix, ploidy) {
    genotype_strings <- precompute_genotype_strings(ploidy)
    genotype_matrix <- matrix(missing_gt(ploidy), nrow = nrow(dosage_matrix), ncol = ncol(dosage_matrix),
                              dimnames = dimnames(dosage_matrix))
    called <- !is.na(dosage_matrix)
    genotype_matrix[called] <- genotype_strings[dosage_matrix[called] + 1]
    genotype_matrix
  }

  matrix_to_character <- function(mat) {
    out <- matrix(as.character(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    out[is.na(mat)] <- "."
    out
  }

  as_numeric_matrix <- function(df, label) {
    char_mat <- as_clean_character_matrix(df)
    numeric_mat <- suppressWarnings(matrix(as.numeric(char_mat), nrow = nrow(char_mat),
                                           ncol = ncol(char_mat), dimnames = dimnames(char_mat)))
    invalid <- !is.na(char_mat) & is.na(numeric_mat)
    if (any(invalid)) {
      stop(label, " contains non-numeric values. Examples: ",
           paste(utils::head(unique(char_mat[invalid]), 5), collapse = ", "))
    }
    numeric_mat
  }

  convert_snp_codes <- function(code_matrix) {
    code_matrix <- as_clean_character_matrix(code_matrix)
    invalid <- !is.na(code_matrix) & !(code_matrix %in% c("0", "1", "2"))
    if (any(invalid)) {
      stop("SNP/INDEL genotype codes must be 0, 1, 2, '-', or NA. Invalid values include: ",
           paste(utils::head(unique(code_matrix[invalid]), 5), collapse = ", "))
    }

    gt_matrix <- matrix("./.", nrow = nrow(code_matrix), ncol = ncol(code_matrix),
                        dimnames = dimnames(code_matrix))
    ud_matrix <- matrix(".", nrow = nrow(code_matrix), ncol = ncol(code_matrix),
                        dimnames = dimnames(code_matrix))

    gt_matrix[!is.na(code_matrix) & code_matrix == "0"] <- "0/0"
    gt_matrix[!is.na(code_matrix) & code_matrix == "1"] <- "1/1"
    gt_matrix[!is.na(code_matrix) & code_matrix == "2"] <- "0/1"
    ud_matrix[!is.na(code_matrix) & code_matrix == "0"] <- "2"
    ud_matrix[!is.na(code_matrix) & code_matrix == "1"] <- "0"
    ud_matrix[!is.na(code_matrix) & code_matrix == "2"] <- "1"

    list(gt = gt_matrix, ud = ud_matrix)
  }

  parse_report <- function(file, ploidy) {
    raw <- suppressMessages(readr::read_csv(file, skip = 6, show_col_types = FALSE))
    raw <- as.data.frame(raw, check.names = FALSE)
    if (nrow(raw) == 0) stop("DArT report file contains no data rows.")

    first_row <- as.character(unlist(raw[1, ], use.names = FALSE))
    is_allele_dose <- length(first_row) >= 5 &&
      identical(first_row[1:5], c("MarkerID", "AvgCountRef", "AvgCountSnp", "Chrom", "ChromPos"))

    if (is_allele_dose) {
      message("DArT report is an Allele Dose report")
      names(raw)[1:5] <- first_row[1:5]
      report <- raw[-1, , drop = FALSE]
      names(report)[names(report) == "MarkerID"] <- "MarkerName"
      sample_cols <- get_sample_cols(report)
      ensure_unique(report$MarkerName, "Report MarkerName")

      dosage_matrix <- as_numeric_matrix(report[, sample_cols, drop = FALSE], "Allele dose report")
      invalid_dosage <- !is.na(dosage_matrix) &
        (dosage_matrix < 0 | dosage_matrix > ploidy | dosage_matrix != round(dosage_matrix))
      if (any(invalid_dosage)) {
        stop("Allele dose values must be integer reference allele counts between 0 and ploidy.")
      }

      rownames(dosage_matrix) <- report$MarkerName
      gt_matrix <- convert_dosage2gt(dosage_matrix, ploidy)
      ud_matrix <- matrix_to_character(dosage_matrix)

      marker_info <- report[, intersect(c("MarkerName", "Chrom", "ChromPos"), names(report)), drop = FALSE]
      rownames(marker_info) <- marker_info$MarkerName

      return(list(type = "allele_dose",
                  marker_info = marker_info,
                  marker_names = report$MarkerName,
                  sample_names = sample_cols,
                  gt = gt_matrix,
                  ud = ud_matrix,
                  variant = NULL,
                  allele_sequence_ref = NULL,
                  allele_sequence_alt = NULL))
    }

    is_snp_1row <- all(c("MarkerName", "AlleleSequenceRef", "AlleleSequenceAlt", "Variant") %in% names(raw))
    is_snp_2row <- all(c("MarkerName", "AlleleSequence", "Variant") %in% names(raw)) &&
      anyDuplicated(raw$MarkerName) > 0

    if (is_snp_1row || is_snp_2row) {
      if (ploidy != 2) {
        stop("SNP/INDEL reports are diploid genotype reports. Use ploidy = 2.")
      }
    }

    if (is_snp_1row) {
      message("DArT report is a SNP/INDEL 1 row report")
      sample_cols <- get_sample_cols(raw)
      ensure_unique(raw$MarkerName, "Report MarkerName")

      code_matrix <- as.matrix(raw[, sample_cols, drop = FALSE])
      rownames(code_matrix) <- raw$MarkerName
      converted <- convert_snp_codes(code_matrix)

      marker_info <- raw[, intersect(c("MarkerName", "Chrom", "ChromPos"), names(raw)), drop = FALSE]
      rownames(marker_info) <- marker_info$MarkerName

      variant <- setNames(trimws(as.character(raw$Variant)), raw$MarkerName)
      allele_sequence_ref <- setNames(as.character(raw$AlleleSequenceRef), raw$MarkerName)
      allele_sequence_alt <- setNames(as.character(raw$AlleleSequenceAlt), raw$MarkerName)

      return(list(type = "snp_1row",
                  marker_info = marker_info,
                  marker_names = raw$MarkerName,
                  sample_names = sample_cols,
                  gt = converted$gt,
                  ud = converted$ud,
                  variant = variant,
                  allele_sequence_ref = allele_sequence_ref,
                  allele_sequence_alt = allele_sequence_alt))
    }

    if (is_snp_2row) {
      message("DArT report is a SNP/INDEL 2 row report")
      sample_cols <- get_sample_cols(raw)
      ref_rows <- is_missing_value(raw$Variant)
      alt_rows <- !ref_rows
      ref_report <- raw[ref_rows, , drop = FALSE]
      alt_report <- raw[alt_rows, , drop = FALSE]
      check_matching_sets(ref_report$MarkerName, alt_report$MarkerName, "SNP ref-row markers", "SNP alt-row markers")
      ensure_unique(ref_report$MarkerName, "SNP ref-row MarkerName")
      ensure_unique(alt_report$MarkerName, "SNP alt-row MarkerName")
      alt_report <- alt_report[match(ref_report$MarkerName, alt_report$MarkerName), , drop = FALSE]

      ref_matrix <- as_clean_character_matrix(ref_report[, sample_cols, drop = FALSE])
      alt_matrix <- as_clean_character_matrix(alt_report[, sample_cols, drop = FALSE])
      invalid <- (!is.na(ref_matrix) & !(ref_matrix %in% c("0", "1"))) |
        (!is.na(alt_matrix) & !(alt_matrix %in% c("0", "1")))
      if (any(invalid)) {
        stop("SNP/INDEL 2 row reports must contain only 0, 1, '-', or NA in sample columns.")
      }

      code_matrix <- matrix(NA_character_, nrow = nrow(ref_matrix), ncol = ncol(ref_matrix),
                            dimnames = list(ref_report$MarkerName, sample_cols))
      called <- !is.na(ref_matrix) & !is.na(alt_matrix)
      code_matrix[called & ref_matrix == "1" & alt_matrix == "0"] <- "0"
      code_matrix[called & ref_matrix == "0" & alt_matrix == "1"] <- "1"
      code_matrix[called & ref_matrix == "1" & alt_matrix == "1"] <- "2"
      code_matrix[called & ref_matrix == "0" & alt_matrix == "0"] <- NA_character_
      converted <- convert_snp_codes(code_matrix)

      marker_info <- data.frame(MarkerName = ref_report$MarkerName, stringsAsFactors = FALSE)
      if (all(c("Chrom", "ChromPos") %in% names(ref_report))) {
        marker_info$Chrom <- ref_report$Chrom
        marker_info$ChromPos <- ref_report$ChromPos
      }
      rownames(marker_info) <- marker_info$MarkerName

      variant <- setNames(trimws(as.character(alt_report$Variant)), alt_report$MarkerName)
      allele_sequence_ref <- setNames(as.character(ref_report$AlleleSequence), ref_report$MarkerName)
      allele_sequence_alt <- setNames(as.character(alt_report$AlleleSequence), alt_report$MarkerName)

      return(list(type = "snp_2row",
                  marker_info = marker_info,
                  marker_names = ref_report$MarkerName,
                  sample_names = sample_cols,
                  gt = converted$gt,
                  ud = converted$ud,
                  variant = variant,
                  allele_sequence_ref = allele_sequence_ref,
                  allele_sequence_alt = allele_sequence_alt))
    }

    stop("Unsupported DArT report format.")
  }

  parse_counts <- function(file) {
    raw <- suppressMessages(readr::read_csv(file, skip = 6, show_col_types = FALSE))
    raw <- as.data.frame(raw, check.names = FALSE)
    if (nrow(raw) == 0) stop("DArT counts file contains no data rows.")

    first_row <- as.character(unlist(raw[1, ], use.names = FALSE))
    if (all(c("MarkerName", "Variant") %in% first_row[seq_len(min(15, length(first_row)))])) {
      message("Counts file contains the counts for the target loci only")
      names(raw)[1:15] <- first_row[1:15]
      counts <- raw[-1, , drop = FALSE]
      sample_cols <- get_sample_cols(counts)

      ref_rows <- is_missing_value(counts$Variant)
      alt_rows <- !ref_rows
      ref_counts <- counts[ref_rows, , drop = FALSE]
      alt_counts <- counts[alt_rows, , drop = FALSE]
      check_matching_sets(ref_counts$MarkerName, alt_counts$MarkerName, "counts ref-row markers", "counts alt-row markers")
      ensure_unique(ref_counts$MarkerName, "Counts ref-row MarkerName")
      ensure_unique(alt_counts$MarkerName, "Counts alt-row MarkerName")
      ref_counts <- ref_counts[match(alt_counts$MarkerName, ref_counts$MarkerName), , drop = FALSE]

      rownames(ref_counts) <- ref_counts$MarkerName
      rownames(alt_counts) <- alt_counts$MarkerName

      bases <- parse_variant_bases(alt_counts$Variant)
      alleles_df <- data.frame(MarkerName = alt_counts$MarkerName,
                               REF = bases$REF,
                               ALT = bases$ALT,
                               CountVariant = trimws(as.character(alt_counts$Variant)),
                               CountAlleleSequenceRef = as.character(ref_counts$AlleleSequence),
                               CountAlleleSequenceAlt = as.character(alt_counts$AlleleSequence),
                               stringsAsFactors = FALSE)
      rownames(alleles_df) <- alleles_df$MarkerName

      return(list(type = "target",
                  alleles_df = alleles_df,
                  marker_names = alt_counts$MarkerName,
                  sample_names = sample_cols,
                  ref_counts = ref_counts[, sample_cols, drop = FALSE],
                  alt_counts = alt_counts[, sample_cols, drop = FALSE]))
    }

    message("Counts file contains the collapsed read counts across all microhaplotypes for the target loci")
    names(raw)[1:5] <- first_row[1:5]
    counts <- raw[-1, , drop = FALSE]
    sample_cols <- get_sample_cols(counts)

    ref_counts <- counts[grepl("Ref$", counts$AlleleID), , drop = FALSE]
    alt_counts <- counts[grepl("Alt$", counts$AlleleID), , drop = FALSE]
    check_matching_sets(ref_counts$CloneID, alt_counts$CloneID, "counts ref-row markers", "counts alt-row markers")
    ensure_unique(ref_counts$CloneID, "Counts ref-row CloneID")
    ensure_unique(alt_counts$CloneID, "Counts alt-row CloneID")
    ref_counts <- ref_counts[match(alt_counts$CloneID, ref_counts$CloneID), , drop = FALSE]

    rownames(ref_counts) <- ref_counts$CloneID
    rownames(alt_counts) <- alt_counts$CloneID
    alleles_df <- data.frame(MarkerName = alt_counts$CloneID,
                             REF = "A",
                             ALT = "B",
                             CountVariant = NA_character_,
                             CountAlleleSequenceRef = NA_character_,
                             CountAlleleSequenceAlt = NA_character_,
                             stringsAsFactors = FALSE)
    rownames(alleles_df) <- alleles_df$MarkerName

    list(type = "collapsed",
         alleles_df = alleles_df,
         marker_names = alt_counts$CloneID,
         sample_names = sample_cols,
         ref_counts = ref_counts[, sample_cols, drop = FALSE],
         alt_counts = alt_counts[, sample_cols, drop = FALSE])
  }

  validate_report_vs_counts <- function(report, counts) {
    if (!is.null(report$variant) && counts$type == "target") {
      report_variant <- report$variant[counts$marker_names]
      count_variant <- counts$alleles_df[counts$marker_names, "CountVariant"]
      mismatch <- !is_missing_value(report_variant) & !is_missing_value(count_variant) &
        trimws(report_variant) != trimws(count_variant)
      if (any(mismatch)) {
        stop("SNP/INDEL report Variant values do not match Counts Variant values. Examples: ",
             paste(utils::head(counts$marker_names[mismatch], 5), collapse = ", "))
      }
    }

    if (!is.null(report$allele_sequence_ref) && counts$type == "target") {
      report_ref <- report$allele_sequence_ref[counts$marker_names]
      report_alt <- report$allele_sequence_alt[counts$marker_names]
      count_ref <- counts$alleles_df[counts$marker_names, "CountAlleleSequenceRef"]
      count_alt <- counts$alleles_df[counts$marker_names, "CountAlleleSequenceAlt"]
      mismatch <- (!is_missing_value(report_ref) & !is_missing_value(count_ref) & report_ref != count_ref) |
        (!is_missing_value(report_alt) & !is_missing_value(count_alt) & report_alt != count_alt)
      if (any(mismatch)) {
        stop("SNP/INDEL report allele sequences do not match Counts allele sequences. Examples: ",
             paste(utils::head(counts$marker_names[mismatch], 5), collapse = ", "))
      }
    }
  }

  ##Get information from DArT Counts and Dosage Report files
  report <- parse_report(dosage_report, ploidy)
  counts <- parse_counts(counts_file)

  check_matching_sets(report$marker_names, counts$marker_names, "report markers", "counts markers")
  check_matching_sets(report$sample_names, counts$sample_names, "report samples", "counts samples")
  validate_report_vs_counts(report, counts)

  marker_order <- report$marker_names
  sample_order <- report$sample_names

  alleles_df <- counts$alleles_df[marker_order, , drop = FALSE]
  marker_info <- report$marker_info[marker_order, , drop = FALSE]
  coordinates <- derive_coordinates(marker_info)
  alleles_df$Chrom <- coordinates$CHROM
  alleles_df$ChromPos <- coordinates$POS

  ref_counts <- counts$ref_counts[marker_order, sample_order, drop = FALSE]
  alt_counts <- counts$alt_counts[marker_order, sample_order, drop = FALSE]
  ref_counts <- as_numeric_matrix(ref_counts, "Reference counts")
  alt_counts <- as_numeric_matrix(alt_counts, "Alternate counts")
  total_counts <- alt_counts + ref_counts

  if (!identical(rownames(report$gt), rownames(ref_counts)) ||
      !identical(colnames(report$gt), colnames(ref_counts))) {
    stop("Internal alignment error: genotype and count matrices are not identically ordered.")
  }

  alleles_df$AltCountsSum <- rowSums(alt_counts)
  alleles_df$RefCountsSum <- rowSums(ref_counts)
  alleles_df$TotalCountSum <- alleles_df$AltCountsSum + alleles_df$RefCountsSum

  vcf_df <- data.frame(
    CHROM = alleles_df$Chrom,
    POS = alleles_df$ChromPos,
    ID = alleles_df$MarkerName,
    REF = alleles_df$REF,
    ALT = alleles_df$ALT,
    QUAL = ".",
    FILTER = ".",
    INFO = NA,
    FORMAT = NA,
    stringsAsFactors = FALSE
  )

  vcf_df$INFO <- paste0("DP=",alleles_df$TotalCountSum,";",
                        "ADS=",alleles_df$RefCountsSum,",",alleles_df$AltCountsSum)
  vcf_df$FORMAT <- paste("GT","UD","DP","RA",sep=":")

  message("Converting dosages to genotype format\n")

  make_vcf_format <- function(gt_matrix, ud_matrix, dp_matrix, ra_matrix) {
    matrix(paste(gt_matrix, ud_matrix, matrix_to_character(dp_matrix), matrix_to_character(ra_matrix), sep = ":"),
           nrow = nrow(gt_matrix), ncol = ncol(gt_matrix), dimnames = dimnames(gt_matrix))
  }

  message("Formatting VCF and generating output file\n")

  # Combine the matrices
  geno_df <- make_vcf_format(report$gt, report$ud, total_counts, ref_counts)

  #Combine the dataframes together
  vcf_df <- cbind(vcf_df,geno_df)

  # Add # to the CHROM column name
  colnames(vcf_df)[1] <- "#CHROM"

  # Write the header to the file
  writeLines(vcf_header, con = output.file)

  # Append the dataframe to the file in tab-separated format
  suppressWarnings(
    write.table(vcf_df, file = output.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
  )
  #Unload all items from memory
  rm(report)
  rm(counts)
  rm(alt_counts)
  rm(ref_counts)
  rm(geno_df)
  rm(vcf_df)
  rm(alleles_df)
  #Clean memory
  gc()

  message("Complete!")
}
