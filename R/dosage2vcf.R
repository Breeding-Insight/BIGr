#' Convert DArTag Dosage and Counts to VCF
#'
#' This function will convert the DArT Dosage Report and Counts files to VCF format
#'
#' This function will convert the Dosage Report and Counts files from DArT into a VCF file.
#' These two files are received directly from DArT for a given sequencing project.
#' The output file will be saved to the location and with the name that is specified.
#' The VCF format is v4.3
#'
#' @param dart.report Path to the DArT dosage report .csv file. Typically contains "Dosage Report" in the file name.
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

  ##Get information from DArT Counts and Dosage Report files

  dosage <- suppressMessages(readr::read_csv(dosage_report,
                                      skip = 6,show_col_types = FALSE))
  colnames(dosage)[1:5] <- dosage[1,1:5]
  dosage <- dosage[-1, ]
  names(dosage)[names(dosage) == "MarkerID"] <- "MarkerName"
  dosage <- as.data.frame(dosage)
  row.names(dosage) <- dosage$MarkerName

  counts <- suppressMessages(readr::read_csv(counts_file,
                                      skip = 6,show_col_types = FALSE))

  #Check that the counts file is in the 2 row format
  #if (anyDuplicated(counts$MarkerName) > 0) { #I am going to just check if there are 2x rows in counts as dosage report
  if (nrow(counts) >= nrow(dosage)*2) {
    message("Note: Counts file is in the 2 row format for Ref and Alt alleles")
  } else {
    stop("Counts file is in single row format. Only two row format (row for ref and alt allele) is currently supported.")
    return()
  }

  #Parse counts file depending on if its the collapsed format or not
  if (all(c("MarkerName", "Variant") %in% counts[1, 1:15])) {
    message("Counts file contains the counts for the target loci only")

    colnames(counts)[1:15] <- counts[1,1:15]
    counts <- counts[-1, ]
    counts <- as.data.frame(counts)

    #Get ref counts dataframe
    ref_counts <- counts[is.na(counts$Variant), ]
    row.names(ref_counts) <- ref_counts$MarkerName
    #Get alt counts dataframe
    alt_counts <- counts[!is.na(counts$Variant), ]
    row.names(alt_counts) <- alt_counts$MarkerName

    #unload counts file
    rm(counts)

    #Get the Ref and Alt allele from the alt_counts file
    ##Note sometimes there are multiple nucleotides, I am assuming this file also includes indels, but make sure this is not an error
    alleles_df <- alt_counts %>%
      rowwise() %>%  # Apply operations to each row
      mutate(Variant = gsub("-:", "",Variant),  # Remove "-:" prefix
             REF = strsplit(Variant, ">")[[1]][1],
             ALT = strsplit(Variant, ">")[[1]][2]) %>%
      ungroup() %>%
      select(MarkerName, REF, ALT)

    #Add the CHROM and POS information to the alleles_df from the dosage report
    alleles_df <- merge(alleles_df, dosage[, c("MarkerName","Chrom", "ChromPos")], by = "MarkerName", all.x = TRUE)
    #Add row name
    row.names(alleles_df) <- alleles_df$MarkerName

  }else{
    message("Counts file contains the collapsed read counts across all microhaplotypes for the target loci")

    colnames(counts)[1:5] <- counts[1,1:5]
    counts <- counts[-1, ]
    counts <- as.data.frame(counts)

    #Get ref counts dataframe
    ref_counts <- counts[grepl("Ref$", counts$AlleleID), ]
    row.names(ref_counts) <- ref_counts$CloneID
    #Get alt counts dataframe
    alt_counts <- counts[grepl("Alt$", counts$AlleleID), ]
    row.names(alt_counts) <- alt_counts$CloneID

    #unload counts file
    rm(counts)

    #Get the Ref and Alt allele from the alt_counts file
    ##Note sometimes there are multiple nucleotides, I am assuming this file also includes indels, but make sure this is not an error
    alleles_df <- data.frame(MarkerName = alt_counts$CloneID,
                             REF = "A",
                             ALT = "B")

    #Add the CHROM and POS information to the alleles_df from the dosage report
    alleles_df <- merge(alleles_df, dosage[, c("MarkerName","Chrom", "ChromPos")], by = "MarkerName", all.x = TRUE)
    #Add row name
    row.names(alleles_df) <- alleles_df$MarkerName

  }

  #Remove the unwanted information from the counts dataframes
  cols_to_remove <- c("MarkerName","AlleleSequence","Variant",
                      "CallRate","OneRatioRef","OneRatioSnp","FreqHomSnp",
                      "FreqHets","PICRef","PICSnp","AvgPIC","FreqHomRef",
                      "AvgCountRef","AvgCountSnp","RatioAvgCountRefAvgCountSnp",
                      "AlleleID","CloneID")
  alt_counts <- alt_counts[, !(colnames(alt_counts) %in% cols_to_remove)]
  ref_counts <- ref_counts[, !(colnames(ref_counts) %in% cols_to_remove)]

  #Ensure that they are in the same order
  ref_counts <- ref_counts[row.names(alt_counts),]

  #Make the total counts dataframe
  total_counts <- alt_counts + ref_counts

  #Get the total ref, total alt, and total read depth for each marker
  alleles_df$AltCountsSum <- rowSums(alt_counts)[rownames(alleles_df)] #Alt
  alleles_df$RefCountsSum <- rowSums(ref_counts)[rownames(alleles_df)] #Ref
  alleles_df$TotalCountSum <- alleles_df$AltCountsSum + alleles_df$RefCountsSum

  #Remove the unwanted information from the dosage dataframe
  d_cols_to_remove <- c("MarkerName","AvgCountRef","AvgCountSnp","Chrom","ChromPos")
  dosage <- dosage[, !(colnames(dosage) %in% d_cols_to_remove)]

  #Make the VCF df
  vcf_df <- data.frame(
    CHROM = alleles_df$Chrom,
    POS = alleles_df$ChromPos,
    ID = alleles_df$MarkerName,
    REF = alleles_df$REF,
    ALT = alleles_df$ALT,
    QUAL = ".",
    FILTER = ".",
    INFO = NA,
    FORMAT = NA
  )

  #Add the INFO column for each SNP
  vcf_df$INFO <- paste0("DP=",alleles_df$TotalCountSum,";",
                        "ADS=",alleles_df$RefCountsSum,",",alleles_df$AltCountsSum)

  #Add the FORMAT label for each SNP
  vcf_df$FORMAT <- paste("GT","UD","DP","RA",sep=":")

  message("Converting dosages to genotype format\n")

  ###Convert genotypes from dosage to gt
  # Precompute genotype strings for all possible dosage values to improve efficiency
  precompute_genotype_strings <- function(ploidy) {
    genotype_strings <- character(ploidy + 1)
    # Generate the genotype string based on the dosage and ploidy
    # Updog uses the ref counts, which is not typical, so this corrects it
    for (dosage in 0:ploidy) {
      ref_count <- dosage
      alt_count <- ploidy - dosage
      genotype_strings[dosage + 1] <- paste(c(rep("0", ref_count), rep("1", alt_count)), collapse = "/")
    }
    return(genotype_strings)
  }

  # Apply the precomputed genotype strings to the matrix
  convert_dosage2gt <- function(dosage_matrix, ploidy) {
    dosage_matrix <- as.matrix(dosage_matrix)
    genotype_strings <- precompute_genotype_strings(ploidy)

    # Handle missing values separately
    genotype_matrix <- matrix(genotype_strings[dosage_matrix + 1], nrow = nrow(dosage_matrix), ncol = ncol(dosage_matrix))
    genotype_matrix[is.na(dosage_matrix)] <- "./." # Handle missing values

    # Retain row and column names
    rownames(genotype_matrix) <- rownames(dosage_matrix)
    colnames(genotype_matrix) <- colnames(dosage_matrix)

    return(genotype_matrix)
  }

  # Convert the dosage matrix to genotypes
  geno_df <- convert_dosage2gt(dosage, ploidy)

  #Combine info from the matrices to form the VCF information for each sample
  # Combine the matrices into a single matrix with elements separated by ":"
  make_vcf_format <- function(..., separator = ":") {
    matrices <- list(...)
    n <- length(matrices)

    # Convert matrices to long form
    long_forms <- lapply(matrices, function(mat) {
      suppressMessages(reshape2::melt(mat, varnames = c("Row", "Col"), value.name = "Value"))
    })

    # Concatenate the elements
    combined_long <- long_forms[[1]]
    combined_long$Combined <- combined_long$Value

    for (i in 2:n) {
      combined_long$Combined <- paste(combined_long$Combined, long_forms[[i]]$Value, sep = separator)
    }

    # Convert back to wide form
    combined_wide <- suppressMessages(reshape2::dcast(combined_long, Row ~ Col, value.var = "Combined"))

    # Restore row and column names
    rownames(combined_wide) <- combined_wide$Row
    combined_wide$Row <- NULL
    colnames(combined_wide) <- colnames(matrices[[1]])

    return(as.matrix(combined_wide))
  }

  message("Formatting VCF and generating output file\n")

  # Combine the matrices
  geno_df <- make_vcf_format(geno_df, dosage, total_counts, ref_counts)

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
  rm(dosage)
  rm(alt_counts)
  rm(ref_counts)
  rm(geno_df)
  rm(vcf_df)
  rm(alleles_df)
  #Clean memory
  gc()

  message("Complete!")
}
