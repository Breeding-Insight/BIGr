#' Format MADC Target Loci Read Counts Into VCF
#'
#' This function will extract the read count information from a MADC file and convert to VCF file format.
#'
#' The DArTag MADC file format is not commonly supported through existing tools. This function
#' will extract the read count information from a MADC file and convert it to a VCF file format for the
#' genotyping panel target markers only
#'
#' @param madc_file Path to MADC file
#' @param output.file output file name and path
#' @return A VCF file v4.3 with the target marker read count information
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom Rdpack reprompt
#' @importFrom reshape2 melt dcast
#' @references
#' Updog R package
#' @export
madc2vcf <- function(madc_file, output.file) {
  #Making the VCF (This is highly dependent on snps being in a format where the SNP IDs are the CHR_POS)

  matrices <- get_countsMADC(madc_file)
  ref_df <- data.frame(matrices[[1]], check.names = FALSE)
  alt_df <- data.frame(matrices[[2]]-matrices[[1]], check.names = FALSE)
  size_df <- data.frame(matrices[[2]], check.names = FALSE)

  #Make the AD column with ref and alt counts
  ad_df <- data.frame(
    mapply(function(ref, alt) {
      paste(ref, alt, sep = ",")
    }, ref_df, alt_df),
    check.names = FALSE
  )
  row.names(ad_df) <- row.names(ref_df)

  #Obtaining Chr and Pos information from the row_names
  new_df <- size_df %>%
    rownames_to_column(var = "row_name") %>%
    separate(row_name, into = c("CHROM", "POS"), sep = "_") %>%
    select(CHROM, POS)

  # Remove leading zeros from the POS column
  new_df$POS <- sub("^0+", "", new_df$POS)

  #Get read count sums
  new_df$TotalRef <- rowSums(ref_df)
  new_df$TotalAlt <- rowSums(alt_df)
  new_df$TotalSize <- rowSums(size_df)

  #Make a header separate from the dataframe
  vcf_header <- c(
    "##fileformat=VCFv4.3",
    paste0("##BIGr_madc2vcf=",packageVersion("BIGr")),
    "##reference=NA",
    "##contig=<ID=NA,length=NA>",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##INFO=<ID=ADS,Number=R,Type=Integer,Description="Depths for the ref and each alt allele in the order listed">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=RA,Number=1,Type=Integer,Description="Reference allele read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
  )

  # Get REF and ALT
  csv <- read.csv(madc_file)
  ref_seq <- csv$AlleleSequence[grep("\\|Ref.*", csv$AlleleID)]
  ref_ord <- csv$CloneID[grep("\\|Ref.*", csv$AlleleID)]
  alt_seq <- csv$AlleleSequence[grep("\\|Alt.*", csv$AlleleID)]
  alt_ord <- csv$CloneID[grep("\\|Alt.*", csv$AlleleID)]

  if(length(ref_seq) != length(alt_seq)) {
    warning("There are missing reference or alternative sequence, the SNP bases could not be recovery.")
  } else {
    if(all(sort(ref_ord) == sort(alt_ord))){
      ref_seq <- ref_seq[order(ref_ord)]
      alt_seq <- alt_seq[order(alt_ord)]

      ref_base <- alt_base <- vector()
      for(i in 1:length(ref_seq)){
        temp_list <- strsplit(c(ref_seq[i], alt_seq[i]), "")
        idx <- temp_list[[1]] != temp_list[[2]]
        ref_base[i] <- temp_list[[1]][idx]
        alt_base[i] <- temp_list[[2]][idx]
      }
    } else {
      warning("There are missing reference or alternative sequence, the SNP bases could not be recovery.")
    }
  }

  #Make the header#Make the VCF df
  vcf_df <- data.frame(
    CHROM = new_df$CHROM,
    POS = new_df$POS,
    ID = row.names(size_df),
    REF = ref_base,
    ALT = alt_base,
    QUAL = ".",
    FILTER = ".",
    INFO = NA,
    FORMAT = NA
  )

  #Add the INFO column for each SNP
  vcf_df$INFO <- paste0("DP=",new_df$TotalSize,";",
                        "ADS=",new_df$TotalRef,",",new_df$TotalAlt)

  #Add the FORMAT label for each SNP
  vcf_df$FORMAT <- paste("DP","RA","AD",sep=":")

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

  # Combine the matrices
  geno_df <- make_vcf_format(as.matrix(size_df), as.matrix(ref_df), as.matrix(ad_df))

  #Combine the dataframes together
  vcf_df <- cbind(vcf_df,geno_df)

  # Add # to the CHROM column name
  colnames(vcf_df)[1] <- "#CHROM"

  # Sort
  vcf_df <- vcf_df[order(vcf_df[,1],as.numeric(as.character(vcf_df[,2]))),]

  # Write the header to the file
  writeLines(vcf_header, con = output.file)

  # Append the dataframe to the file in tab-separated format
  suppressWarnings(
    write.table(vcf_df, file = output.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
  )
  #Unload all items from memory
  rm(matrices)
  rm(ref_df)
  rm(alt_df)
  rm(size_df)
  rm(ad_df)
  rm(vcf_df)
  rm(geno_df)
}
