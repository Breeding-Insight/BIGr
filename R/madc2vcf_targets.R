#' Format MADC Target Loci Read Counts Into VCF
#'
#' Convert DArTag MADC target read counts to a VCF
#'
#' @description
#' Parses a DArTag **MADC** report and writes a **VCF v4.3** containing per-target
#' read counts for the panel’s target loci. This is useful because MADC is not
#' widely supported by general-purpose tools, while VCF is.
#'
#' @details
#' **What this function does**
#' - Runs basic sanity checks on the MADC file (column presence, fixed allele IDs,
#'   IUPAC/ambiguous bases, lowercase bases, indels).
#' - Extracts reference and total read counts per sample and target.
#' - Derives `AD` (ref,alt) by subtraction (alt = total − ref).
#' - If `get_REF_ALT = TRUE`, attempts to recover true REF/ALT bases by comparing
#'   the Ref/Alt probe sequences; targets with >1 polymorphism are discarded.
#' - Optionally accepts a `markers_info` CSV to supply `CHROM`, `POS`, `REF`, `ALT`
#'   (and `Type`, `Indel_pos` when indels are present), bypassing sequence-based
#'   inference.
#'
#' **Output VCF layout**
#' - `INFO` fields:
#'   * `DP`   — total depth across all samples for the locus
#'   * `ADS`  — total counts across samples in the order `ref,alt`
#' - `FORMAT` fields (per sample):
#'   * `DP`   — total reads (ref + alt)
#'   * `RA`   — reads supporting the reference allele
#'   * `AD`   — `"ref,alt"` counts
#'
#' **Strand handling**
#' If a target ID appears in `botloci_file`, its probe sequences are reverse-
#' complemented prior to base comparison so that REF/ALT are reported in the
#' top-strand genomic orientation.
#'
#' **Sanity check behavior**
#' - If required columns or fixed IDs are missing, the function `stop()`s.
#' - If IUPAC/lowercase/indels are detected and `markers_info` is **not**
#'   provided, the function `stop()`s with a diagnostic message explaining what to fix.
#'
#' @param madc_file character. Path to the input MADC CSV file.
#' @param output.file character. Path to the output VCF file to write.
#' @param botloci_file character. Path to a plain-text file listing target IDs
#'   designed on the **bottom** strand (one ID per line). Required when
#'   `get_REF_ALT = TRUE` and `markers_info` is not provided.
#' @param markers_info character or `NULL`. Optional path to a CSV providing target
#'   metadata. Required columns: `BI_markerID, Chr, Pos, Ref, Alt`. If indels are
#'   present, also require `Type, Indel_pos`. When supplied, these values populate
#'   `#CHROM, POS, REF, ALT` in the VCF directly.
#' @param get_REF_ALT logical (default `FALSE`). If `TRUE`, attempts to infer REF/ALT
#'   bases from the Ref/Alt probe sequences in the MADC file (with strand correction
#'   using `botloci_file`). Targets with more than one difference between Ref/Alt
#'   sequences are removed.
#'
#' @return (Invisibly) returns the path to `output.file`. The side effect is a
#'   **VCF v4.3** written to disk containing one row per target and columns for all
#'   samples in the MADC file.
#'
#' @section Dependencies:
#' Uses **dplyr**, **tidyr**, **tibble**, **reshape2**, **Biostrings** and base
#' **utils**. Helper functions expected in this package: `check_madc_sanity()`,
#' `get_countsMADC()`, `get_counts()`, and `check_botloci()`.
#'
#' @examples
#' # Example files shipped with the package
#' madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
#' bot_file  <- system.file("example_SNPs_DArTag-probe-design_f180bp.botloci",
#'                          package = "BIGr")
#' out_vcf <- tempfile(fileext = ".vcf")
#'
#' # Convert MADC to VCF (attempting to recover REF/ALT from probe sequences)
#' \dontrun{
#' madc2vcf_targets(
#'   madc_file    = madc_file,
#'   output.file  = out_vcf,
#'   botloci_file = bot_file,
#'   get_REF_ALT  = TRUE
#' )
#' }
#'
#' # Clean up (example)
#' unlink(out_vcf)
#'
#' @seealso
#' `check_madc_sanity()`, `get_countsMADC()`, `check_botloci()`
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom reshape2 melt dcast
#' @importFrom utils write.table
#' @importFrom Biostrings DNAString reverseComplement
#'
#' @export
madc2vcf_targets <- function(madc_file,
                             output.file,
                             botloci_file,
                             markers_info = NULL,
                             get_REF_ALT = FALSE) {

  # MADC checks
  report <- read.csv(madc_file)
  checks <- check_madc_sanity(report)

  messages_results <- mapply(function(check, message) {
    if (check)  message[1] else message[2]
  }, checks$checks, checks$messages)

  if(any(!(checks$checks[c("Columns", "FixAlleleIDs")]))){
    idx <- which(!(checks$checks[c("Columns", "FixAlleleIDs")]))
    stop(paste("The MADC file does not pass the sanity checks:\n",
               paste(messages_results[c("Columns", "FixAlleleIDs")[idx]], collapse = "\n")))
  }

  if(any(checks$checks[c("IUPACcodes", "LowerCase", "Indels")])){
    idx <- which((checks$checks[c("IUPACcodes", "LowerCase", "Indels")]))
    if(is.null(markers_info)) stop(paste(messages_results[c("IUPACcodes", "LowerCase", "Indels")[idx]], collapse = "\n"))
  }

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
  if(is.null(markers_info)){
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

    # Get REF and ALT
    if(get_REF_ALT){
      if(is.null(botloci_file)) stop("Please provide the botloci file to recover the reference and alternative bases.")
      csv <- get_counts(madc_file)
      # Keep only the ones that have alt and ref
      csv <- csv[which(csv$CloneID %in% rownames(ad_df)),]

      # Get reverse complement the tag is present in botloci
      botloci <- read.table(botloci_file, header = FALSE)

      # Check if the botloci file marker IDs match with the MADC file
      checked_botloci <- check_botloci(botloci, csv)
      botloci <- checked_botloci[[1]]
      csv <- checked_botloci[[2]]

      # FIXED: Store original sequences before any transformation
      csv$OriginalAlleleSequence <- csv$AlleleSequence

      # Apply reverse complement to sequences for bottom strand markers
      idx <- which(csv$CloneID %in% botloci[,1])
      csv$AlleleSequence[idx] <- sapply(csv$AlleleSequence[idx], function(sequence) as.character(reverseComplement(DNAString(sequence))))

      ref_seq <- csv$AlleleSequence[grep("\\|Ref.*", csv$AlleleID)]
      ref_ord <- csv$CloneID[grep("\\|Ref.*", csv$AlleleID)]
      alt_seq <- csv$AlleleSequence[grep("\\|Alt.*", csv$AlleleID)]
      alt_ord <- csv$CloneID[grep("\\|Alt.*", csv$AlleleID)]

      # FIXED: Get original sequences for SNP calling
      orig_ref_seq <- csv$OriginalAlleleSequence[grep("\\|Ref.*", csv$AlleleID)]
      orig_alt_seq <- csv$OriginalAlleleSequence[grep("\\|Alt.*", csv$AlleleID)]

      if(all(sort(ref_ord) == sort(alt_ord))){
        # Order sequences consistently
        ref_seq <- ref_seq[order(ref_ord)]
        alt_seq <- alt_seq[order(alt_ord)]
        orig_ref_seq <- orig_ref_seq[order(ref_ord)]
        orig_alt_seq <- orig_alt_seq[order(alt_ord)]
        ordered_clone_ids <- sort(ref_ord)

        ref_base <- alt_base <- vector()
        for(i in 1:length(orig_ref_seq)){
          # FIXED: Use original sequences for SNP calling
          temp_list <- strsplit(c(orig_ref_seq[i], orig_alt_seq[i]), "")
          idx_diff <- which(temp_list[[1]] != temp_list[[2]])

          if(length(idx_diff) > 1) { # If finds more than one polymorphism between Ref and Alt sequences
            ref_base[i] <- NA
            alt_base[i] <- NA
          } else if(length(idx_diff) == 1) {
            orig_ref_base <- temp_list[[1]][idx_diff]
            orig_alt_base <- temp_list[[2]][idx_diff]

            # FIXED: Apply reverse complement to bases only if marker is in botloci
            if(ordered_clone_ids[i] %in% botloci[,1]) {
              ref_base[i] <- as.character(reverseComplement(DNAString(orig_ref_base)))
              alt_base[i] <- as.character(reverseComplement(DNAString(orig_alt_base)))
            } else {
              ref_base[i] <- orig_ref_base
              alt_base[i] <- orig_alt_base
            }
          } else {
            # No differences found
            ref_base[i] <- NA
            alt_base[i] <- NA
          }
        }
      } else {
        warning("There are missing reference or alternative sequence, the SNP bases could not be recovery.")
        ref_base <- "."
        alt_base <- "."
      }

    } else {
      ref_base <- "."
      alt_base <- "."
    }
  } else {
    # Verify markers_info file
    df <- read.csv(markers_info)
    if(checks$checks["Indels"]){
      if(!all(c("BI_markerID","Chr","Pos","Ref","Alt","Type", "Indel_pos") %in% colnames(df)))
        stop("The markers_info dataframe must contain the following columns: BI_markerID, CHROM, POS, REF, ALT, Type, Indel_pos")
    }
    if(!all(c("BI_markerID","Chr","Pos","Ref","Alt") %in% colnames(df)))
      stop("The markers_info dataframe must contain the following columns: BI_markerID, CHROM, POS, REF, ALT")

    if(!all(rownames(ad_df)%in% df$BI_markerID))
      warning("Not all MADC CloneID was found in the markers_info file. These markers will be removed.")

    matched <- df[match(rownames(ad_df), df$BI_markerID),]

    new_df <- data.frame(
      CHROM = matched$Chr,
      POS = matched$Pos
    )

    #Get read count sums
    new_df$TotalRef <- rowSums(ref_df)
    new_df$TotalAlt <- rowSums(alt_df)
    new_df$TotalSize <- rowSums(size_df)

    ref_base <- matched$Ref
    alt_base <- matched$Alt
  }

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

  if(sum(is.na(vcf_df$REF)) >1) {
    warning(paste("Markers removed because of presence of more than one polymorphism between ref and alt sequences:",sum(is.na(vcf_df$REF))))
    vcf_df <- vcf_df[-which(is.na(vcf_df$REF)),]
  }

  # Write the header to the file
  writeLines(vcf_header, con = output.file)

  # Append the dataframe to the file in tab-separated format
  suppressWarnings(
    write.table(vcf_df, file = output.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
  )
}
