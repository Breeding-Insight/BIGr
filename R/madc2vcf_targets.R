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
#'   metadata. Required columns: `CloneID, Chr, Pos, Ref, Alt`. This file is required in
#'   case your MADC CloneID column doesn't have the format CHR_POS. If indels are
#'   present, columns `Type, Indel_pos` are also required.
#' @param get_REF_ALT logical (default `FALSE`). If `TRUE`, attempts to infer REF/ALT
#'   bases from the Ref/Alt probe sequences in the MADC file (with strand correction
#'   using `botloci_file`). Targets with more than one difference between Ref/Alt
#'   sequences are removed.
#' @param collapse_matches_counts logical (default `FALSE`). If `TRUE`, counts for targets with identical `CHROM_POS` are summed together. This is useful when the MADC file contains multiple rows per target (e.g., due to multiple alleles or technical replicates) and you want to aggregate them into a single entry per unique target.
#' @param verbose logical (default `FALSE`). If `TRUE`, prints detailed messages about
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
                             get_REF_ALT = FALSE,
                             collapse_matches_counts = FALSE,
                             verbose = TRUE) {

  vmsg("Checking inputs", verbose = verbose, level = 0, type = ">>")

  # Input checks
  if(!file.exists(madc_file)) stop("The MADC file does not exist.")
  if(!is.character(output.file)) stop("output.file must be a character string.")
  if(get_REF_ALT && is.null(botloci_file)) stop("Please provide the botloci file to recover the reference and alternative bases.")
  if(get_REF_ALT && !file.exists(botloci_file)) stop("The botloci file does not exist.")
  if(!is.null(markers_info) && !file.exists(markers_info)) stop("The markers_info file does not exist.")
  if(!is.null(markers_info) && !is.character(markers_info)) stop("markers_info must be a character string or NULL.")
  if(!is.logical(get_REF_ALT)) stop("get_REF_ALT must be a logical value (TRUE or FALSE).")
  if(!is.logical(verbose)) stop("verbose must be a logical value (TRUE or FALSE).")

  # Create a VCF header line with metadata about the command and its parameters
  bigr_meta <- paste0('##BIGrCommandLine.madc2vcf_targets=<ID=madc2vcf_targets,Version="',
                      packageVersion("BIGr"), '",Data="',
                      Sys.time(),'", CommandLine="> madc2vcf_targets(',deparse(substitute(madc)),', ',
                      "output.file= ", output.file, ', ',
                      "botloci_file= ", botloci_file, ', ',
                      "markers_info= ", markers_info, ', ',
                      "get_REF_ALT= ", get_REF_ALT, ', ',
                      "verbose= ", verbose,')">')

  # MADC checks
  report <- read.csv(madc_file)
  checks <- check_madc_sanity(report)

  messages_results <- mapply(function(check, message) {
    if (check)  message[1] else message[2]
  }, checks$checks, checks$messages)

  for(i in seq_along(messages_results))
    vmsg(messages_results[i], verbose = verbose, level = 1, type = ">>")

  if(any(!(checks$checks[c("Columns", "FixAlleleIDs")]))){
    idx <- which(!(checks$checks[c("Columns", "FixAlleleIDs")]))
    stop(paste("The MADC file does not pass the sanity checks:\n",
               paste(messages_results[c("Columns", "FixAlleleIDs")[idx]], collapse = "\n")))
  }

  if(any(checks$checks[c("IUPACcodes", "LowerCase", "Indels")])){
    idx <- which((checks$checks[c("IUPACcodes", "LowerCase", "Indels")]))
    if(is.null(markers_info)) stop("Please provide a markers_info file to proceed. The MADC file does not pass the sanity checks:\n",
                                  paste(messages_results[c("IUPACcodes", "LowerCase", "Indels")[idx]], collapse = "\n"))
    else vmsg("MADC file has some issues (IUPAC codes, lowercase bases, indels), but a markers_info file is provided, so proceeding with VCF generation.", verbose = verbose, level = 1, type = ">>")
  }

  # Check marker names compatibility between MADC and botloci
  if(!is.null(botloci_file)){
    botloci <- read.csv(botloci_file, header = F)
    checked_botloci <- check_botloci(botloci, report)
    botloci <- checked_botloci[[1]]
    report <- checked_botloci[[2]]
  }

  vmsg("Input checks done", verbose = verbose, level = 1, type = ">>")

  vmsg("Extracting depth information", verbose = verbose, level = 0, type = ">>")

  matrices <- get_countsMADC(madc_object = report, collapse_matches_counts = collapse_matches_counts, verbose = verbose)
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

  vmsg("Depth information extracted", verbose = verbose, level = 1, type = ">>")

  #Obtaining Chr and Pos information from the row_names
  if(is.null(markers_info)){
    vmsg("No markers_info file provided. Attempting to recover CHROM and POS from CloneID...", verbose = verbose, level = 0, type = ">>")
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

    vmsg("CHROM and POS recovered from CloneID", verbose = verbose, level = 1, type = ">>")
    # Get REF and ALT
    if(get_REF_ALT){
      vmsg("get_REF_ALT = TRUE. Attempting to recover REF and ALT bases from probe sequences...", verbose = verbose, level = 0, type = ">>")
      csv <- get_counts(madc_object = report, collapse_matches_counts = collapse_matches_counts, verbose = FALSE)
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

        more_poly <- no_diff <- 0
        ref_base <- alt_base <- vector()
        for(i in seq_along(orig_ref_seq)){
          # FIXED: Use original sequences for SNP calling
          temp_list <- strsplit(c(orig_ref_seq[i], orig_alt_seq[i]), "")
          idx_diff <- which(temp_list[[1]] != temp_list[[2]])

          if(length(idx_diff) > 1) { # If finds more than one polymorphism between Ref and Alt sequences
            ref_base[i] <- NA
            alt_base[i] <- NA
            more_poly <- more_poly + 1
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
            no_diff <- no_diff + 1
          }
        }
        if(more_poly > 0) vmsg(paste(more_poly, "markers removed because more than one polymorphism was found between Ref and Alt sequences"), verbose = verbose, level = 2, type = ">>")
        if(no_diff > 0) vmsg(paste(no_diff, "markers removed because no differences were found between Ref and Alt sequences"), verbose = verbose, level = 2, type = ">>")

      } else {
        ref_base <- "."
        alt_base <- "."
        vmsg(paste("REF and ALT bases could not be recovered because of missing reference or alternative sequences"), verbose = verbose, level = 1, type = ">>")
      }
    } else {
      ref_base <- "."
      alt_base <- "."
      vmsg(paste("REF and ALT bases not recovered because get_REF_ALT = FALSE"), verbose = verbose, level = 1, type = ">>")
    }
  } else {
    vmsg("markers_info file provided. Using CHROM, POS, REF and ALT from the file.", verbose = verbose, level = 0, type = ">>")
    # Verify markers_info file
    df <- read.csv(markers_info)

    # Accept either CloneID or BI_markerID as the marker ID column
    if ("BI_markerID" %in% colnames(df)) {
      id_col <- "BI_markerID"
    } else if ("CloneID" %in% colnames(df)) {
      id_col <- "CloneID"
    } else {
      stop("The markers_info file must contain a marker ID column named either 'CloneID' or 'BI_markerID'.")
    }

    if(checks$checks["Indels"]){
      vmsg("Indels detected in MADC file. Checking for required columns in markers_info...", verbose = verbose, level = 1, type = ">>")
      if(!all(c(id_col,"Chr","Pos","Ref","Alt","Type", "Indel_pos") %in% colnames(df)))
        stop(paste0("The markers_info dataframe must contain the following columns: ", id_col, ", Chr, Pos, Ref, Alt, Type, Indel_pos"))
    }
    if(!all(c(id_col,"Chr","Pos","Ref","Alt") %in% colnames(df)))
      stop(paste0("The markers_info dataframe must contain the following columns: ", id_col, ", Chr, Pos, Ref, Alt"))

    if(!all(rownames(ad_df) %in% df[[id_col]])){
      miss_CloneIDs <- rownames(ad_df)[!rownames(ad_df) %in% df[[id_col]]]
      vmsg(paste("Not all MADC CloneID was found in the markers_info file. These markers will be removed:", paste(miss_CloneIDs, collapse = " ")), verbose = verbose, level = 2, type = ">>")
      warning("Not all MADC CloneID was found in the markers_info file. These markers will be removed.")
    }
    matched <- df[match(rownames(ad_df), df[[id_col]]),]

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

  vmsg("CHROM, POS, REF and ALT columns prepared", verbose = verbose, level = 1, type = ">>")

  vmsg("Preparing VCF dataframe", verbose = verbose, level = 0, type = ">>")
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
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    bigr_meta
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
    vmsg("Sample columns formatted for VCF", verbose = verbose, level = 1, type = ">>")

    return(as.matrix(combined_wide))
  }

  # Combine the matrices
  geno_df <- make_vcf_format(as.matrix(size_df), as.matrix(ref_df), as.matrix(ad_df))

  #Combine the dataframes together
  vcf_df <- cbind(vcf_df,geno_df)
  vmsg("VCF dataframe prepared", verbose = verbose, level = 1, type = ">>")

  # Add # to the CHROM column name
  colnames(vcf_df)[1] <- "#CHROM"

  # Sort
  vcf_df <- vcf_df[order(vcf_df[,1],as.numeric(as.character(vcf_df[,2]))),]

  if(sum(is.na(vcf_df$REF)) >1) {
    vmsg(paste(sum(is.na(vcf_df$REF)), "markers removed because of presence of more than one polymorphism between ref and alt sequences."), verbose = verbose, level = 1, type = ">>")
    warning(paste("Markers removed because of presence of more than one polymorphism between ref and alt sequences:",sum(is.na(vcf_df$REF))))
    vcf_df <- vcf_df[-which(is.na(vcf_df$REF)),]
  }

  # Write the header to the file
  writeLines(vcf_header, con = output.file)

  # Append the dataframe to the file in tab-separated format
  suppressWarnings(
    write.table(vcf_df, file = output.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
  )
  vmsg(paste("VCF file written to", output.file), verbose = verbose, level = 0, type = ">>")
}
