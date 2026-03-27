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
#' - Runs basic sanity checks on the MADC file via `check_madc_sanity()` (column
#'   presence, fixed allele IDs, IUPAC/ambiguous bases, lowercase bases, indels,
#'   chromosome/position format, all-NA rows/columns, Ref/Alt sequence presence).
#' - Extracts reference and total read counts per sample and target.
#' - Derives `AD` (ref,alt) by subtraction (alt = total − ref).
#' - If `get_REF_ALT = TRUE`, recovers REF/ALT bases either from `markers_info`
#'   (when `Ref`/`Alt` columns are present) or by comparing the Ref/Alt probe
#'   sequences in the MADC file (with strand correction via `botloci_file`).
#'   Targets with >1 polymorphism between sequences are discarded.
#' - Optionally accepts a `markers_info` CSV to supply `CHROM`, `POS`, `REF`,
#'   `ALT`, bypassing sequence-based inference.
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
#' **Sanity check behaviour and requirements**
#'
#' The function always stops if required columns (`CloneID`, `AlleleID`,
#' `AlleleSequence`) are missing.
#'
#' For the remaining checks the required inputs depend on the combination of
#' check result and `get_REF_ALT`:
#'
#' | Check | Status | `get_REF_ALT` | Required |
#' |---|---|---|---|
#' | **IUPAC codes** | detected | `TRUE` | `markers_info` with `Ref`/`Alt` |
#' | | detected | `FALSE` | — |
#' | | not detected | `TRUE` | `botloci_file` **or** `markers_info` with `Ref`/`Alt` |
#' | | not detected | `FALSE` | — |
#' | **Indels** | detected | `TRUE` | `markers_info` with `Ref`/`Alt` |
#' | | detected | `FALSE` | — |
#' | | not detected | `TRUE` | `botloci_file` **or** `markers_info` with `Ref`/`Alt` |
#' | | not detected | `FALSE` | — |
#' | **ChromPos** | valid | `TRUE` | `botloci_file` **or** `markers_info` with `Ref`/`Alt` |
#' | | valid | `FALSE` | — |
#' | | invalid | `TRUE` | `markers_info` with `Chr`/`Pos`/`Ref`/`Alt` **or** `markers_info` with `Chr`/`Pos` + `botloci_file` |
#' | | invalid | `FALSE` | `markers_info` with `Chr`/`Pos` |
#' | **FixAlleleIDs** | fixed | `TRUE` | `botloci_file` **or** `markers_info` with `Ref`/`Alt` |
#' | | fixed | `FALSE` | — |
#' | | not fixed | `TRUE` | `markers_info` with `Ref`/`Alt` |
#' | | not fixed | `FALSE` | — |
#'
#' @param madc_file character. Path to the input MADC CSV file.
#' @param output.file character. Path to the output VCF file to write.
#' @param botloci_file character or `NULL` (default `NULL`). Path to a plain-text
#'   file listing target IDs designed on the **bottom** strand (one ID per line).
#'   Used for strand-correcting probe sequences when `get_REF_ALT = TRUE` and
#'   `markers_info` does not supply `Ref` and `Alt` columns. Also required when
#'   `ChromPos` is invalid and `markers_info` does not provide `Ref`/`Alt`.
#' @param markers_info character or `NULL`. Optional path to a CSV providing target
#'   metadata. Accepted columns:
#'   - `CloneID` or `BI_markerID` (required as marker identifier);
#'   - `Chr`, `Pos` — required when `CloneID` does not follow the `Chr_Pos` format;
#'   - `Ref`, `Alt` — required when `get_REF_ALT = TRUE` and probe-sequence
#'     inference is not possible (IUPAC codes, indels, or unfixed allele IDs).
#' @param get_REF_ALT logical (default `FALSE`). If `TRUE`, attempts to recover
#'   REF/ALT bases. The source is chosen automatically: `markers_info` `Ref`/`Alt`
#'   columns take priority; otherwise probe sequences from the MADC are compared
#'   (with `botloci_file` for strand correction). Targets with more than one
#'   difference between Ref/Alt sequences are removed. When `FALSE`, REF and ALT
#'   are set to `"."` in the output VCF.
#' @param collapse_matches_counts logical (default `FALSE`). If `TRUE`, counts for
#'   `|AltMatch` and `|RefMatch` rows are summed into their corresponding `|Ref`
#'   and `|Alt` rows before building the matrices. Useful when the MADC contains
#'   multiple allele rows per target that should be aggregated.
#' @param verbose logical (default `TRUE`). If `TRUE`, prints detailed progress
#'   messages about each processing step.
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
                             botloci_file = NULL,
                             markers_info = NULL,
                             get_REF_ALT = FALSE,
                             collapse_matches_counts = FALSE,
                             verbose = TRUE) {

  vmsg("Checking inputs", verbose = verbose, level = 0, type = ">>")

  # Input checks
  if(!(file.exists(madc_file) | url_exists(madc_file))) stop("The MADC file does not exist.")
  if(!is.character(output.file)) stop("output.file must be a character string.")
  if(!is.null(markers_info) && !is.character(markers_info)) stop("markers_info must be a character string or NULL.")
  if(!is.null(markers_info) && !file.exists(markers_info) && !url_exists(markers_info)) stop("The markers_info file does not exist.")
  if(!is.logical(get_REF_ALT)) stop("get_REF_ALT must be a logical value (TRUE or FALSE).")
  if(!is.logical(verbose)) stop("verbose must be a logical value (TRUE or FALSE).")

  # Create a VCF header line with metadata about the command and its parameters
  bigr_meta <- paste0('##BIGrCommandLine.madc2vcf_targets=<ID=madc2vcf_targets,Version="',
                      packageVersion("BIGr"), '",Data="',
                      Sys.time(),'", CommandLine="> madc2vcf_targets(',deparse(substitute(madc_file)),', ',
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

  if(any(!(checks$checks[c("Columns")]))){
    idx <- which(!(checks$checks[c("Columns")]))
    if(length(idx) > 0)
    stop(paste("The MADC file does not pass the sanity checks:\n",
               paste(messages_results[c("Columns")[idx]], collapse = "\n")))
  }

  if(any(checks$checks[c("IUPACcodes", "Indels")]) && get_REF_ALT){
    idx <- which((checks$checks[c("IUPACcodes", "Indels")]))
    if(is.null(markers_info)) stop(paste("Please provide a markers_info file to proceed. The MADC file does not pass the sanity checks:",
                                  paste(messages_results[c("IUPACcodes", "Indels")[idx]], collapse = "\n")))
    else vmsg("MADC file has IUPAC codes and/or indels, but a markers_info file is provided, so proceeding with VCF generation.", verbose = verbose, level = 1, type = ">>")
  }

  if(checks$checks["LowerCase"]){
    vmsg("MADC Allele Sequences presented lower case characters. They were converted to upper case.", verbose = verbose, level = 1)
    report$AlleleSequence <- toupper(report$AlleleSequence)
  }

  # ---- Validate botloci and pre-process CloneIDs based on get_REF_ALT logic ----
  mi_df          <- NULL   # markers_info data frame (loaded once, reused below)
  mi_has_ref_alt <- FALSE  # TRUE when markers_info provides Ref and Alt columns
  botloci        <- NULL   # botloci data frame (set when needed)

  # Check whether markers_info is present and contains Ref + Alt columns
  if(!is.null(markers_info)) {
    mi_df          <- read.csv(markers_info)
    mi_has_ref_alt <- all(c("Ref", "Alt") %in% colnames(mi_df))
  }

  if(!checks$checks["FixAlleleIDs"]){
    vmsg("MADC file has not been processed by HapApp.", verbose = verbose, level = 1)
    if(get_REF_ALT){
      if(!mi_has_ref_alt) stop("MADC file has not been processed by HapApp. BIGr only provide results if get_REF_ALT is set to FALSE or if is TRUE but a marker_info with REF and ALT information is provided.")
    }
    # The check points to FALSE if the 6 initial rows exist or if there are no fixed allele ID (aka _0001, _0002, etc)
    n <- nrow(report)
    idx <- seq_len(min(6L, n))
    first_col_vals <- report[[1]][idx]
    all_blank_or_star <- all(first_col_vals %in% c("", "*"), na.rm = TRUE)
    # Also require that both _0001 and _0002 appear in AlleleID
    if(all_blank_or_star) {
      colnames(report) <- report[7,]
      report <- report[-c(1:7),]
    }
  }

  if(checks$checks["allNArow"]){
    idx <- apply(report, 1, function(x) all(is.na(x) | x == ""))
    report <- report[!idx, ]
    vmsg("MADC contains rows with all NA values. Rows %s will be removed.", verbose = verbose, level = 1, type = ">>", paste(which(idx), collapse = ", "))
  }

  if(checks$checks["allNAcol"]){
    idx <- apply(report, 2, function(x) all(is.na(x) | x == ""))
    report <- report[, !idx]
    vmsg("MADC contains columns with all NA values. Columns %s will be removed.", verbose = verbose, level = 1,  type = ">>", paste0(which(idx), collapse = ","))
  }

  if(!isTRUE(checks$checks["ChromPos"])) {
    if(is.null(markers_info)){
      stop("CloneID column does not follow the 'Chr_Pos'. ",
           "Please provide a markers_info file with at least 'CloneID'/'BI_markerID', ",
           "'Chr', and 'Pos' columns.")
    } else {

      if(!all(c("Chr", "Pos") %in% colnames(mi_df)))
        stop("CloneID column does not follow the 'Chr_Pos' format. ",
             "markers_info must contain at least 'Chr' and 'Pos' columns to remap marker IDs.")

    }
  }

  if(get_REF_ALT) {

    if(mi_has_ref_alt) {
      # markers_info supplies REF and ALT — no botloci required
      vmsg("markers_info contains Ref and Alt columns. REF and ALT will be taken from markers_info.",
           verbose = verbose, level = 1, type = ">>")

    } else {
      if(checks$checks["Indels"])
        stop("Indels detected in MADC file. Since get_REF_ALT = TRUE, a markers_info file with REF/ALT information is required.")

      # REF/ALT must be extracted from probe sequences — botloci is required
      if(is.null(botloci_file) || (!file.exists(botloci_file) && !url_exists(botloci_file)))
        stop("get_REF_ALT = TRUE but no markers_info file with Ref and Alt columns was provided neither a botloci_file to extrat REF/ALT from probe sequences. Please provide one of the these files or set get_REF_ALT to FALSE.")

      # Validate that CloneIDs match botloci marker names (after any remapping)
      botloci         <- read.table(botloci_file, header = FALSE)
      checked_botloci <- check_botloci(botloci, report, ChromPos = checks$checks["ChromPos"], mi_df = mi_df, verbose = verbose)
      botloci         <- checked_botloci[[1]]
      report          <- checked_botloci[[2]]

    }
  }

  # Throw message if OtherAlleles are present
  if(checks$checks["OtherAlleles"]) {
    vmsg("AlleleID contains alleles other than Ref and Alt. These will be ignored in the VCF output. Use function madc2vcf_all to include them.", verbose = verbose, level = 1, type = ">>")
  }

  # Make sure counts are numeric
  count_cols <- setdiff(colnames(report), c("CloneID", "AlleleID", "AlleleSequence"))
  report[count_cols] <- lapply(report[count_cols], function(x) as.numeric(as.character(x)))

  vmsg("Input checks done!", verbose = verbose, level = 1, type = ">>")

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

  if(get_REF_ALT && mi_has_ref_alt) {
    vmsg("Using markers_info for CHROM, POS, REF and ALT.", verbose = verbose, level = 0, type = ">>")

    if(is.null(mi_df)) mi_df <- read.csv(markers_info)
    id_col <- if ("BI_markerID" %in% colnames(mi_df)) "BI_markerID" else
              if ("CloneID"     %in% colnames(mi_df)) "CloneID"     else
      stop("The markers_info file must contain a marker ID column named either 'CloneID' or 'BI_markerID'.")

    if(checks$checks["Indels"])
      vmsg("Indels detected in MADC file. But it is okay because Ref and Alt are provided in markers_info.",
           verbose = verbose, level = 1, type = ">>")

    if(!all(c(id_col, "Chr", "Pos", "Ref", "Alt") %in% colnames(mi_df)))
      stop(paste0("The markers_info dataframe must contain the following columns: ",
                  id_col, ", Chr, Pos, Ref, Alt"))

    if(!all(rownames(ad_df) %in% mi_df[[id_col]])) {
      miss_CloneIDs <- rownames(ad_df)[!rownames(ad_df) %in% mi_df[[id_col]]]
      if(length(miss_CloneIDs) == nrow(ad_df)) stop("None of the MADC CloneID could be found in the markers_info CloneID or BI_markerID. Please make sure they match.")
      vmsg(paste("Not all MADC CloneID was found in the markers_info file. These markers will be removed:",
                 paste(miss_CloneIDs, collapse = " ")), verbose = verbose, level = 2, type = ">>")
      warning("Not all MADC CloneID was found in the markers_info file. These markers will be removed.")
    }
    matched <- mi_df[match(rownames(ad_df), mi_df[[id_col]]), ]

    new_df <- data.frame(CHROM = matched$Chr, POS = matched$Pos)
    new_df$TotalRef  <- rowSums(ref_df)
    new_df$TotalAlt  <- rowSums(alt_df)
    new_df$TotalSize <- rowSums(size_df)

    ref_base <- matched$Ref
    alt_base <- matched$Alt

  } else if(!is.null(markers_info) && !get_REF_ALT) {
    vmsg("markers_info file provided. Using CHROM and POS from the file.", verbose = verbose, level = 0, type = ">>")

    if(is.null(mi_df)) mi_df <- read.csv(markers_info)
    id_col <- if ("BI_markerID" %in% colnames(mi_df)) "BI_markerID" else
              if ("CloneID"     %in% colnames(mi_df)) "CloneID"     else
      stop("The markers_info file must contain a marker ID column named either 'CloneID' or 'BI_markerID'.")

    if(checks$checks["Indels"])
      vmsg("Indels detected in MADC file. Since get_REF_ALT = FALSE, Type and Indel_pos are not required in markers_info.",
           verbose = verbose, level = 1, type = ">>")

    if(!all(c(id_col, "Chr", "Pos") %in% colnames(mi_df)))
      stop(paste0("The markers_info dataframe must contain the following columns: ", id_col, ", Chr, Pos"))

    if(!all(rownames(ad_df) %in% mi_df[[id_col]])) {
      miss_CloneIDs <- rownames(ad_df)[!rownames(ad_df) %in% mi_df[[id_col]]]
      vmsg(paste("Not all MADC CloneID was found in the markers_info file. These markers will be removed:",
                 paste(miss_CloneIDs, collapse = " ")), verbose = verbose, level = 2, type = ">>")
      warning("Not all MADC CloneID was found in the markers_info file. These markers will be removed.")
    }
    matched <- mi_df[match(rownames(ad_df), mi_df[[id_col]]), ]

    new_df <- data.frame(CHROM = matched$Chr, POS = matched$Pos)
    new_df$TotalRef  <- rowSums(ref_df)
    new_df$TotalAlt  <- rowSums(alt_df)
    new_df$TotalSize <- rowSums(size_df)

    ref_base <- "."
    alt_base <- "."
    vmsg("REF and ALT not recovered (get_REF_ALT = FALSE).", verbose = verbose, level = 1, type = ">>")

  } else {
    vmsg(ifelse(get_REF_ALT,
                "Recovering CHROM and POS from CloneID for probe-sequence REF/ALT extraction...",
                "No markers_info file provided. Recovering CHROM and POS from CloneID..."),
         verbose = verbose, level = 0, type = ">>")

    # Split on the last underscore to handle chromosome names containing underscores
    # (e.g. Chr_01_000123456). When ChromPos was FALSE, check_botloci already
    # remapped CloneIDs to Chr_PaddedPos, so this split is always valid.
    new_df <- size_df %>%
      rownames_to_column(var = "row_name") %>%
      separate(row_name, into = c("CHROM", "POS"), sep = "_(?=[^_]*$)") %>%
      select(CHROM, POS)
    new_df$POS <- sub("^0+", "", new_df$POS)
    vmsg("CHROM and POS recovered from CloneID.", verbose = verbose, level = 1, type = ">>")

    new_df$TotalRef  <- rowSums(ref_df)
    new_df$TotalAlt  <- rowSums(alt_df)
    new_df$TotalSize <- rowSums(size_df)

    if(get_REF_ALT) {
      vmsg("get_REF_ALT = TRUE. Attempting to recover REF and ALT bases from probe sequences...",
           verbose = verbose, level = 0, type = ">>")

      csv <- get_counts(madc_object = report, collapse_matches_counts = collapse_matches_counts, verbose = FALSE)
      csv <- csv[which(csv$CloneID %in% rownames(ad_df)), ]

      ref_ord      <- csv$CloneID[grep("\\|Ref.*", csv$AlleleID)]
      alt_ord      <- csv$CloneID[grep("\\|Alt.*", csv$AlleleID)]
      orig_ref_seq <- csv$AlleleSequence[grep("\\|Ref.*", csv$AlleleID)]
      orig_alt_seq <- csv$AlleleSequence[grep("\\|Alt.*", csv$AlleleID)]

      if(all(sort(ref_ord) == sort(alt_ord))) {
        # Key sequences by CloneID, then reorder to MADC row order so that
        # loop index i always corresponds to rownames(size_df)[i].
        ref_seq_by_id <- setNames(orig_ref_seq, ref_ord)
        alt_seq_by_id <- setNames(orig_alt_seq, alt_ord)
        madc_ids      <- rownames(size_df)
        orig_ref_seq  <- ref_seq_by_id[madc_ids]
        orig_alt_seq  <- alt_seq_by_id[madc_ids]

        more_poly <- no_diff <- 0
        ref_base <- alt_base <- rep(NA_character_, length(madc_ids))
        names(ref_base) <- names(alt_base) <- madc_ids
        for(i in seq_along(madc_ids)) {
          if(is.na(orig_ref_seq[i]) || is.na(orig_alt_seq[i])) next
          temp_list <- strsplit(c(orig_ref_seq[i], orig_alt_seq[i]), "")
          if(length(temp_list[[1]]) != length(temp_list[[2]]))
            stop(paste0("Marker '", madc_ids[i], "' has Ref and Alt probe sequences of different lengths ",
                        "(", length(temp_list[[1]]), " vs ", length(temp_list[[2]]), "). ",
                        "This should not happen for SNP markers. ",
                        "If this is an indel, please provide a markers_info file with Ref and Alt columns."))
          idx_diff  <- which(temp_list[[1]] != temp_list[[2]])

          if(length(idx_diff) > 1) {
            ref_base[i] <- NA
            alt_base[i] <- NA
            more_poly <- more_poly + 1
          } else if(length(idx_diff) == 1) {
            orig_ref_base <- temp_list[[1]][idx_diff]
            orig_alt_base <- temp_list[[2]][idx_diff]
            if(madc_ids[i] %in% botloci[, 1]) {
              ref_base[i] <- as.character(reverseComplement(DNAString(orig_ref_base)))
              alt_base[i] <- as.character(reverseComplement(DNAString(orig_alt_base)))
            } else {
              ref_base[i] <- orig_ref_base
              alt_base[i] <- orig_alt_base
            }
          } else {
            ref_base[i] <- NA
            alt_base[i] <- NA
            no_diff <- no_diff + 1
          }
        }
        if(more_poly > 0) vmsg(paste(more_poly, "markers removed: more than one polymorphism between Ref and Alt sequences."), verbose = verbose, level = 2, type = ">>")
        if(no_diff   > 0) vmsg(paste(no_diff,   "markers removed: no differences found between Ref and Alt sequences."),       verbose = verbose, level = 2, type = ">>")

      } else {
        ref_base <- "."
        alt_base <- "."
        vmsg("REF and ALT bases could not be recovered: missing reference or alternative sequences.",
             verbose = verbose, level = 1, type = ">>")
      }

    } else {
      # ── get_REF_ALT = FALSE, no markers_info ─────────────────────────
      ref_base <- "."
      alt_base <- "."
      vmsg("REF and ALT not recovered (get_REF_ALT = FALSE).", verbose = verbose, level = 1, type = ">>")
    }
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
  m_size <- melt(as.matrix(size_df), varnames = c("Row", "Col"), value.name = "Value")
  m_ref  <- melt(as.matrix(ref_df),  varnames = c("Row", "Col"), value.name = "Value")
  m_ad   <- melt(as.matrix(ad_df),   varnames = c("Row", "Col"), value.name = "Value")

  combined_long <- m_size
  combined_long$Combined <- paste(m_size$Value, m_ref$Value, m_ad$Value, sep = ":")

  combined_wide <- suppressMessages(dcast(combined_long, Row ~ Col, value.var = "Combined"))
  rownames(combined_wide) <- combined_wide$Row
  combined_wide$Row <- NULL
  colnames(combined_wide) <- colnames(size_df)

  geno_df <- as.matrix(combined_wide)
  vmsg("Sample columns formatted for VCF", verbose = verbose, level = 1, type = ">>")

  #Combine the dataframes together
  vcf_df <- cbind(vcf_df,geno_df)
  vmsg("VCF dataframe prepared", verbose = verbose, level = 1, type = ">>")

  # Add # to the CHROM column name
  colnames(vcf_df)[1] <- "#CHROM"

  # Sort
  vcf_df <- vcf_df[order(vcf_df[,1],as.numeric(as.character(vcf_df[,2]))),]

  # Remove markers with NA CHROM/POS (unmatched in markers_info, Case 3)
  na_coord <- is.na(vcf_df[, 1]) | is.na(vcf_df$POS)
  if(any(na_coord)) {
    vmsg(paste(sum(na_coord), "markers removed: no matching entry found in markers_info."), verbose = verbose, level = 1, type = ">>")
    warning(paste(sum(na_coord), "markers removed: no matching entry found in markers_info."))
    vcf_df <- vcf_df[!na_coord, ]
  }

  if(sum(is.na(vcf_df$REF)) > 0) {
    vmsg(
      paste(
        sum(is.na(vcf_df$REF)),
        "markers removed because REF could not be reliably determined (e.g., multiple polymorphisms or no difference between probe sequences)."
      ),
      verbose = verbose,
      level = 1,
      type = ">>"
    )
    warning(
      paste(
        "Markers removed because REF could not be reliably determined (e.g., multiple polymorphisms or no difference between probe sequences):",
        sum(is.na(vcf_df$REF))
      )
    )
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
