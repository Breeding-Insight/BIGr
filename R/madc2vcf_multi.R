#' Convert MADC file to VCF using polyRAD for multiallelic genotyping
#'
#' This function converts a DArTag fixed allele ID MADC file to a VCF
#' containing multiallelic markers based on the microhaplotypes using
#' the polyRAD package's `readDArTag`, `IterateHWE` population model
#' and `RADdata2VCF` pipeline.
#'
#' @param madc_file character. Path or URL to the input MADC CSV file.
#' @param botloci_file character. Path or URL to the botloci file listing target
#'   IDs designed on the bottom strand.
#' @param outfile character. Path for the output VCF file.
#' @param markers_info character or NULL. Optional path or URL to a CSV file
#'   with marker metadata. Required when CloneIDs do not follow the
#'   \code{Chr_Pos} format; must contain \code{CloneID} (or
#'   \code{BI_markerID}), \code{Chr}, and \code{Pos} columns.
#' @param ploidy integer. Ploidy level of the samples passed to \code{taxaPloidy}.
#'   Default is 2.
#' @param vcf_header character or NULL. Optional path or URL to a VCF header file
#'   with the \code{##reference} and \code{##contig} lines of the reference
#'   genome, such as the VCF header files provided with the microhaplotype
#'   databases. Its \code{##contig} lines replace the ones written by polyRAD,
#'   which have no contig lengths. The other lines of the file are not used. A
#'   warning is given if a chromosome in the output is not listed in the
#'   \code{##contig} lines.
#' @param reference character or NULL. Optional reference genome ID written to
#'   the \code{##reference} line (e.g. \code{"GCF_020283565.1"}). Use it when no
#'   \code{vcf_header} file is available, or to replace the \code{##reference}
#'   line of \code{vcf_header}. If neither \code{vcf_header} nor
#'   \code{reference} is given, a warning is given and the header has no
#'   \code{##reference} line and only the \code{##contig} lines written by polyRAD.
#' @param align_indels logical. Whether to align the alleles of markers with
#'   indels before loading them into polyRAD (see "Markers with indels" below).
#'   If \code{FALSE}, polyRAD cuts the alleles of each marker to the length of
#'   its shortest allele, as in BIGr < 0.11.0. Default is \code{TRUE}.
#' @param asSNPs logical. Passed to \code{polyRAD::RADdata2VCF}. If
#'   \code{FALSE} (default), each VCF record is a marker and its alleles are
#'   the full haplotype sequences. If \code{TRUE}, each record is a variable
#'   site (SNP or indel) of a marker, with the read counts and genotypes of the
#'   haplotypes summed per site allele. Sites with more than two alleles are
#'   kept as multiallelic records.
#' @param verbose logical. Whether to print progress messages. Default is TRUE.
#'
#' @return Invisible NULL. Writes a VCF file to \code{outfile}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reads the MADC file and runs \code{check_madc_sanity}.
#'   \item Validates the botloci file against MADC CloneIDs using
#'         \code{check_botloci}, fixing any padding mismatches automatically.
#'   \item Converts lowercase bases in AlleleSequence to uppercase if detected.
#'   \item Removes all-NA rows and columns if detected.
#'   \item Aligns the alleles of markers with indels (see below), unless
#'         \code{align_indels = FALSE}.
#'   \item Writes the corrected data to a temporary file and passes it to
#'         \code{polyRAD::readDArTag}, with the bottom strand marker IDs from
#'         the botloci file. polyRAD reverse complements the sequences of these
#'         markers, so all VCF sequences are on the top strand of the reference.
#'   \item Sets the position of each marker to the start of its amplicon: the
#'         target SNP position (from the CloneID or \code{markers_info}) minus
#'         the position of the target SNP in the Ref sequence, plus 1.
#'   \item Estimates overdispersion with \code{polyRAD::TestOverdispersion} and
#'         calls \code{polyRAD::IterateHWE}, then exports the result with
#'         \code{polyRAD::RADdata2VCF}, as haplotypes or, with
#'         \code{asSNPs = TRUE}, as SNPs.
#'   \item Adds back the MADC information that polyRAD does not keep (see below).
#'   \item Adds the reference genome lines from \code{vcf_header} and/or
#'         \code{reference} and a \code{##BIGrCommandLine} line to the VCF header.
#' }
#'
#' **Information added to the polyRAD VCF**
#'
#' \code{polyRAD::RADdata2VCF} writes the full amplicon sequences as alleles,
#' sets the positions to the start of the amplicons and does not keep the MADC
#' IDs. The marker of each record is taken from the INFO \code{LU} field written
#' by polyRAD (index of the marker in the RADdata object), and its REF and ALT
#' sequences (or their reverse complements) are matched to the MADC
#' \code{AlleleSequence} of that marker only, to add:
#' \itemize{
#'   \item ID column: the target SNP ID (marker CloneID).
#'   \item INFO \code{targetSNP}: the same target SNP ID.
#'   \item INFO \code{DArTstrand}: \code{-} for markers in the botloci file,
#'         \code{+} otherwise.
#'   \item INFO \code{hapID}: AlleleIDs of the REF and ALT alleles, in order,
#'         separated by commas (\code{.} if a sequence is not found).
#' }
#' With \code{asSNPs = TRUE}, a marker gives one record per variable site, so
#' only the record at the target SNP position (or, if there is none, such as
#' for an indel written from the base before it, the first record whose REF
#' spans it) gets the target SNP ID in the ID column. The other records get a
#' \code{Chr_Pos} ID, with the position zero-padded like the target SNP ID.
#' \code{targetSNP} and \code{DArTstrand} are added to all records, and
#' \code{hapID} is not added. The genotypes are called on the haplotypes as
#' with \code{asSNPs = FALSE}; only for writing the SNP records, the bases
#' missing at the start or end of a shorter allele are set to \code{N}, so they
#' are not written as deletions. At a SNP that a shorter allele does not reach,
#' that allele is written as an \code{N} ALT allele. An indel at the first
#' base of an amplicon is written with \code{N} as the base before it, and
#' \code{POS} 1 bp before the amplicon start.
#'
#' The \code{##SAMPLE} header lines written by polyRAD are removed. A warning
#' is given if REF sequences are not found among the MADC alleles of their marker.
#'
#' **Markers with indels**
#'
#' polyRAD needs all alleles of a marker to have the same length, and
#' otherwise cuts them to the length of the shortest allele, losing the end of
#' the longer alleles. For markers whose alleles differ in length, each allele
#' is first aligned to the marker's Ref allele (global alignment with
#' \code{pwalign::pairwiseAlignment}) and the gaps are filled with \code{-}, so
#' all alleles have the same length and no bases are lost: \code{-} for a
#' deletion, and \code{.} in the columns of an insertion for the alleles
#' without it, including Ref. polyRAD removes them when writing the VCF, so the
#' VCF has the full allele sequences, and does not count \code{.} in the SNP
#' positions of \code{asSNPs = TRUE}. As polyRAD counts gaps as differences, a
#' shorter allele whose bases all match a longer allele stays a separate allele.
#'
#' With \code{align_indels = FALSE}, polyRAD cuts the alleles instead. Alleles
#' that differ only after the cut are merged, and if the target SNP is cut
#' away, Ref and Alt become the same allele: a warning gives the number of
#' these markers, which are still placed using the MADC Ref and Alt sequences.
#'
#' **Sanity check behaviour and requirements**
#'
#' The function always stops if IUPAC codes, unpaired Ref/Alt sequences, or
#' unfixed AlleleIDs are detected (see \code{check_madc_sanity}). For the
#' remaining checks the required inputs are:
#'
#' | Check | Status | Required |
#' |---|---|---|
#' | **Indels** | detected | `botloci_file` |
#' | | not detected | `botloci_file` |
#' | **ChromPos** | valid | `botloci_file` |
#' | | invalid | `markers_info` with `Chr`/`Pos` + `botloci_file` |
#'
#' @importFrom utils read.csv write.csv read.table
#'
#' @export
madc2vcf_multi <- function(madc_file,
                           botloci_file,
                           outfile,
                           markers_info = NULL,
                           ploidy = 2L,
                           vcf_header = NULL,
                           reference = NULL,
                           align_indels = TRUE,
                           asSNPs = FALSE,
                           verbose = TRUE) {

  vmsg("Running BIGr madc2vcf_multi", verbose = verbose, level = 0, type = ">>")
  vmsg("madc_file    : %s", verbose = verbose, level = 1, type = ">>", madc_file)
  vmsg("botloci_file : %s", verbose = verbose, level = 1, type = ">>", botloci_file)
  vmsg("markers_info : %s", verbose = verbose, level = 1, type = ">>", if (is.null(markers_info)) "NULL" else markers_info)
  vmsg("outfile      : %s", verbose = verbose, level = 1, type = ">>", outfile)
  vmsg("ploidy       : %s", verbose = verbose, level = 1, type = ">>", ploidy)
  vmsg("vcf_header   : %s", verbose = verbose, level = 1, type = ">>", if (is.null(vcf_header)) "NULL" else vcf_header)
  vmsg("reference    : %s", verbose = verbose, level = 1, type = ">>", if (is.null(reference)) "NULL" else reference)
  vmsg("align_indels : %s", verbose = verbose, level = 1, type = ">>", align_indels)
  vmsg("asSNPs       : %s", verbose = verbose, level = 1, type = ">>", asSNPs)

  vmsg("Checking inputs", verbose = verbose, level = 0, type = ">>")

  if (!(file.exists(madc_file) | url_exists(madc_file)))       stop("MADC file not found. Please provide a valid path or URL.")
  if (!(file.exists(botloci_file) | url_exists(botloci_file))) stop("Botloci file not found. Please provide a valid path or URL.")
  if (!is.null(markers_info) && !(file.exists(markers_info) | url_exists(markers_info))) stop("markers_info file not found. Please provide a valid path or URL.")
  if (!is.numeric(ploidy) || ploidy < 1) stop("ploidy must be a positive integer.")
  if (!is.logical(align_indels) || length(align_indels) != 1L || is.na(align_indels))
    stop("align_indels must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(asSNPs) || length(asSNPs) != 1L || is.na(asSNPs))
    stop("asSNPs must be TRUE or FALSE.", call. = FALSE)
  # Create a VCF header line with metadata about the command and its parameters
  bigr_meta <- paste0('##BIGrCommandLine.madc2vcf_multi=<ID=madc2vcf_multi,Version="',
                      packageVersion("BIGr"), '",Data="',
                      Sys.time(),'",CommandLine="', vcf_quote_escape(paste0('> madc2vcf_multi(',deparse1(substitute(madc_file)),', ',
                      "botloci_file= ", botloci_file, ', ',
                      "outfile= ", outfile, ', ',
                      "markers_info= ", markers_info, ', ',
                      "ploidy= ", ploidy, ', ',
                      "vcf_header= ", vcf_header, ', ',
                      "reference= ", reference, ', ',
                      "align_indels= ", align_indels, ', ',
                      "asSNPs= ", asSNPs, ', ',
                      "verbose= ", verbose,')')), '">')

  # Reference genome lines for the VCF header
  ref_lines <- vcf_reference_lines(vcf_header, reference)

  # ---- Load markers_info if provided ----
  mi_df <- if (!is.null(markers_info)) read.csv(markers_info) else NULL

  # ---- Read and sanity-check ----
  report <- read.csv(madc_file, check.names = FALSE)
  checks <- check_madc_sanity(report)

  messages_results <- mapply(function(check, message) {
    if (isTRUE(check)) message[1] else message[2]
  }, checks$checks, checks$messages)

  for (i in seq_along(messages_results))
    vmsg(messages_results[i], verbose = verbose, level = 1, type = ">>")

  if (!checks$checks["Columns"])
    stop("The MADC file is missing required columns (CloneID, AlleleID, AlleleSequence)")

  if (checks$checks["IUPACcodes"])
    stop("MADC Allele Sequences contain different IUPAC (non-ATCG) codes in REF and ALT sequences. Please run HapApp to clean MADC file before using this function.")

  if (!isTRUE(checks$checks["RefAltSeqs"]))
    stop("Not all Ref sequences have a corresponding Alt or vice versa. Please provide a complete MADC file before using this function.")

  if (!isTRUE(checks$checks["FixAlleleIDs"]))
    stop("The MADC file does not have fixed AlleleIDs. Please process the MADC file through HapApp before using this function.")

  if (!isTRUE(checks$checks["ChromPos"])) {
    if (is.null(markers_info))
      stop("CloneID column does not follow the 'Chr_Pos' format. ",
           "Please provide a markers_info file with at least 'CloneID'/'Marker_ID'/'BI_markerID', ",
           "'Chr', and 'Pos' columns.")
    if (!all(c("Chr", "Pos") %in% colnames(mi_df)))
      stop("CloneID column does not follow the 'Chr_Pos' format. ",
           "markers_info must contain at least 'Chr' and 'Pos' columns to remap marker IDs.")
  }

  # ---- Check botloci vs MADC CloneIDs ----
  vmsg("Checking botloci file", verbose = verbose, level = 0, type = ">>")
  cloneids_before <- report$CloneID
  botloci_df      <- read.table(botloci_file, header = FALSE)
  checked_botloci <- check_botloci(botloci_df, report, ChromPos = checks$checks["ChromPos"], mi_df = mi_df, verbose = verbose)
  botloci_df      <- checked_botloci[[1]]
  report          <- checked_botloci[[2]]
  mi_df           <- checked_botloci[[3]]
  cloneid_changed  <- !identical(report$CloneID, cloneids_before)
  # polyRAD::readDArTag() takes the bottom strand marker IDs as a character
  # vector, not a file path (a path matches no marker, so all are set to top)
  botloci_ids <- unique(trimws(as.character(botloci_df$V1)))
  vmsg("%s of %s MADC markers are on the bottom strand (botloci)", verbose = verbose, level = 1, type = ">>",
       sum(unique(report$CloneID) %in% botloci_ids), length(unique(report$CloneID)))

  # ---- Corrections: only create a temp file if needed ----
  need_temp <- isTRUE(checks$checks["allNArow"]) || isTRUE(checks$checks["allNAcol"]) || cloneid_changed

  # Sample names with special characters: cleaned and passed to polyRAD in a temp file
  if (isTRUE(checks$checks["SpecialChar"])) {
    report <- clean_sample_names(report, verbose = verbose)
    need_temp <- TRUE
  }

  # polyRAD requires integer read counts. Counts stored as e.g. "12.0" are read
  # as double, so coerce whole-number count columns to integer and pass a temp file.
  count_cols <- setdiff(colnames(report), c("AlleleID", "CloneID", "AlleleSequence", "readCountSum"))
  non_int <- count_cols[!vapply(report[count_cols], is.integer, logical(1))]
  if (length(non_int) > 0) {
    vals <- unlist(report[non_int], use.names = FALSE)
    if (!is.numeric(vals) || any(vals %% 1 != 0, na.rm = TRUE))
      stop("MADC read count columns must contain whole numbers. Non-integer values found in: ",
           paste(utils::head(non_int, 5), collapse = ", "),
           if (length(non_int) > 5) ", ..." else "")
    vmsg("Converting %s read count column(s) stored as decimals to integer",
         verbose = verbose, level = 1, type = ">>", length(non_int))
    report[non_int] <- lapply(report[non_int], as.integer)
    need_temp <- TRUE
  }

  # MADC sequences before the alignment, used to place the target SNP of
  # markers whose Ref and Alt become identical in polyRAD (align_indels = FALSE)
  madc_seqs <- stats::setNames(toupper(report$AlleleSequence), report$AlleleID)

  # polyRAD needs all alleles of a marker to have the same length, and cuts
  # them to the shortest allele otherwise. Align the alleles of markers with
  # indels to their Ref allele and fill the gaps instead, so no bases are lost.
  # polyRAD removes the gaps when writing the VCF.
  if (align_indels) {
    aligned <- align_madc_indels(clone_id = report$CloneID, allele_id = report$AlleleID, seqs = report$AlleleSequence)
    if (aligned$n_markers > 0) {
      vmsg("Aligning the alleles of %s marker(s) with indels to their Ref allele (gaps filled with '-' and '.')",
           verbose = verbose, level = 1, type = ">>", aligned$n_markers)
      report$AlleleSequence <- aligned$seqs
      need_temp <- TRUE
    }
  } else {
    n_cut <- sum(tapply(nchar(report$AlleleSequence), report$CloneID, function(x) length(unique(x)) > 1))
    if (n_cut > 0)
      vmsg("align_indels = FALSE: polyRAD will cut the alleles of %s marker(s) to the length of their shortest allele",
           verbose = verbose, level = 1, type = ">>", n_cut)
  }

  if (need_temp) {
    if (checks$checks["LowerCase"]) {
      vmsg("MADC Allele Sequences contain lowercase characters. Converting to uppercase",
           verbose = verbose, level = 1, type = ">>")
      report$AlleleSequence <- toupper(report$AlleleSequence)
    }

    if (checks$checks["allNArow"]) {
      idx <- apply(report, 1, function(x) all(is.na(x) | x == ""))
      vmsg("Removing %s all-NA row(s)", verbose = verbose, level = 1, type = ">>", sum(idx))
      report <- report[!idx, ]
    }

    if (checks$checks["allNAcol"]) {
      idx <- apply(report, 2, function(x) all(is.na(x) | x == ""))
      vmsg("Removing %s all-NA column(s)", verbose = verbose, level = 1, type = ">>", sum(idx))
      report <- report[, !idx]
    }

    tmp_madc <- tempfile(fileext = ".csv")
    on.exit(unlink(tmp_madc), add = TRUE)
    write.csv(report, tmp_madc, row.names = FALSE, quote = TRUE)
    input_file <- tmp_madc
  } else {
    if (checks$checks["LowerCase"])
      vmsg("MADC Allele Sequences contain lowercase characters. polyRAD will handle them",
           verbose = verbose, level = 1, type = ">>")
    input_file <- madc_file
  }

  vmsg("Loading MADC into polyRAD", verbose = verbose, level = 0, type = ">>")

  if (!requireNamespace("polyRAD", quietly = TRUE)) {
    stop("Package 'polyRAD' is required for madc2vcf_multi(). Please install it with install.packages('polyRAD').", call. = FALSE)
  }
  raddat <- polyRAD::readDArTag(
    file              = input_file,
    botloci           = botloci_ids,
    n.header.rows     = 0L,
    sample.name.row   = 1,
    trim.sample.names = "",
    taxaPloidy        = as.integer(ploidy)
  )

  # polyRAD derives Chr/Pos by splitting the CloneID at its last "_". When
  # CloneIDs do not follow Chr_Pos (and were not remapped by check_botloci),
  # take Chr/Pos from markers_info instead so RADdata2VCF gets valid positions.
  loc_ids <- rownames(raddat$locTable)
  target_pos <- suppressWarnings(as.integer(sub("^.+_", "", loc_ids)))
  if (!isTRUE(checks$checks["ChromPos"]) && !is.null(mi_df) &&
      (anyNA(raddat$locTable$Pos) || !cloneid_changed)) {
    use_col <- pick_markers_info_id_col(mi_df, loc_ids)
    idx <- match(loc_ids, mi_df[[use_col]])
    found <- !is.na(idx)
    vmsg("Using Chr and Pos from markers_info (%s column) for %s of %s loci",
         verbose = verbose, level = 1, type = ">>", use_col, sum(found), length(loc_ids))
    raddat$locTable$Chr[found] <- as.character(mi_df$Chr[idx[found]])
    target_pos[found] <- as.integer(mi_df$Pos[idx[found]])
  }

  # POS of each record is the amplicon start: the target SNP position minus the
  # position of the target SNP in the Ref sequence, counting only bases (not the
  # gaps added by align_madc_indels()). polyRAD counts the gaps too, so its Pos
  # is replaced.
  snp_in_ref <- target_index_in_ref(raddat)
  # Without the alignment, polyRAD may cut the target SNP away, leaving Ref and
  # Alt identical. Place these markers with the MADC Ref and Alt sequences.
  lost <- is.na(snp_in_ref) & !is.na(target_pos)
  if (any(lost)) {
    snp_in_ref[lost] <- target_index_from_madc(
      clone_id = report$CloneID, allele_id = report$AlleleID, seqs = unname(madc_seqs[report$AlleleID]),
      loci = loc_ids[lost], bottom = loc_ids[lost] %in% botloci_ids,
      ref_len = nchar(gsub("[-.]", "", as.character(raddat$locTable$Ref[lost]))))
    n_lost <- sum(!is.na(snp_in_ref[lost]))
    if (n_lost > 0)
      warning(n_lost, " marker(s) have identical Ref and Alt alleles in polyRAD, because their target SNP ",
              "is outside the length of their shortest allele",
              if (!align_indels) " (align_indels = FALSE)", ". Their VCF records do not include the target SNP.",
              call. = FALSE)
  }
  ok <- !is.na(snp_in_ref) & !is.na(target_pos)
  raddat$locTable$Pos[ok] <- target_pos[ok] - snp_in_ref[ok] + 1L
  raddat$locTable$Pos[is.na(target_pos)] <- NA_integer_
  bad_pos <- loc_ids[is.na(raddat$locTable$Pos)]
  if (length(bad_pos) > 0)
    stop("Could not determine the position of ", length(bad_pos), " marker(s) (e.g. ",
         paste(utils::head(bad_pos, 5), collapse = ", "),
         "). Make sure they are in markers_info with valid 'Chr' and 'Pos' values.")

  # polyRAD returns optimal = NA when the best value is at the edge of the
  # tested range, so widen the range until an interior optimum is found.
  to_test <- seq(6, 20, by = 2)  # polyRAD default
  my_ovdisp <- NA
  best <- NA
  for (attempt in 1:4) {
    res <- run_test_od(raddat, to_test, verbose = verbose)
    my_ovdisp <- res$optimal
    if (!is.na(my_ovdisp)) break
    best <- res$best
    if (is.na(best)) break
    if (best == min(to_test) && min(to_test) > 0.5) {
      to_test <- seq(0.5, min(to_test) + 1, by = 0.5)
    } else if (best == max(to_test)) {
      to_test <- seq(max(to_test) - 1, max(to_test) * 3, length.out = 20)
    } else {
      break
    }
  }
  if (is.na(my_ovdisp)) my_ovdisp <- best
  if (is.na(my_ovdisp))
    stop("Could not estimate the overdispersion parameter with polyRAD::TestOverdispersion().")

  vmsg("Running HWE iteration (overdispersion = %s)", verbose = verbose, level = 0, type = ">>", my_ovdisp)

  raddat_hwe <- polyRAD::IterateHWE(raddat, overdispersion = my_ovdisp)

  vmsg("Writing VCF to %s", verbose = verbose, level = 0, type = ">>", outfile)

  if (asSNPs) raddat_hwe <- anchor_leading_gaps(mask_terminal_gaps(raddat_hwe))
  polyRAD::RADdata2VCF(raddat_hwe, file = outfile, asSNPs = asSNPs, hindhe = FALSE)

  # polyRAD drops the MADC IDs. Add the target SNP ID, strand and (for
  # haplotype records) haplotype IDs back, and add the reference genome lines
  # and the BIGr command line to the header.
  vmsg("Adding target SNP IDs, strand%s to the VCF", verbose = verbose, level = 0, type = ">>",
       if (asSNPs) "" else " and haplotype IDs")
  hap_info <- madc_hap_lookup(report, botloci_ids, loci = rownames(raddat_hwe$locTable), cut = !align_indels)
  if (asSNPs) hap_info$target_pos <- stats::setNames(target_pos, loc_ids)
  meta <- update_vcf_reference_lines(read_vcf_meta(outfile), ref_lines)
  meta <- add_hap_info_header(meta, hapID = !asSNPs)
  n_missing <- 0L
  replace_vcf_meta(outfile, c(meta, bigr_meta), transform = function(lines) {
    res <- add_hap_info_records(lines, hap_info)
    n_missing <<- n_missing + res$n_missing
    res$lines
  })
  if (n_missing > 0)
    warning(n_missing, " VCF record(s) have a REF sequence not found among the MADC alleles of their marker. ",
            "Their REF hapID is left as '.'.", call. = FALSE)
  check_vcf_contigs(raddat_hwe$locTable$Chr, ref_lines$contigs)

  vmsg("Done!", verbose = verbose, level = 0, type = ">>")

  invisible(NULL)
}

#' Read the meta lines (starting with ##) of a VCF file
#'
#' @param file Path to an uncompressed VCF file.
#' @return Character vector with the meta lines.
#' @noRd
read_vcf_meta <- function(file) {
  con <- file(file, "r")
  on.exit(close(con))
  meta <- character(0)
  repeat {
    lines <- readLines(con, n = 1000L)
    if (length(lines) == 0) break
    is_meta <- startsWith(lines, "##")
    meta <- c(meta, lines[is_meta])
    if (!all(is_meta)) break
  }
  meta
}

#' Replace the meta lines (starting with ##) of a VCF file
#'
#' The #CHROM line and the records are copied unchanged, in chunks, so the file
#' is not read into memory at once.
#'
#' @param file Path to an uncompressed VCF file.
#' @param meta Character vector with the new meta lines.
#' @param transform Optional function to also edit the variant records (the
#'   lines after the #CHROM line) while the file is copied. The records are
#'   read in chunks of up to 10000 lines; for each chunk, `transform` receives
#'   a character vector with the record lines, unsplit (tab-separated), and
#'   must return the character vector of lines to write in their place. The
#'   #CHROM line is never passed to it. With the default NULL, the records are
#'   copied unchanged.
#' @return No return value, called for side effects.
#' @noRd
replace_vcf_meta <- function(file, meta, transform = NULL) {
  tmp <- tempfile(fileext = ".vcf")
  con_in <- file(file, "r")
  con_out <- file(tmp, "w")
  writeLines(meta, con_out)
  repeat {
    lines <- readLines(con_in, n = 10000L)
    if (length(lines) == 0) break
    lines <- lines[!startsWith(lines, "##")]
    if (!is.null(transform) && length(lines) > 0) {
      is_rec <- !startsWith(lines, "#")
      if (any(is_rec)) lines <- c(lines[!is_rec], transform(lines[is_rec]))
    }
    if (length(lines) > 0) writeLines(lines, con_out)
  }
  close(con_in)
  close(con_out)
  file.copy(tmp, file, overwrite = TRUE)
  unlink(tmp)
  invisible()
}

#' Build the lookup tables used to add the MADC information to the polyRAD VCF
#'
#' The MADC sequences and their reverse complements are mapped to the
#' AlleleID, within each marker (CloneID), so a sequence shared by two markers
#' (e.g. overlapping amplicons) is not assigned to the wrong marker.
#'
#' The gaps added by `align_madc_indels()` ("-" and ".") are removed. With
#' `cut = TRUE` (no alignment), the sequences of each marker are cut to the
#' length of its shortest allele, as `polyRAD::readDArTag()` does. Alleles that
#' become identical are merged by polyRAD; the last of their AlleleIDs (in MADC
#' order) is used for them.
#'
#' @param report MADC data.frame with fixed AlleleIDs (`CloneID`, `AlleleID`,
#'   `AlleleSequence`), after the CloneID corrections of `check_botloci()`.
#' @param botloci Character vector with the IDs of the bottom strand markers,
#'   matching the CloneIDs of `report`.
#' @param loci Character vector with the marker IDs in the order of the
#'   RADdata object (`rownames(locTable)`), used to read the marker of each VCF
#'   record from its INFO `LU` field. If NULL, the marker is found from the
#'   REF sequence.
#' @param cut Logical. Whether to cut the sequences of each marker to the
#'   length of its shortest allele.
#' @return A list with `seq_to_id` (AlleleIDs named by sequence),
#'   `locus_seq_to_id` (AlleleIDs named by "CloneID sequence"), `botloci` and
#'   `loci`.
#' @noRd
madc_hap_lookup <- function(report, botloci, loci = NULL, cut = FALSE) {
  seqs <- gsub("[-.]", "", toupper(report$AlleleSequence))  # gaps from align_madc_indels()
  if (cut) {
    min_len <- stats::ave(nchar(seqs), report$CloneID, FUN = min)
    seqs <- substr(seqs, 1, min_len)
  }
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs)))
  keys <- c(rbind(seqs, rc))
  ids <- rep(as.character(report$AlleleID), each = 2)
  clone <- rep(as.character(report$CloneID), each = 2)
  keep <- !duplicated(keys, fromLast = TRUE)
  locus_keys <- paste(clone, keys)
  keep_locus <- !duplicated(locus_keys, fromLast = TRUE)
  list(seq_to_id = stats::setNames(ids[keep], keys[keep]),
       locus_seq_to_id = stats::setNames(ids[keep_locus], locus_keys[keep_locus]),
       botloci = as.character(botloci),
       loci = if (is.null(loci)) NULL else as.character(loci))
}

#' Add the targetSNP, DArTstrand and hapID INFO lines to a VCF header
#'
#' The ##SAMPLE lines written by polyRAD are removed.
#'
#' @param meta Meta lines (starting with ##) of the VCF header.
#' @param hapID Logical. Whether to add the hapID INFO line (haplotype records).
#' @return Character vector with the updated meta lines.
#' @noRd
add_hap_info_header <- function(meta, hapID = TRUE) {
  meta <- meta[!startsWith(meta, "##SAMPLE=")]
  meta <- meta[!grepl("^##INFO=<ID=(targetSNP|DArTstrand|hapID),", meta)]
  new_info <- c(
    '##INFO=<ID=targetSNP,Number=1,Type=String,Description="Target SNP ID (MADC CloneID) of the DArTag amplicon">',
    '##INFO=<ID=DArTstrand,Number=1,Type=String,Description="DArTag amplicon strand">',
    '##INFO=<ID=hapID,Number=R,Type=String,Description="Haplotype IDs (MADC AlleleIDs) of REF and ALT alleles">'
  )
  if (!hapID) new_info <- new_info[-3]
  pos <- utils::tail(grep("^##INFO=", meta), 1)
  if (length(pos) == 0) pos <- utils::tail(grep("^##(fileformat|reference|contig)=", meta), 1)
  if (length(pos) == 0) pos <- 0
  append(meta, new_info, after = pos)
}

#' Add the MADC information to polyRAD VCF records
#'
#' The marker of each record is read from the INFO `LU` field written by
#' polyRAD (index of the marker in the RADdata object), or, when `hap_info$loci`
#' is NULL or `LU` is missing, from the CloneID of the MADC allele matching the
#' REF sequence. The ID column is set to this target SNP ID and the targetSNP,
#' DArTstrand (`-` if the marker is in the botloci file) and hapID INFO fields
#' are added. REF and ALT sequences are matched only to the AlleleIDs of their
#' own marker. POS, REF and ALT are not changed.
#'
#' When `hap_info$target_pos` is given (SNP records, `asSNPs = TRUE`), hapID is
#' not added and the ID column is the target SNP ID only for the record at the
#' target SNP position (or, if there is none, the first record whose REF spans
#' it); the other records get `Chr_Pos` IDs, with the position zero-padded to
#' the width used in the target SNP ID.
#'
#' @param lines Record lines of the VCF.
#' @param hap_info List returned by `madc_hap_lookup()`.
#' @return A list with `lines` (updated record lines) and `n_missing` (number
#'   of records whose REF sequence was not found in the MADC for their marker).
#' @noRd
add_hap_info_records <- function(lines, hap_info) {
  # Split only the 8 fixed columns; FORMAT and sample columns are kept as is
  fixed_end <- attr(regexpr("^([^\t]*\t){8}", lines), "match.length")
  if (any(fixed_end < 0)) stop("VCF records with less than 9 columns found.")
  fixed <- do.call(rbind, strsplit(substr(lines, 1, fixed_end - 1), "\t", fixed = TRUE))
  rest <- substring(lines, fixed_end + 1)

  ref <- toupper(fixed[, 4])

  # Marker of each record: from polyRAD's LU index, else from the REF sequence
  target <- rep(NA_character_, nrow(fixed))
  if (!is.null(hap_info$loci)) {
    lu <- suppressWarnings(as.integer(sub("^(.*;)?LU=([0-9]+)(;.*)?$", "\\2", fixed[, 8])))
    ok <- !is.na(lu) & lu >= 1 & lu <= length(hap_info$loci)
    target[ok] <- hap_info$loci[lu[ok]]
  }
  by_seq <- is.na(target)
  target[by_seq] <- sub("\\|.*$", "", unname(hap_info$seq_to_id[ref[by_seq]]))
  target[is.na(target)] <- "."
  bottom <- target %in% hap_info$botloci

  ref_id <- unname(hap_info$locus_seq_to_id[paste(target, ref)])
  found <- !is.na(ref_id)

  # ALT alleles: one row per allele, then the IDs are pasted back per record
  alts <- strsplit(toupper(fixed[, 5]), ",", fixed = TRUE)
  rec <- rep(seq_along(alts), lengths(alts))
  alt <- unlist(alts, use.names = FALSE)
  is_missing <- alt == "."
  alt_id <- unname(hap_info$locus_seq_to_id[paste(target[rec], alt)])
  alt_id[is.na(alt_id)] <- "."
  hap_alt <- tapply(alt_id[!is_missing], factor(rec[!is_missing], levels = seq_along(alts)),
                    paste, collapse = ",")

  hap_id <- ifelse(found, ref_id, ".")
  has_alt <- !is.na(hap_alt)
  hap_id[has_alt] <- paste(hap_id[has_alt], hap_alt[has_alt], sep = ",")

  has_target <- target != "."
  strand <- paste0("targetSNP=", target, ";DArTstrand=", ifelse(bottom, "-", "+"))
  if (!is.null(hap_info$target_pos)) {
    pos <- suppressWarnings(as.integer(fixed[, 2]))
    tpos <- unname(hap_info$target_pos[target])
    at_target <- has_target & !is.na(tpos) & !is.na(pos) & pos == tpos
    # An indel (written from the base before it) or a site merged into an
    # indel does not start at the target position: use the first record that
    # spans it, for markers without a record at the target position
    spans <- has_target & !is.na(tpos) & !is.na(pos) & pos < tpos & pos + nchar(fixed[, 4]) - 1L >= tpos &
      !target %in% target[at_target]
    at_target <- at_target | (spans & !duplicated(ifelse(spans, target, NA_character_), incomparables = NA))
    suffix <- sub("^.*_", "", target)
    width <- ifelse(grepl("^[0-9]+$", suffix), nchar(suffix), 0L)
    snp_id <- paste0(fixed[, 1], "_", sprintf("%0*d", width, pos))
    fixed[, 3] <- ifelse(at_target, target, snp_id)
    info <- strand
    found <- rep(TRUE, nrow(fixed))
  } else {
    fixed[has_target, 3] <- target[has_target]
    info <- paste0(strand, ";hapID=", hap_id)
  }
  fixed[, 8] <- ifelse(fixed[, 8] %in% c(".", ""), info, paste(fixed[, 8], info, sep = ";"))

  list(lines = paste(do.call(paste, c(as.data.frame(fixed), sep = "\t")), rest, sep = "\t"),
       n_missing = sum(!found))
}

#' Align the alleles of MADC markers with indels
#'
#' `polyRAD::readDArTag()` needs all alleles of a marker to have the same
#' length and otherwise cuts them to the length of the shortest allele, which
#' loses the end of the longer alleles. For each marker whose alleles differ in
#' length, every allele is aligned to the marker's Ref allele (global
#' alignment) and the pairwise alignments are combined, so all alleles get the
#' same columns:
#' \itemize{
#'   \item a deletion in an allele becomes "-" in that allele;
#'   \item an insertion in an allele adds columns that are "." in the Ref and in
#'         the alleles without that insertion. Different insertions at the same
#'         Ref position are left-aligned and padded with ".".
#' }
#' This is polyRAD's notation: "." is a column with no Ref base, so it is not
#' counted in the SNP positions of `RADdata2VCF(asSNPs = TRUE)`. polyRAD
#' treats "-" and "." as different from any base when comparing alleles and
#' removes them when writing the VCF, so the VCF has the full sequences. Markers whose alleles
#' all have the same length are not changed.
#'
#' @param clone_id,allele_id,seqs MADC CloneID, AlleleID and AlleleSequence
#'   columns.
#' @return A list with `seqs` (the sequences, aligned for markers with indels)
#'   and `n_markers` (number of markers aligned).
#' @noRd
align_madc_indels <- function(clone_id, allele_id, seqs) {
  seqs <- as.character(seqs)
  lens <- nchar(seqs)
  n_len <- tapply(lens, clone_id, function(x) length(unique(x)))
  to_align <- names(n_len)[n_len > 1]
  if (length(to_align) == 0) return(list(seqs = seqs, n_markers = 0L))

  submat <- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = FALSE)
  for (cl in to_align) {
    rows <- which(clone_id == cl)
    ref_row <- rows[grepl("\\|Ref(_0*1)?$", allele_id[rows])][1]
    if (is.na(ref_row)) ref_row <- rows[1]
    ref <- toupper(seqs[ref_row])
    al <- toupper(seqs[rows])
    aln <- pwalign::pairwiseAlignment(Biostrings::DNAStringSet(al), Biostrings::DNAString(ref),
                                      type = "global", substitutionMatrix = submat,
                                      gapOpening = 5, gapExtension = 2)
    pat <- strsplit(as.character(pwalign::alignedPattern(aln)), "")
    sub <- strsplit(as.character(pwalign::alignedSubject(aln)), "")
    n_ref <- nchar(ref)

    # For each allele: bases inserted before each Ref position (1..n_ref + 1)
    # and the allele base at each Ref position ("-" if deleted)
    parts <- lapply(seq_along(al), function(k) {
      ref_idx <- cumsum(sub[[k]] != "-")   # Ref positions consumed so far
      is_ins <- sub[[k]] == "-"
      ins_at <- ref_idx[is_ins] + 1L       # inserted before this Ref position
      ins <- vapply(seq_len(n_ref + 1L), function(j)
        paste(pat[[k]][is_ins][ins_at == j], collapse = ""), character(1))
      list(ins = ins, base = pat[[k]][!is_ins])
    })
    max_ins <- do.call(pmax, lapply(parts, function(p) nchar(p$ins)))
    # Insertion columns have no Ref base: their gaps are "." (polyRAD does not
    # count "." when it numbers the Ref positions of SNPs)
    pad <- function(x, n) paste0(x, strrep(".", n - nchar(x)))
    seqs[rows] <- vapply(parts, function(p) {
      paste0(paste0(pad(p$ins[seq_len(n_ref)], max_ins[seq_len(n_ref)]), p$base, collapse = ""),
             pad(p$ins[n_ref + 1L], max_ins[n_ref + 1L]))
    }, character(1))
  }
  list(seqs = seqs, n_markers = length(to_align))
}

#' Position of the target SNP in the Ref sequence of each RADdata locus
#'
#' The target SNP is the first column where the Ref and Alt alleles differ, as
#' in `polyRAD::readDArTag()`, but the position counts only bases, not the "-"
#' and "." added by `align_madc_indels()`.
#'
#' @param raddat RADdata object from `polyRAD::readDArTag()`.
#' @return Integer vector with one value per locus (NA if the Ref or Alt allele
#'   is not found).
#' @noRd
target_index_in_ref <- function(raddat) {
  ref <- as.character(raddat$locTable$Ref)
  if (length(ref) == 0) return(rep(NA_integer_, nrow(raddat$locTable)))
  al_names <- colnames(raddat$alleleDepth)
  is_alt <- grepl("\\|Alt(_0*2)?$", al_names)
  alt <- rep(NA_character_, length(ref))
  alt[raddat$alleles2loc[is_alt]] <- as.character(raddat$alleleNucleotides[is_alt])
  mapply(function(r, a) {
    if (is.na(r) || is.na(a) || nchar(r) != nchar(a)) return(NA_integer_)
    r <- strsplit(r, "")[[1]]
    d <- which(r != strsplit(a, "")[[1]])[1]
    if (is.na(d)) return(NA_integer_)
    as.integer(sum(!r[seq_len(d - 1L)] %in% c("-", ".")) + 1L)
  }, ref, alt, USE.NAMES = FALSE)
}

#' Mark the missing ends of shorter alleles as unknown bases
#'
#' Used only for the SNP export (`asSNPs = TRUE`), after the genotypes are
#' called, so it does not change which alleles polyRAD keeps or their
#' genotypes. The bases missing at the start or end of a shorter allele (its
#' leading and trailing "-") are set to "N". polyRAD does not count "N" as a
#' variant, so these ends are not written as deletions and the SNPs inside
#' them get their own records; the alleles with "N" at a SNP are written as an
#' "N" ALT allele of that record. The "." of insertion columns (no Ref base)
#' are kept.
#'
#' @param raddat RADdata object.
#' @return The RADdata object with the terminal "-" of each allele set to "N".
#' @noRd
mask_terminal_gaps <- function(raddat) {
  nuc <- as.character(raddat$alleleNucleotides)
  has_gap <- grepl("^[-.]*-|-[-.]*$", nuc)
  if (!any(has_gap)) return(raddat)
  new_nuc <- raddat$alleleNucleotides
  new_nuc[has_gap] <- vapply(nuc[has_gap], function(x) {
    ch <- strsplit(x, "")[[1]]
    is_gap <- ch %in% c("-", ".")
    bases <- which(!is_gap)
    if (length(bases) == 0) return(x)
    ends <- seq_along(ch) < bases[1] | seq_along(ch) > bases[length(bases)]
    ch[ends & ch == "-"] <- "N"
    paste(ch, collapse = "")
  }, character(1), USE.NAMES = FALSE)
  attributes(new_nuc) <- attributes(raddat$alleleNucleotides)
  raddat$alleleNucleotides <- new_nuc
  raddat
}

#' Add an anchor column to markers whose alleles start with a gap
#'
#' With `asSNPs = TRUE`, polyRAD writes an indel with the base before it
#' (VCF style) and fails ("basic_string") when the indel is in the first
#' column, e.g. a shorter allele of a bottom strand marker, which misses the
#' first bases after the reverse complement. For these markers an "N" column is
#' added before the first column of all alleles and of Ref, and Pos is moved 1
#' bp back, so the indel is written with "N" as its first base. "N" matches any
#' base, so it adds no variant site.
#'
#' @param raddat RADdata object.
#' @return The RADdata object with the anchor column added where needed.
#' @noRd
anchor_leading_gaps <- function(raddat) {
  nuc <- as.character(raddat$alleleNucleotides)
  lead_gap <- substr(nuc, 1, 1) %in% c("-", ".")
  loci <- unique(raddat$alleles2loc[lead_gap])
  if (length(loci) == 0) return(raddat)
  al <- raddat$alleles2loc %in% loci
  new_nuc <- raddat$alleleNucleotides
  new_nuc[al] <- paste0("N", nuc[al])
  attributes(new_nuc) <- attributes(raddat$alleleNucleotides)
  raddat$alleleNucleotides <- new_nuc
  raddat$locTable$Ref[loci] <- paste0("N", raddat$locTable$Ref[loci])
  raddat$locTable$Pos[loci] <- raddat$locTable$Pos[loci] - 1L
  raddat
}

#' Position of the target SNP from the MADC Ref and Alt sequences
#'
#' Used for markers whose Ref and Alt alleles are identical in polyRAD, because
#' polyRAD cut the target SNP away (alleles cut to the shortest allele). The
#' target SNP is the first position where the MADC Ref and Alt differ. polyRAD
#' cuts the end of the MADC sequences and then reverse complements the bottom
#' strand markers, so for these the position is counted from the end of the
#' cut sequence (it can be 0 or negative: the target is before the amplicon
#' start in the VCF). For top strand markers it can be larger than the cut
#' sequence.
#'
#' @param clone_id,allele_id,seqs MADC CloneID, AlleleID and AlleleSequence
#'   columns (without alignment gaps).
#' @param loci Marker IDs.
#' @param bottom Logical, one per marker: is it a bottom strand marker?
#' @param ref_len Length of the Ref sequence of each marker in polyRAD.
#' @return Integer vector with one value per marker (NA if the Ref or Alt
#'   allele is not found or they are identical).
#' @noRd
target_index_from_madc <- function(clone_id, allele_id, seqs, loci, bottom, ref_len) {
  first_row <- function(pattern) vapply(loci, function(L) which(clone_id == L & grepl(pattern, allele_id))[1], integer(1))
  ref <- seqs[first_row("\\|Ref(_0*1)?$")]
  alt <- seqs[first_row("\\|Alt(_0*2)?$")]
  idx <- mapply(function(r, a) {
    if (is.na(r) || is.na(a)) return(NA_integer_)
    n <- min(nchar(r), nchar(a))
    which(strsplit(substr(r, 1, n), "")[[1]] != strsplit(substr(a, 1, n), "")[[1]])[1]
  }, ref, alt, USE.NAMES = FALSE)
  as.integer(ifelse(bottom, ref_len - idx + 1L, idx))
}

#' Run polyRAD::TestOverdispersion() for one set of values
#'
#' polyRAD (<= 2.0.0) errors instead of returning NA when the best value is the
#' minimum of `to_test` ("missing value where TRUE/FALSE needed"), so the call
#' is wrapped in tryCatch and the best value is parsed from its printed
#' "Optimal value is X." line.
#'
#' @param raddat RADdata object.
#' @param to_test Numeric vector of overdispersion values to test.
#' @param verbose Logical. Whether to print the polyRAD output.
#' @return List with `optimal` (polyRAD's optimum, NA when it is at the edge of
#'   `to_test` or polyRAD errored) and `best` (the best tested value, NA if it
#'   could not be parsed).
#' @noRd
run_test_od <- function(raddat, to_test, verbose = TRUE) {
  od <- NULL
  printed <- utils::capture.output(
    od <- tryCatch(polyRAD::TestOverdispersion(raddat, to_test = to_test),
                   error = function(e) e)
  )
  if (verbose) cat(printed, sep = "\n")
  hit <- regmatches(printed, regexpr("(?<=Optimal value is )[0-9.eE+-]+(?=\\.$)", printed, perl = TRUE))
  best <- if (length(hit)) as.numeric(hit[length(hit)]) else NA_real_
  if (inherits(od, "error")) {
    if (is.na(best)) stop(od)
    optimal <- NA_real_
  } else {
    optimal <- od$optimal
  }
  list(optimal = optimal, best = best)
}
