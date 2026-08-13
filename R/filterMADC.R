#' Filter MADC Files
#'
#' Filter and process MADC files to remove low quality microhaplotypes
#'
#' @details
#' This function filters a fixed allele ID MADC file (one processed through HapApp so that
#' AlleleIDs carry unique \code{|Ref_0001}/\code{|Alt_0002}-style suffixes). The input is
#' validated with \code{\link{check_madc_sanity}}; raw DArT MADC files (or any file without
#' fixed AlleleIDs) are rejected with an error.
#'
#' Filtering is designed so that a locus never loses its target \code{|Ref}/\code{|Alt} pair:
#' the per-locus read-depth window (\code{min.locus.depth}/\code{max.locus.depth}) removes
#' whole loci (Ref and Alt together), and the other filters
#' (\code{max.mhaps.per.loci}, \code{min.ind.with.reads}, \code{target.only}) only prune
#' non-target \code{|RefMatch}/\code{|AltMatch}/\code{|Other} mhaps. The filtered data can
#' optionally be saved to a file.
#'
#'@import dplyr
#'@importFrom utils read.csv write.csv
#'
#'@param madc_file Path to the fixed allele ID MADC file to be filtered
#'@param min.locus.depth Minimum mean read depth per sample at a locus. Depth is the sum of all mhap reads at the locus (\code{|Ref}, \code{|Alt}, \code{|RefMatch}, \code{|AltMatch}, and \code{|Other}) divided by the number of samples. Loci below this threshold are removed entirely (Ref and Alt together). If NULL, no minimum is applied.
#'@param max.locus.depth Maximum mean read depth per sample at a locus (computed as for \code{min.locus.depth}). Loci above this threshold are removed entirely, e.g. to exclude paralogous over-amplification. If NULL, no maximum is applied.
#'@param max.mhaps.per.loci Maximum number of mhaps per target loci. At loci whose total number of mhaps exceeds the \code{max.mhaps.per.loci} threshold, only the target \code{|Ref} and \code{|Alt} alleles are retained (the \code{|RefMatch}, \code{|AltMatch}, and \code{|Other} alleles are removed).
#'@param min.reads.per.site Minimum number of reads for a site (mhap in a sample) to count toward \code{min.ind.with.reads}.
#'@param min.ind.with.reads Minimum number of individuals with \code{min.reads.per.site} reads required to retain a mhap. Target \code{|Ref}/\code{|Alt} alleles are always retained; only non-target (\code{|RefMatch}/\code{|AltMatch}/\code{|Other}) mhaps are removed by this filter.
#'@param target.only Logical indicating whether to retain only the target \code{|Ref} and \code{|Alt} alleles (dropping \code{|RefMatch}, \code{|AltMatch}, and \code{|Other})
#@param plot.summary Logical indicating whether to plot summary statistics
#'@param output.file Path to save the filtered data (if NULL, data will not be saved)
#'
#'@return data.frame or saved csv file
#'
#'@examples
#' #Example
#'
#' #Example MADC
#' madc_file <- system.file("example_MADC_FixedAlleleID.csv", package="BIGr")
#'
#' #Remove mhaps exceeding 3 per target region including the ref and alt target mhaps
#' filtered_df <- filterMADC(madc_file,
#'                          min.mean.reads = NULL,
#'                          max.mean.reads = NULL,
#'                          max.mhaps.per.loci = 3,
#'                          min.reads.per.site = 1,
#'                          min.ind.with.reads = NULL,
#'                          target.only = FALSE,
#'                          output.file = NULL)
#'
#'
#'
#'@export
filterMADC <- function(madc_file,
                       min.locus.depth = NULL,
                       max.locus.depth = NULL,
                       max.mhaps.per.loci = NULL,
                       min.reads.per.site = 1,
                       min.ind.with.reads = NULL,
                       target.only = FALSE,
                       #plot.summary = FALSE,
                       output.file = NULL) {


  #Read and sanity-check (fixed allele ID MADC only)
  report <- read.csv(madc_file, check.names = FALSE)
  checks <- check_madc_sanity(report)

  #Surface all sanity checks as informational messages (non-blocking)
  msgs <- mapply(function(check, message) if (isTRUE(check)) message[1] else message[2],
                 checks$checks, checks$messages)
  for (i in seq_along(msgs)) message(msgs[i])

  #Hard stops: only accept fixed allele ID MADC files with the required columns
  if (!isTRUE(checks$checks[["Columns"]]))
    stop("The MADC file is missing required columns (CloneID, AlleleID, AlleleSequence)")

  if (!isTRUE(checks$checks[["FixAlleleIDs"]]))
    stop("The MADC file does not have fixed AlleleIDs. Please process the MADC file through HapApp before using this function.")

  #Fixed allele ID MADC has a standard layout: columns 1:3 are the ID columns
  #(AlleleID, CloneID, AlleleSequence) and columns 4:ncol are numeric samples.
  filtered_df <- report

  ## Filtering

  #Per-locus depth window (runs first, on the full mhap set so paralogous
  #over-amplification is measured before any trimming). Depth is the mean reads
  #per sample at a locus = (sum of all mhap reads at the locus) / number of samples.
  #Loci outside [min.locus.depth, max.locus.depth] are removed whole (Ref and Alt
  #together), so retained loci always keep their Ref/Alt pair intact.
  if (!is.null(min.locus.depth) || !is.null(max.locus.depth)) {
    n_samples <- ncol(filtered_df) - 3
    row_totals <- rowSums(filtered_df[, -c(1:3), drop = FALSE], na.rm = TRUE)
    locus_depth <- tapply(row_totals, filtered_df$CloneID, sum) / n_samples

    keep_mask <- rep(TRUE, length(locus_depth))
    if (!is.null(min.locus.depth)) {
      message("Filtering for minimum mean read depth per sample per locus")
      keep_mask <- keep_mask & (locus_depth >= min.locus.depth)
    }
    if (!is.null(max.locus.depth)) {
      message("Filtering for maximum mean read depth per sample per locus")
      keep_mask <- keep_mask & (locus_depth <= max.locus.depth)
    }
    keep_loci <- names(locus_depth)[keep_mask]
    filtered_df <- filtered_df[filtered_df$CloneID %in% keep_loci, ]
  }

  #Remove refmatch, altmatch, and other if wanted (retain only the target Ref and Alt haplotypes)
  if (target.only) {
    message("Retaining target markers only")
    filtered_df <- filtered_df[grepl("\\|(Ref|Alt)_", filtered_df$AlleleID), ]
  }

  #Max mhaps per loci
  if (!is.null(max.mhaps.per.loci)) {
    message("Filtering for maximum number of matching mhaps per target loci")
    #Count ALL mhaps per locus; at loci exceeding the max, retain only the target Ref/Alt alleles
    clone_ids_to_target <- filtered_df %>%
      group_by(CloneID) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      filter(Count > max.mhaps.per.loci) %>%
      pull(CloneID)

    filtered_df <- filtered_df %>%
      filter(
        !( # "keep rows that DO NOT match both conditions"
          CloneID %in% clone_ids_to_target &  # Condition 1: locus is over the threshold
            !grepl("\\|(Ref|Alt)_", AlleleID) # Condition 2: AlleleID is NOT a target Ref/Alt
        )
      )
  }

  #Min individuals with reads
  if (!is.null(min.ind.with.reads)) {
    message("Filtering for minimum number of individuals with reads per site")
    message(paste0("Minimum number of individuals with reads per site: ", min.ind.with.reads))
    message(paste0("Minimum number of reads per site: ", min.reads.per.site))

    #Getting colnames
    cols_to_check <- colnames(filtered_df)[-(1:3)]

    filtered_df <- filtered_df %>%
      rowwise() %>%  # Process data row by row
      mutate(
        # For each row, count how many of the 'cols_to_check' meet the criterion
        qualifying_sites_count = sum(
          c_across(all_of(cols_to_check)) >= min.reads.per.site,
          na.rm = TRUE # Treats NAs in data as not meeting the criterion
        )
      ) %>%
      ungroup() %>% # Always ungroup after rowwise operations
      # Always retain target Ref/Alt; only prune non-target (Match/Other) mhaps that
      # fail the 'min.ind.with.reads' threshold, so a locus never loses its Ref/Alt pair
      filter(grepl("\\|(Ref|Alt)_", AlleleID) | qualifying_sites_count >= min.ind.with.reads) %>%
      # Optionally, remove the temporary count column if it's no longer needed
      select(-qualifying_sites_count)
  }

  #Plots
  #if (plot.summary) {
  #  message("Plotting summary statistics")
  #  #Plot mean read depth
  #  mean_reads <- rowMeans(filtered_df[, -c(1:3)], na.rm = TRUE)
  #  hist(mean_reads, main = "Mean Read Depth", xlab = "Mean Reads", ylab = "Frequency")

  #  #Plot number of Altmatch and Refmatch mhaps per target loci
  #  altmatch_counts <- filtered_df %>%
  #    filter(grepl("\\|AltMatch", AlleleID)) %>%
  #    group_by(CloneID) %>%
  #    summarise(Count = n(), .groups = 'drop')

  #  refmatch_counts <- filtered_df %>%
  #    filter(grepl("\\|RefMatch", AlleleID)) %>%
  #    group_by(CloneID) %>%
  #    summarise(Count = n(), .groups = 'drop')

  #  barplot(cbind(altmatch_counts$Count, refmatch_counts$Count), beside = TRUE,
  #          names.arg = altmatch_counts$CloneID, main = "Number of AltMatch and RefMatch Mhaps",
  #          xlab = "Clone ID", ylab = "Count")

    #Plot density of number of CloneID per site on a marker distribution plot

  #}

  #Save the output to disk if file name provided
  if (!is.null(output.file)) {
    message("Saving filtered data to file")
    write.csv(filtered_df, paste0(output.file,".csv"), row.names = FALSE)
  } else {
    message("No output file provided. Returning filtered data.")
    return(filtered_df)
  }

}
