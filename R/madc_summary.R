#' Summarize a fixed allele ID MADC
#'
#' Computes read-only summary tables for a fixed allele ID MADC file (one
#' processed through HapApp): a per-marker microhaplotype frequency
#' distribution, per-sample and per-marker metric tables (read depth, missing
#' rate, number of mHaps, off-target read fraction), a missingness threshold
#' sweep, an overall summary, and - when a sample metadata table is supplied -
#' the per-sample metrics aggregated by category.
#'
#' @details
#' The input is validated with [check_madc_sanity()]; raw DArT MADC files (or
#' any file without fixed AlleleIDs) are rejected with an error. A locus x sample
#' cell is treated as missing when its total read depth is below `min.depth`.
#' When `target.only = TRUE`, every metric is restricted to the target
#' `|Ref`/`|Alt` alleles (`|RefMatch`/`|AltMatch`/`|Other` are ignored).
#'
#' The per-sample, per-marker, and by-category tables report both `depth` (total
#' reads at the locus - how well the region amplifies) and `target_depth` (reads
#' on the designed `|Ref`/`|Alt` alleles - how well the target markers work),
#' together with `offtarget_fraction`, the share of reads coming from the
#' `|RefMatch`/`|AltMatch`/`|Other` match alleles (i.e. the additional
#' information contributed by non-target haplotypes).
#'
#' @param madc Path to a fixed allele ID MADC file, or an already-read
#'   data.frame.
#' @param target.only Logical; restrict all metrics to the target `|Ref`/`|Alt`
#'   alleles. Default `FALSE`.
#' @param min.depth Minimum total read depth for a locus x sample cell to be
#'   considered present (below this it is missing). Default `10`.
#' @param metadata Optional data.frame mapping samples to a grouping category.
#'   The sample-id column is auto-detected (`sample`/`Sample`/`ID`/`SampleID`) or
#'   taken as the first column; ids are matched to the MADC sample column names.
#' @param group.col Name of the category column in `metadata` (e.g. species or
#'   population). Required to produce the `by_category` table.
#' @param markers_info Optional marker lookup table with an id column
#'   (`CloneID`/`Marker_ID`/`BI_markerID`) plus `Chr` and `Pos`, used when the
#'   CloneIDs are not in `Chr_Pos` format.
#' @param mhap.cap Largest explicit row in the microhaplotype-frequency table;
#'   markers with more mHaps are pooled into a `">mhap.cap"` bucket. Default `10`.
#' @param miss.thresholds Numeric vector of missingness fractions for the
#'   threshold sweep. Default `c(0,5,...,100)/100`.
#' @param mhap.min.reads Minimum reads for a mHap to count as present in a sample
#'   when tallying per-sample mHap counts. Default `1`.
#' @param output.file If not `NULL`, each table is written to
#'   `paste0(output.file, "_<table>.csv")`.
#' @param verbose Logical; print progress/validation messages. Default `TRUE`.
#'
#' @return A named list of data.frames: `mhap_freq`, `per_sample`, `per_marker`,
#'   `missingness`, `overall`, and (when `metadata`/`group.col` are given)
#'   `by_category`.
#'
#' @examples
#' madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
#' tabs <- madc_summary(madc_file, min.depth = 10)
#' tabs$mhap_freq
#' head(tabs$per_sample)
#'
#' @seealso [madc_plot()], [filterMADC()], [check_madc_sanity()]
#' @export
madc_summary <- function(madc,
                         target.only     = FALSE,
                         min.depth       = 10,
                         metadata        = NULL,
                         group.col       = NULL,
                         markers_info    = NULL,
                         mhap.cap        = 10,
                         miss.thresholds = c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100) / 100,
                         mhap.min.reads  = 1,
                         output.file     = NULL,
                         verbose         = TRUE) {

  report <- .read_and_check_madc(madc, verbose = verbose)
  m <- .madc_metrics(report, min.depth = min.depth, target.only = target.only,
                     mhap.min.reads = mhap.min.reads, markers_info = markers_info,
                     verbose = FALSE)

  ## --- per-marker microhaplotype frequency distribution (paralog-flag table) ---
  cnt <- m$n_mhaps_marker[m$n_mhaps_marker >= 2]
  labels <- c(as.character(2:mhap.cap), paste0(">", mhap.cap))
  n_markers <- c(vapply(2:mhap.cap, function(k) sum(cnt == k), integer(1)),
                 sum(cnt > mhap.cap))
  mhap_freq <- data.frame(n_mhaps   = factor(labels, levels = labels),
                          n_markers = n_markers,
                          stringsAsFactors = FALSE)

  ## --- per-sample metrics ---
  # depth = total reads at the locus (how well the region amplifies);
  # target_depth = reads on the designed |Ref/|Alt alleles (how well the target
  # markers work); offtarget_fraction = share of reads from Match/Other alleles.
  per_sample <- data.frame(
    sample             = m$samples,
    depth              = colMeans(m$depth_total),
    target_depth       = colMeans(m$ontarget_depth),
    missing_rate       = colMeans(m$missing_mask),
    n_mhaps            = colMeans(m$n_mhaps_present),
    offtarget_fraction = colMeans(m$offtarget_frac, na.rm = TRUE),
    stringsAsFactors   = FALSE)
  rownames(per_sample) <- NULL

  ## --- per-marker metrics ---
  per_marker <- data.frame(
    CloneID            = m$markers,
    Chr                = m$chr,
    Pos                = m$pos,
    depth              = rowMeans(m$depth_total),
    target_depth       = rowMeans(m$ontarget_depth),
    missing_rate       = rowMeans(m$missing_mask),
    n_mhaps            = as.integer(m$n_mhaps_marker),
    offtarget_fraction = rowMeans(m$offtarget_frac, na.rm = TRUE),
    stringsAsFactors   = FALSE)
  rownames(per_marker) <- NULL
  if (all(is.na(per_marker$Chr)) && all(is.na(per_marker$Pos))) {
    per_marker$Chr <- NULL
    per_marker$Pos <- NULL
  }

  ## --- missingness threshold sweep ---
  locus_miss <- rowMeans(m$missing_mask)
  samp_miss  <- colMeans(m$missing_mask)
  missingness <- data.frame(
    threshold_pct    = miss.thresholds * 100,
    loci_removed     = vapply(miss.thresholds, function(t) sum(locus_miss >  t), integer(1)),
    loci_retained    = vapply(miss.thresholds, function(t) sum(locus_miss <= t), integer(1)),
    samples_removed  = vapply(miss.thresholds, function(t) sum(samp_miss  >  t), integer(1)),
    samples_retained = vapply(miss.thresholds, function(t) sum(samp_miss  <= t), integer(1)),
    stringsAsFactors = FALSE)

  ## --- overall one-row summary ---
  overall <- data.frame(
    n_samples               = length(m$samples),
    n_markers               = length(m$markers),
    n_mhaps_total           = sum(m$n_mhaps_marker),
    mean_depth              = mean(m$depth_total),
    overall_missing_rate    = mean(m$missing_mask),
    mean_offtarget_fraction = mean(m$offtarget_frac, na.rm = TRUE),
    median_mhaps_per_marker = stats::median(m$n_mhaps_marker),
    stringsAsFactors        = FALSE)

  out <- list(mhap_freq = mhap_freq, per_sample = per_sample,
              per_marker = per_marker, missingness = missingness, overall = overall)

  ## --- optional by-category aggregation ---
  grp <- .madc_group(metadata, group.col, m$samples)
  if (!is.null(grp) && !all(is.na(grp))) {
    keep <- !is.na(grp)
    g    <- grp[keep]
    cats <- sort(unique(g))
    agg  <- function(v, f) vapply(cats, function(cc) f(v[keep][g == cc]), numeric(1))
    out$by_category <- data.frame(
      group                   = cats,
      n_samples               = as.integer(table(g)[cats]),
      depth_mean              = agg(per_sample$depth, mean),
      target_depth_mean       = agg(per_sample$target_depth, mean),
      missing_rate_mean       = agg(per_sample$missing_rate, mean),
      n_mhaps_mean            = agg(per_sample$n_mhaps, mean),
      offtarget_fraction_mean = agg(per_sample$offtarget_fraction,
                                    function(x) mean(x, na.rm = TRUE)),
      stringsAsFactors        = FALSE)
  }

  ## --- optional CSV output ---
  if (!is.null(output.file)) {
    for (nm in names(out))
      utils::write.csv(out[[nm]], paste0(output.file, "_", nm, ".csv"), row.names = FALSE)
    if (verbose) message("Wrote ", length(out), " summary tables to ",
                         output.file, "_*.csv")
  }

  out
}
