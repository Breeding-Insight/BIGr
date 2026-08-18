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
#' @param depth.thresholds Numeric vector of read-depth cutoffs for the
#'   `depth_sweep` table. Default `c(1, 5, 10, 20, 50, 100)`.
#' @param mhap.min.reads Minimum reads for a mHap to count as present in a sample
#'   when tallying per-sample mHap counts. Default `1`.
#' @param min.locus.depth,max.locus.depth Marker-QC read-depth window (mean reads
#'   per sample at the locus). A marker below `min.locus.depth` is `FAIL`ed (low
#'   depth); one above `max.locus.depth` is `FLAG`ged (over-amplified). Names/units
#'   match [filterMADC()], so a `FAIL` on `min.locus.depth` is the marker
#'   `filterMADC` would drop at the same value. `min.locus.depth` defaults to
#'   `min.depth`; `max.locus.depth` is off (`NULL`) by default. Set either to
#'   `NULL` to disable that criterion.
#' @param max.missing Marker-QC maximum per-marker missing rate; markers above it
#'   are `FAIL`ed (catches the high-dropout case a mean-depth cut alone misses).
#'   `NULL` (default) disables it.
#' @param max.offtarget Marker-QC maximum read-weighted off-target fraction;
#'   markers above it are `FLAG`ged. `NULL` (default) disables it.
#' @param max.mhaps.per.loci Marker-QC maximum number of defined mHaps; markers
#'   above it are `FLAG`ged (paralog-suspect). `NULL` (default) disables it.
#' @param output.file If not `NULL`, the tables are written here: as
#'   `paste0(output.file, "_<table>.csv")` when `output.format = "csv"`, or as a
#'   single `.xlsx` workbook (one sheet per table) when `output.format = "xlsx"`.
#' @param output.format One of `"csv"` (default; one CSV per table) or `"xlsx"`
#'   (a single multi-sheet workbook, requires the `writexl` package).
#' @param verbose Logical; print progress/validation messages. Default `TRUE`.
#'
#' @return A named list of data.frames: `mhap_freq`, `per_sample`, `per_marker`,
#'   `missingness`, `depth_sweep`, `marker_qc`, `overall`, and (when
#'   `metadata`/`group.col` are given) `by_category`. Depth is summarized by both
#'   the mean and robust statistics (`median_depth`, `depth_mad`, `depth_q10`,
#'   `depth_q90`); off-target is reported both as the mean per-locus ratio
#'   (`offtarget_fraction`) and read-weighted (`offtarget_fraction_weighted`); and
#'   mHap counts are split into `n_mhaps_defined` (allele rows in the panel) and
#'   `n_mhaps_observed` (alleles actually supported by reads).
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
                         depth.thresholds = c(1, 5, 10, 20, 50, 100),
                         mhap.min.reads  = 1,
                         min.locus.depth    = min.depth,
                         max.locus.depth    = NULL,
                         max.missing        = NULL,
                         max.offtarget      = NULL,
                         max.mhaps.per.loci = NULL,
                         output.file     = NULL,
                         output.format   = c("csv", "xlsx"),
                         verbose         = TRUE) {

  output.format <- match.arg(output.format)
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

  ## --- shared helpers -------------------------------------------------------
  # depth       = total reads at the locus (how well the region amplifies);
  # target_depth= reads on the designed |Ref/|Alt alleles (how well the targets
  #               work); offtarget_fraction is the mean of the per-locus off-target
  #               ratios, offtarget_fraction_weighted is the read-weighted share
  #               sum(offtarget)/sum(total) (the "% of reads off target" number).
  q10 <- function(v) stats::quantile(v, 0.10, names = FALSE, na.rm = TRUE)
  q90 <- function(v) stats::quantile(v, 0.90, names = FALSE, na.rm = TRUE)
  safe_ratio <- function(a, b) { r <- a / b; r[!is.finite(r)] <- NA; r }

  ## --- per-sample metrics ---
  per_sample <- data.frame(
    sample                      = m$samples,
    depth                       = colMeans(m$depth_total),
    median_depth                = apply(m$depth_total, 2, stats::median),
    depth_mad                   = apply(m$depth_total, 2, stats::mad),
    depth_q10                   = apply(m$depth_total, 2, q10),
    depth_q90                   = apply(m$depth_total, 2, q90),
    target_depth                = colMeans(m$ontarget_depth),
    total_reads                 = colSums(m$depth_total),
    target_reads                = colSums(m$ontarget_depth),
    missing_rate                = colMeans(m$missing_mask),
    call_rate                   = 1 - colMeans(m$missing_mask),
    n_mhaps_observed            = colMeans(m$n_mhaps_present),
    offtarget_fraction          = colMeans(m$offtarget_frac, na.rm = TRUE),
    offtarget_fraction_weighted = safe_ratio(colSums(m$offtarget), colSums(m$depth_total)),
    stringsAsFactors            = FALSE)
  rownames(per_sample) <- NULL

  ## --- per-marker metrics ---
  per_marker <- data.frame(
    CloneID                     = m$markers,
    Chr                         = m$chr,
    Pos                         = m$pos,
    depth                       = rowMeans(m$depth_total),
    median_depth                = apply(m$depth_total, 1, stats::median),
    depth_mad                   = apply(m$depth_total, 1, stats::mad),
    depth_q10                   = apply(m$depth_total, 1, q10),
    depth_q90                   = apply(m$depth_total, 1, q90),
    target_depth                = rowMeans(m$ontarget_depth),
    missing_rate                = rowMeans(m$missing_mask),
    call_rate                   = 1 - rowMeans(m$missing_mask),
    n_mhaps_defined             = as.integer(m$n_mhaps_defined),
    n_mhaps_observed            = as.integer(m$n_mhaps_observed),
    offtarget_fraction          = rowMeans(m$offtarget_frac, na.rm = TRUE),
    offtarget_fraction_weighted = safe_ratio(rowSums(m$offtarget), rowSums(m$depth_total)),
    stringsAsFactors            = FALSE)
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

  ## --- depth threshold sweep (mirrors the missingness sweep) ---
  cell_depth <- m$depth_total
  depth_sweep <- data.frame(
    depth_threshold            = depth.thresholds,
    cells_passing              = vapply(depth.thresholds, function(t) mean(cell_depth >= t), numeric(1)),
    loci_passing_90pct_samples = vapply(depth.thresholds, function(t)
                                          sum(rowMeans(cell_depth >= t) >= 0.90), integer(1)),
    loci_passing_95pct_samples = vapply(depth.thresholds, function(t)
                                          sum(rowMeans(cell_depth >= t) >= 0.95), integer(1)),
    samples_passing_90pct_loci = vapply(depth.thresholds, function(t)
                                          sum(colMeans(cell_depth >= t) >= 0.90), integer(1)),
    stringsAsFactors           = FALSE)

  ## --- per-marker QC status (advisory; does NOT filter) ---
  # Thresholds reuse filterMADC's vocabulary so a marker FAILing `min.locus.depth`
  # here is the same marker filterMADC drops at that threshold (both key off the
  # mean per-sample locus depth). Each threshold set to NULL disables its criterion.
  n_mk      <- length(m$markers)
  qc_reason <- character(n_mk)
  fail <- logical(n_mk); flag <- logical(n_mk)
  add_qc <- function(cond, label, is_fail) {
    cond[is.na(cond)] <- FALSE
    if (is_fail) fail <<- fail | cond else flag <<- flag | cond
    qc_reason[cond] <<- ifelse(nzchar(qc_reason[cond]),
                               paste0(qc_reason[cond], "; ", label), label)
  }
  if (!is.null(min.locus.depth))
    add_qc(per_marker$depth < min.locus.depth,       sprintf("low depth (<%g)", min.locus.depth), TRUE)
  if (!is.null(max.missing))
    add_qc(per_marker$missing_rate > max.missing,    sprintf("high missing (>%g)", max.missing), TRUE)
  if (!is.null(max.locus.depth))
    add_qc(per_marker$depth > max.locus.depth,       sprintf("over-amplified (>%g)", max.locus.depth), FALSE)
  if (!is.null(max.offtarget))
    add_qc(per_marker$offtarget_fraction_weighted > max.offtarget,
                                                     sprintf("high off-target (>%g)", max.offtarget), FALSE)
  if (!is.null(max.mhaps.per.loci))
    add_qc(per_marker$n_mhaps_defined > max.mhaps.per.loci,
                                                     sprintf("high mHaps (>%g, paralog-suspect)", max.mhaps.per.loci), FALSE)
  marker_qc <- data.frame(
    CloneID                     = m$markers,
    Chr                         = m$chr,
    Pos                         = m$pos,
    mean_depth                  = per_marker$depth,
    call_rate                   = 1 - rowMeans(m$missing_mask),
    offtarget_fraction_weighted = per_marker$offtarget_fraction_weighted,
    n_mhaps_defined             = as.integer(m$n_mhaps_defined),
    n_mhaps_observed            = as.integer(m$n_mhaps_observed),
    qc_status                   = ifelse(fail, "FAIL", ifelse(flag, "FLAG", "PASS")),
    qc_reason                   = qc_reason,
    stringsAsFactors            = FALSE)
  rownames(marker_qc) <- NULL
  if (all(is.na(marker_qc$Chr)) && all(is.na(marker_qc$Pos))) {
    marker_qc$Chr <- NULL; marker_qc$Pos <- NULL
  }

  ## --- overall one-row summary ---
  overall <- data.frame(
    n_samples                        = length(m$samples),
    n_markers                        = length(m$markers),
    n_mhaps_total                    = sum(m$n_mhaps_marker),
    mean_depth                       = mean(m$depth_total),
    median_depth                     = stats::median(m$depth_total),
    overall_missing_rate             = mean(m$missing_mask),
    overall_call_rate                = 1 - mean(m$missing_mask),
    mean_offtarget_fraction          = mean(m$offtarget_frac, na.rm = TRUE),
    read_weighted_offtarget_fraction = safe_ratio(sum(m$offtarget), sum(m$depth_total)),
    median_mhaps_per_marker          = stats::median(m$n_mhaps_marker),
    stringsAsFactors                 = FALSE)

  out <- list(mhap_freq = mhap_freq, per_sample = per_sample,
              per_marker = per_marker, missingness = missingness,
              depth_sweep = depth_sweep, marker_qc = marker_qc, overall = overall)

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

  ## --- optional output: per-table CSVs or a single multi-sheet workbook ---
  if (!is.null(output.file)) {
    if (output.format == "xlsx") {
      if (!requireNamespace("writexl", quietly = TRUE))
        stop("output.format = 'xlsx' requires the 'writexl' package. Install it with install.packages('writexl').")
      path <- if (grepl("\\.xlsx$", output.file, ignore.case = TRUE)) output.file else paste0(output.file, ".xlsx")
      writexl::write_xlsx(out, path)
      if (verbose) message("Wrote ", length(out), " summary tables as sheets to ", path)
    } else {
      for (tab in names(out))
        utils::write.csv(out[[tab]], paste0(output.file, "_", tab, ".csv"), row.names = FALSE)
      if (verbose) message("Wrote ", length(out), " summary tables to ",
                           output.file, "_*.csv")
    }
  }

  out
}
