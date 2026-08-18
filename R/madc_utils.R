# Internal shared helpers for madc_plot() and madc_summary().
# Not exported. External packages are fully qualified (ggplot2::, grid::, ...);
# NSE symbols used in ggplot2 aes() are registered in utils::globalVariables()
# in R/utils.R.

#' Read (if needed) and validate a fixed allele ID MADC
#'
#' Accepts a file path or an already-read data.frame, runs check_madc_sanity(),
#' surfaces the checks as messages, and hard-stops unless the file has the
#' required columns and fixed AlleleIDs.
#'
#' @keywords internal
#' @noRd
.read_and_check_madc <- function(madc, verbose = TRUE) {
  if (is.data.frame(madc)) {
    report <- madc
  } else if (is.character(madc) && length(madc) == 1L) {
    if (!file.exists(madc)) stop("MADC file not found: ", madc)
    report <- utils::read.csv(madc, check.names = FALSE)
  } else {
    stop("`madc` must be a single file path or a data.frame.")
  }

  checks <- check_madc_sanity(report)

  msgs <- mapply(function(check, message) if (isTRUE(check)) message[1] else message[2],
                 checks$checks, checks$messages)
  if (verbose) for (i in seq_along(msgs)) message(msgs[i])

  if (!isTRUE(checks$checks[["Columns"]]))
    stop("The MADC file is missing required columns (CloneID, AlleleID, AlleleSequence)")
  if (!isTRUE(checks$checks[["FixAlleleIDs"]]))
    stop("The MADC file does not have fixed AlleleIDs. Please process the MADC file through HapApp before using this function.")

  report
}

#' Compute the shared per-locus / per-sample MADC metric set once
#'
#' Everything downstream (plots and tables) reads from this so visuals and
#' tables never disagree. All matrices are returned on a master CloneID index
#' (file order of unique CloneIDs); markers that get_countsMADC drops (no Ref or
#' no Alt) still appear, with NA in the ratio matrices.
#'
#' @keywords internal
#' @noRd
.madc_metrics <- function(report, min.depth = 10, target.only = FALSE,
                          mhap.min.reads = 1, markers_info = NULL, verbose = FALSE) {

  samp_cols <- names(report)[-(1:3)]
  if (length(samp_cols) < 1L) stop("MADC has no sample columns.")

  X <- data.matrix(report[samp_cols])
  storage.mode(X) <- "double"
  X[is.na(X)] <- 0

  clone  <- as.character(report$CloneID)
  master <- unique(clone)

  is_target <- grepl("\\|(Ref|Alt)_", report$AlleleID)
  is_off    <- grepl("RefMatch_|AltMatch_|[|]Other", report$AlleleID)
  use_rows  <- if (target.only) is_target else rep(TRUE, nrow(report))

  # --- reindex helpers onto the master locus set ---
  reidx <- function(m, fill = NA_real_) {
    out <- matrix(fill, length(master), length(samp_cols),
                  dimnames = list(master, samp_cols))
    if (!is.null(m) && nrow(m) > 0) {
      common <- intersect(rownames(m), master)
      if (length(common)) out[common, colnames(m)] <- m[common, , drop = FALSE]
    }
    out
  }

  # --- Ref/Alt/size/ratio from the existing count engine (paired markers only) ---
  cm <- suppressWarnings(get_countsMADC(madc_object = report,
                                        collapse_matches_counts = !target.only,
                                        verbose = FALSE))
  size_p <- cm$size_matrix
  ref_p  <- cm$ref_matrix
  alt_p  <- size_p - ref_p
  ratio_p <- alt_p / size_p
  ratio_p[!is.finite(ratio_p)] <- NA

  size_matrix <- reidx(size_p)
  ref_matrix  <- reidx(ref_p)
  alt_matrix  <- reidx(alt_p)
  alt_ratio   <- reidx(ratio_p)

  # --- total read depth per (locus x sample), over the used rows ---
  depth_total <- reidx(rowsum(X[use_rows, , drop = FALSE], clone[use_rows]), fill = 0)

  # --- missing mask (TRUE = missing) ---
  missing_mask <- depth_total < min.depth

  # --- number of mHaps present per (locus x sample) ---
  present <- (X[use_rows, , drop = FALSE] >= mhap.min.reads) * 1
  n_mhaps_present <- reidx(rowsum(present, clone[use_rows]), fill = 0)

  # --- panel-level number of mHaps per marker ---
  # n_mhaps_defined  = number of AlleleID rows defined for the marker in the panel
  #                    (independent of read counts).
  # n_mhaps_observed = number of those alleles actually supported (>= mhap.min.reads
  #                    in at least one sample) - empirical run evidence.
  tab <- table(clone[use_rows])
  n_mhaps_marker <- stats::setNames(as.integer(tab[master]), master)
  n_mhaps_marker[is.na(n_mhaps_marker)] <- 0L

  allele_obs <- rowSums(present) >= 1
  obs_tab <- table(clone[use_rows][allele_obs])
  n_mhaps_observed <- stats::setNames(as.integer(obs_tab[master]), master)
  n_mhaps_observed[is.na(n_mhaps_observed)] <- 0L

  # --- off-target reads and fraction ---
  if (!target.only && any(is_off)) {
    offtarget <- reidx(rowsum(X[is_off, , drop = FALSE], clone[is_off]), fill = 0)
  } else {
    offtarget <- matrix(0, length(master), length(samp_cols),
                        dimnames = list(master, samp_cols))
  }
  offtarget_frac <- offtarget / depth_total
  offtarget_frac[!is.finite(offtarget_frac)] <- NA
  # on-target reads = total minus off-target (i.e. the target |Ref/|Alt reads)
  ontarget_depth <- depth_total - offtarget

  # --- chromosome / position ---
  cp <- .madc_chrpos(master, markers_info = markers_info)

  list(ref_matrix = ref_matrix, size_matrix = size_matrix, alt_matrix = alt_matrix,
       alt_ratio = alt_ratio, depth_total = depth_total, ontarget_depth = ontarget_depth,
       missing_mask = missing_mask,
       n_mhaps_present = n_mhaps_present, n_mhaps_marker = n_mhaps_marker,
       n_mhaps_defined = n_mhaps_marker, n_mhaps_observed = n_mhaps_observed,
       offtarget = offtarget, offtarget_frac = offtarget_frac,
       markers = master, chr = cp$chr, pos = cp$pos, samples = samp_cols,
       mhap.min.reads = mhap.min.reads,
       target.only = target.only, min.depth = min.depth)
}

#' Extract Chr/Pos from CloneIDs (furthest-right underscore = Pos) or a lookup
#'
#' @keywords internal
#' @noRd
.madc_chrpos <- function(cloneids, markers_info = NULL) {
  if (!is.null(markers_info)) {
    idcol <- pick_markers_info_id_col(markers_info, cloneids)
    if (!all(c("Chr", "Pos") %in% names(markers_info)))
      stop("`markers_info` must contain 'Chr' and 'Pos' columns.")
    mm  <- match(cloneids, markers_info[[idcol]])
    chr <- as.character(markers_info$Chr[mm])
    pos <- suppressWarnings(as.numeric(markers_info$Pos[mm]))
    return(list(chr = chr, pos = pos))
  }
  parts <- strsplit(cloneids, "_")
  chr <- vapply(parts, function(p)
    if (length(p) >= 2) paste(p[-length(p)], collapse = "_") else NA_character_,
    character(1))
  pos <- suppressWarnings(as.numeric(vapply(parts, function(p)
    if (length(p) >= 2) p[length(p)] else NA_character_, character(1))))
  list(chr = chr, pos = pos)
}

#' Resolve a per-sample grouping vector from a metadata table
#'
#' The sample-id column is auto-detected (sample/Sample/ID/SampleID/...) or the
#' first column; values are matched to the MADC sample column names.
#'
#' @keywords internal
#' @noRd
.madc_group <- function(metadata, group.col, samples) {
  if (is.null(metadata)) return(NULL)
  if (!is.data.frame(metadata)) stop("`metadata` must be a data.frame.")
  if (is.null(group.col) || !group.col %in% names(metadata))
    stop("`group.col` must be a column name in `metadata`.")
  id_candidates <- intersect(c("sample", "Sample", "ID", "SampleID", "sample_id",
                               "Sample_ID"), names(metadata))
  id_col <- if (length(id_candidates)) id_candidates[1] else names(metadata)[1]
  grp <- as.character(metadata[[group.col]][match(samples, metadata[[id_col]])])
  n_match <- sum(!is.na(grp))
  if (all(is.na(grp))) {
    warning("No `metadata` sample IDs matched the MADC sample columns; grouping skipped.")
  } else if (n_match < length(samples)) {
    message(sprintf("Metadata matched %d of %d MADC samples; %d unmatched sample(s) grouped as NA.",
                    n_match, length(samples), length(samples) - n_match))
  }
  grp
}

#' Shared minimal ggplot theme for MADC plots
#'
#' @keywords internal
#' @noRd
.madc_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position  = "right")
}

#' Arrange ggplots/grobs into a layout using base grid (no gridExtra/patchwork)
#'
#' Returns a drawable gTree (does not draw to a device itself).
#'
#' @keywords internal
#' @noRd
.arrange_grobs <- function(grobs, nrow, ncol, heights = NULL, widths = NULL,
                           align_widths = FALSE) {
  hts <- if (is.null(heights)) grid::unit(rep(1, nrow), "null") else heights
  wds <- if (is.null(widths))  grid::unit(rep(1, ncol), "null") else widths

  # Draw on an explicit throwaway device so neither ggplotGrob() text-measurement
  # nor the grid drawing leaves a stray Rplots.pdf. The layout is captured with
  # grid.grab(wrap.grobs = TRUE) - wrapping avoids gtables' shared internal grob
  # names colliding when several panels are combined.
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 7, height = 7, units = "in", res = 100)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)

  gl <- lapply(grobs, function(g) if (grid::is.grob(g)) g else ggplot2::ggplotGrob(g))
  if (align_widths) {                 # line up stacked panels on their x-axes
    has_w <- vapply(gl, function(g) !is.null(g$widths), logical(1))
    if (all(has_w)) {
      maxw <- do.call(grid::unit.pmax, lapply(gl, function(g) g$widths))
      gl <- lapply(gl, function(g) { g$widths <- maxw; g })
    }
  }

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(nrow, ncol, heights = hts, widths = wds)))
  for (i in seq_along(gl)) {
    r  <- ((i - 1) %/% ncol) + 1
    cc <- ((i - 1) %% ncol) + 1
    grid::pushViewport(grid::viewport(layout.pos.row = r, layout.pos.col = cc))
    grid::grid.draw(gl[[i]])
    grid::upViewport()
  }
  grid::upViewport()
  grid::grid.grab(wrap.grobs = TRUE)
}

#' Save a ggplot or an assembled grob to PNG/PDF
#'
#' @keywords internal
#' @noRd
.madc_save <- function(x, file, width = 8, height = 6, dpi = 300) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "pdf") {
    grDevices::pdf(file, width = width, height = height)
  } else {
    grDevices::png(file, width = width, height = height, units = "in", res = dpi)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  if (inherits(x, "ggplot")) print(x) else { grid::grid.newpage(); grid::grid.draw(x) }
  invisible(file)
}
