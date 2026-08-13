#' Plot a fixed allele ID MADC
#'
#' Visualize a fixed allele ID MADC file (one processed through HapApp) with a
#' choice of plots selected via `plot.type`: a PCA of alternative read ratios, a
#' linear marker-distribution genome plot, a markers x samples read-depth/#mHaps
#' heatmap, a per-sample missing-data boxplot, and a `circos` circular plot with
#' concentric marker-position, #mHaps, and read-depth tracks. Requesting one type
#' returns that plot; requesting several (non-circos) returns a multi-panel figure.
#'
#' @details
#' The input is validated with [check_madc_sanity()]; raw DArT MADC files are
#' rejected with an error. All plots honor `target.only` (restrict to `|Ref`/
#' `|Alt`) and `min.depth` (locus x sample cells below this depth are treated as
#' missing). Plotting uses only `ggplot2` and base `grid`; the `circos` plot
#' additionally requires the suggested `circlize` package and, being a base
#' graphics figure, must be requested on its own (it cannot be embedded in a
#' multi-panel).
#'
#' @param madc Path to a fixed allele ID MADC file, or an already-read data.frame.
#' @param plot.type One or more of `"pca"`, `"marker"`, `"heatmap"`, `"missing"`,
#'   `"circos"`. A single value returns that plot; several (excluding `"circos"`)
#'   return a multi-panel figure. `"circos"` must be requested on its own.
#' @param target.only Logical; restrict to target `|Ref`/`|Alt` alleles. Default `FALSE`.
#' @param min.depth Minimum locus x sample read depth to be considered present. Default `10`.
#' @param metadata Optional sample-metadata data.frame (see [madc_summary()]) used
#'   to color the PCA / group the missing-data boxplot.
#' @param group.col Category column in `metadata` (e.g. species/population).
#' @param markers_info Optional marker lookup (id + `Chr` + `Pos`) for the marker/circos plot.
#' @param pc.x,pc.y Principal components to plot (PCA). Defaults `1`, `2`.
#' @param loci.miss.max,sample.miss.max Drop loci/samples with missingness above
#'   these fractions before PCA. Defaults `0.5`.
#' @param fill Heatmap fill metric: `"depth"` (default) or `"mhaps"`.
#' @param facet.chrom Logical; facet the heatmap by chromosome. Default `FALSE`.
#' @param max.loci,max.samples Optional caps; larger heatmaps are subsampled (with a note).
#' @param mhap.min.reads Minimum reads for a mHap to count as present. Default `1`.
#' @param output.file If not `NULL`, the figure is also saved here (`.png`/`.pdf`).
#' @param width,height,dpi Saved-figure dimensions (inches) and resolution.
#' @param verbose Logical; print progress/validation messages. Default `TRUE`.
#'
#' @return For a single `plot.type`, a `ggplot`; for several, an (invisible) list
#'   with `plots` (named list) and `panel` (the assembled multi-panel grob).
#'   `"circos"` draws to the active device and returns `NULL` invisibly.
#'
#' @examples
#' madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
#' p <- madc_plot(madc_file, plot.type = "missing")
#' \donttest{
#' madc_plot(madc_file, plot.type = c("pca", "heatmap"))
#' }
#'
#' @seealso [madc_summary()], [filterMADC()]
#' @export
madc_plot <- function(madc,
                      plot.type      = c("pca", "marker", "heatmap", "missing"),
                      target.only    = FALSE,
                      min.depth      = 10,
                      metadata       = NULL,
                      group.col      = NULL,
                      markers_info   = NULL,
                      pc.x           = 1,
                      pc.y           = 2,
                      loci.miss.max  = 0.5,
                      sample.miss.max = 0.5,
                      fill           = c("depth", "mhaps"),
                      facet.chrom    = FALSE,
                      max.loci       = NULL,
                      max.samples    = NULL,
                      mhap.min.reads = 1,
                      output.file    = NULL,
                      width = 8, height = 6, dpi = 300,
                      verbose        = TRUE) {

  plot.type <- match.arg(plot.type, c("pca", "marker", "heatmap", "missing", "circos"),
                         several.ok = TRUE)
  fill   <- match.arg(fill)

  report <- .read_and_check_madc(madc, verbose = verbose)
  m <- .madc_metrics(report, min.depth = min.depth, target.only = target.only,
                     mhap.min.reads = mhap.min.reads, markers_info = markers_info,
                     verbose = FALSE)
  grp <- .madc_group(metadata, group.col, m$samples)

  # circos is a standalone base-graphics (circlize) figure - it cannot be embedded
  # in a grid multi-panel, so it must be requested on its own.
  if ("circos" %in% plot.type) {
    if (length(plot.type) > 1)
      stop("plot.type = 'circos' must be requested on its own (it renders as a standalone circular figure).")
    save_to <- output.file
    # non-interactive with no output.file: capture to a temp file rather than
    # letting base graphics leave a stray Rplots.pdf.
    if (is.null(save_to) && !interactive()) {
      save_to <- tempfile(fileext = ".png")
      if (verbose) message("No graphics device / output.file; circos written to ", save_to)
    }
    if (!is.null(save_to)) {
      ext <- tolower(tools::file_ext(save_to))
      if (ext == "pdf") grDevices::pdf(save_to, width = width, height = height)
      else grDevices::png(save_to, width = width, height = height, units = "in", res = dpi)
      on.exit(grDevices::dev.off(), add = TRUE)
    }
    .madc_plot_circos(m)
    return(invisible(save_to))
  }

  build <- list(
    pca     = function() .madc_plot_pca(m, pc.x, pc.y, loci.miss.max, sample.miss.max, grp),
    marker  = function() .madc_plot_marker(m, facet.chrom = facet.chrom),
    heatmap = function() .madc_plot_heatmap(m, fill, facet.chrom, max.loci, max.samples, verbose),
    missing = function() .madc_plot_missing(m, grp)
  )

  plots <- lapply(plot.type, function(pt) build[[pt]]())
  names(plots) <- plot.type

  if (length(plots) == 1L) {
    result <- plots[[1]]
    if (!is.null(output.file)) {
      .madc_save(result, output.file, width, height, dpi)
      return(invisible(result))   # saved to file: don't auto-print the object
    }
    return(result)
  }

  n <- length(plots)
  ncol <- 2L
  nrow <- ceiling(n / ncol)
  panel <- .arrange_grobs(plots, nrow = nrow, ncol = ncol)
  if (!is.null(output.file)) .madc_save(panel, output.file, width, height, dpi)
  invisible(list(plots = plots, panel = panel))
}

# ---- PCA -------------------------------------------------------------------
#' @keywords internal
#' @noRd
.madc_plot_pca <- function(m, pc.x = 1, pc.y = 2,
                           loci.miss.max = 0.5, sample.miss.max = 0.5, grp = NULL) {
  mat <- t(m$alt_ratio)                                   # samples x loci
  mat <- mat[, colMeans(is.na(mat)) <= loci.miss.max, drop = FALSE]
  mat <- mat[rowMeans(is.na(mat)) <= sample.miss.max, , drop = FALSE]
  if (nrow(mat) < 2 || ncol(mat) < 2)
    stop("Too few samples/loci pass the missingness filters for PCA.")

  # per-locus mean imputation
  col_means <- colMeans(mat, na.rm = TRUE)
  na_idx <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(na_idx)) mat[na_idx] <- col_means[na_idx[, "col"]]
  # drop zero-variance loci
  keep <- apply(mat, 2, function(x) stats::sd(x) > 0)
  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) < 2) stop("No informative (non-constant) loci remain for PCA.")

  pc <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  df <- data.frame(sample = rownames(mat),
                   PCx = pc$x[, pc.x], PCy = pc$x[, pc.y],
                   stringsAsFactors = FALSE)
  if (!is.null(grp)) df$group <- grp[match(df$sample, m$samples)]

  aes_pca <- if (!is.null(grp))
    ggplot2::aes(x = PCx, y = PCy, color = group) else
    ggplot2::aes(x = PCx, y = PCy)

  ggplot2::ggplot(df, aes_pca) +
    ggplot2::geom_point(size = 2.5, alpha = 0.9) +
    ggplot2::labs(title = "PCA of alternative read ratios",
                  x = sprintf("PC%d (%.1f%%)", pc.x, 100 * ve[pc.x]),
                  y = sprintf("PC%d (%.1f%%)", pc.y, 100 * ve[pc.y]),
                  color = NULL) +
    .madc_theme()
}

# ---- Marker distribution (linear, positions only) --------------------------
#' @keywords internal
#' @noRd
.madc_plot_marker <- function(m, facet.chrom = FALSE) {
  df <- .madc_marker_df(m)
  chr_levels <- rev(sort(unique(df$chr)))
  df$chr <- factor(df$chr, levels = chr_levels)
  backbone <- data.frame(chr = factor(chr_levels, levels = chr_levels),
                         len = vapply(chr_levels, function(c) max(df$pos[df$chr == c]),
                                      numeric(1)), stringsAsFactors = FALSE)
  mb <- function(b) format(round(b / 1e6, 2), big.mark = ",", trim = TRUE)

  ggplot2::ggplot(df) +
    ggplot2::geom_segment(data = backbone,
      ggplot2::aes(x = 0, xend = len, y = chr, yend = chr),
      color = "grey75", linewidth = 3) +
    ggplot2::geom_point(ggplot2::aes(x = pos, y = chr), shape = 124, size = 4,
                        color = "steelblue") +
    ggplot2::scale_x_continuous("Position (Mb)", labels = mb,
                                expand = ggplot2::expansion(mult = c(0.01, 0.03))) +
    ggplot2::labs(title = "Marker distribution", y = NULL) + .madc_theme()
}

# ---- Circos (circlize, three concentric tracks) ----------------------------
#' @keywords internal
#' @noRd
.madc_plot_circos <- function(m, paralog.flag = 10) {
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("plot.type = 'circos' requires the 'circlize' package. Install it with install.packages('circlize').")

  df <- .madc_marker_df(m)
  df$chr <- factor(df$chr, levels = sort(unique(df$chr)))
  df <- df[order(df$chr, df$pos), , drop = FALSE]

  # per-sector bar half-width scaled to marker spacing (so bars are visible but not overlapping)
  bar_hw <- max(vapply(split(df$pos, df$chr), function(p)
    if (length(p) > 1) stats::median(diff(sort(p))) * 0.35 else diff(range(df$pos)) * 0.01,
    numeric(1)))

  dq <- stats::quantile(df$depth, c(0, 0.5, 1), na.rm = TRUE)
  if (dq[2] <= dq[1]) dq[2] <- dq[1] + 1e-6
  if (dq[3] <= dq[2]) dq[3] <- dq[2] + 1e-6
  col_depth <- circlize::colorRamp2(dq, c("#2166ac", "#f7f7f7", "#b2182b"))

  # bar drawer: a rectangle from 0 to value at each marker position
  bars <- function(x, y, cols) circlize::circos.rect(
    x - bar_hw, 0, x + bar_hw, y, col = cols, border = NA)

  op <- graphics::par(mar = c(1, 1, 2, 1)); on.exit(graphics::par(op), add = TRUE)
  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  circlize::circos.par(gap.degree = if (nlevels(df$chr) == 1) 14 else 4,
                       start.degree = 90, cell.padding = c(0, 0, 0, 0),
                       track.margin = c(0.004, 0.008),
                       points.overflow.warning = FALSE)
  circlize::circos.initialize(sectors = df$chr, x = df$pos)

  # Track 1: chromosome band + Mb axis + chromosome label + marker ticks
  circlize::circos.track(
    sectors = df$chr, x = df$pos, ylim = c(0, 1), track.height = 0.07, bg.border = NA,
    panel.fun = function(x, y) {
      xl <- circlize::CELL_META$xlim
      circlize::circos.rect(xl[1], 0.25, xl[2], 0.75, col = "grey70", border = NA)
      circlize::circos.segments(x, rep(0.1, length(x)), x, rep(0.9, length(x)),
                                col = "grey25", lwd = 0.7)
      circlize::circos.text(circlize::CELL_META$xcenter, 2.2,
                            circlize::CELL_META$sector.index,
                            facing = "bending.inside", niceFacing = TRUE, cex = 0.8)
      at <- pretty(xl, n = 5); at <- at[at >= xl[1] & at <= xl[2]]
      circlize::circos.axis(h = "top", major.at = at,
                            labels = paste0(round(at / 1e6, 1)),
                            labels.cex = 0.45, major.tick.length = 0.5, lwd = 0.6)
    })

  # Track 2: # microhaplotypes per marker (bar height = count; paralog-suspect flagged red)
  circlize::circos.track(
    sectors = df$chr, x = df$pos, y = df$n_mhaps, ylim = c(0, max(df$n_mhaps)),
    track.height = 0.22, bg.border = "grey90",
    panel.fun = function(x, y)
      bars(x, y, ifelse(y >= paralog.flag, "#b2182b", "#4575b4")))
  circlize::circos.yaxis("left", sector.index = levels(df$chr)[1], track.index = 2,
                         labels.cex = 0.4, lwd = 0.5)

  # Track 3: mean read depth per marker (bar height + RdBu color = depth)
  circlize::circos.track(
    sectors = df$chr, x = df$pos, y = df$depth, ylim = c(0, max(df$depth)),
    track.height = 0.22, bg.border = "grey90",
    panel.fun = function(x, y) bars(x, y, col_depth(y)))
  circlize::circos.yaxis("left", sector.index = levels(df$chr)[1], track.index = 3,
                         labels.cex = 0.4, lwd = 0.5)

  graphics::title("Marker distribution (circos)", cex.main = 1)
  graphics::legend("bottomleft", inset = c(0, 0), bty = "n", cex = 0.7, border = NA,
                   fill = c("#4575b4", "#b2182b"),
                   legend = c("# mHaps", sprintf("# mHaps >= %d (paralog-suspect)", paralog.flag)))
  graphics::legend("bottomright", inset = c(0, 0), bty = "n", cex = 0.7, border = NA,
                   fill = c("#2166ac", "#f7f7f7", "#b2182b"),
                   legend = c("depth low", "depth mid", "depth high"))
  invisible(NULL)
}

# shared marker data.frame (positions + per-marker #mHaps and mean depth)
#' @keywords internal
#' @noRd
.madc_marker_df <- function(m) {
  if (all(is.na(m$pos)))
    stop("CloneIDs are not in Chr_Pos format; supply `markers_info` with Chr/Pos for the marker/circos plot.")
  df <- data.frame(CloneID = m$markers, chr = m$chr, pos = m$pos,
                   n_mhaps = as.integer(m$n_mhaps_marker),
                   depth = rowMeans(m$depth_total), stringsAsFactors = FALSE)
  df[!is.na(df$pos), , drop = FALSE]
}

# ---- Heatmap ---------------------------------------------------------------
#' @keywords internal
#' @noRd
.madc_plot_heatmap <- function(m, fill = "depth", facet.chrom = FALSE,
                               max.loci = NULL, max.samples = NULL, verbose = TRUE) {
  mat <- if (fill == "mhaps") m$n_mhaps_present else m$depth_total
  fill_lab <- if (fill == "mhaps") "# mHaps" else "Read depth"

  # order loci by chr + pos
  ord <- order(m$chr, m$pos, na.last = TRUE)
  mat <- mat[ord, , drop = FALSE]
  chr_ord <- m$chr[ord]

  # optional subsample for very large panels
  if (!is.null(max.loci) && nrow(mat) > max.loci) {
    idx <- sort(sample.int(nrow(mat), max.loci))
    mat <- mat[idx, , drop = FALSE]; chr_ord <- chr_ord[idx]
    if (verbose) message("Heatmap: subsampled to ", max.loci, " of ",
                         length(ord), " loci.")
  }
  if (!is.null(max.samples) && ncol(mat) > max.samples) {
    idx <- sort(sample.int(ncol(mat), max.samples))
    mat <- mat[, idx, drop = FALSE]
    if (verbose) message("Heatmap: subsampled to ", max.samples, " samples.")
  }

  long <- reshape2::melt(mat, varnames = c("marker", "sample"), value.name = "value")
  long$marker <- factor(long$marker, levels = rownames(mat))
  if (facet.chrom) long$chr <- chr_ord[match(long$marker, rownames(mat))]

  # 0 (failed / no data) -> black; positive values -> a colorblind-safe diverging
  # scale (RdBu) centered at the median so both low- and high-depth markers stand out.
  pos_mid <- stats::median(long$value[long$value > 0], na.rm = TRUE)
  if (!is.finite(pos_mid)) pos_mid <- 1
  long$value[long$value == 0] <- NA

  p <- ggplot2::ggplot(long, ggplot2::aes(x = sample, y = marker, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(fill_lab, low = "#2166ac", mid = "#f7f7f7",
                                  high = "#b2182b", midpoint = pos_mid,
                                  na.value = "black") +
    ggplot2::labs(title = paste0("MADC heatmap (", fill_lab, ")"), x = NULL, y = NULL) +
    .madc_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (nrow(mat) > 60) p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  if (ncol(mat) > 60) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  if (facet.chrom) p <- p + ggplot2::facet_wrap(~ chr, scales = "free_y")
  p
}

# ---- Missing-data boxplot --------------------------------------------------
#' @keywords internal
#' @noRd
.madc_plot_missing <- function(m, grp = NULL) {
  df <- data.frame(sample = m$samples,
                   missing_rate = colMeans(m$missing_mask),
                   stringsAsFactors = FALSE)
  if (!is.null(grp)) df$group <- grp[match(df$sample, m$samples)]

  if (!is.null(grp) && !all(is.na(df$group))) {
    df <- df[!is.na(df$group), , drop = FALSE]
    p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = missing_rate, fill = group)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
      ggplot2::labs(x = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none")
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = "all samples", y = missing_rate)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey80") +
      ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
      ggplot2::labs(x = NULL)
  }
  p + ggplot2::labs(title = "Per-sample missing rate", y = "Missing rate") +
    ggplot2::ylim(0, NA) + .madc_theme()
}
