#' Plot a fixed allele ID MADC
#'
#' Visualize a fixed allele ID MADC file (one processed through HapApp) with a
#' choice of plots selected via `plot.type`: a PCA of alternative read ratios, a
#' linear marker-distribution genome plot, a markers x samples read-depth/#mHaps
#' heatmap, a per-sample missing-data boxplot, an allele read-ratio balance
#' diagnostic, a marker mean-depth (uniformity) distribution, and a `circos`
#' circular plot with concentric marker-position, #mHaps, and read-depth tracks.
#' Requesting one type returns that plot; requesting several (non-circos) returns
#' a multi-panel figure.
#'
#' @details
#' The input is validated with [check_madc_sanity()]; raw DArT MADC files are
#' rejected with an error. All plots honor `target.only` (restrict to `|Ref`/
#' `|Alt`) and `min.depth` (locus x sample cells below this depth are treated as
#' missing). Plotting uses only `ggplot2` and base `grid`; the `circos` plot
#' additionally requires the suggested `circlize` package and, being a base
#' graphics figure, must be requested on its own (it cannot be embedded in a
#' multi-panel). The `"pca"` plot is a sample/data-interpretation diagnostic
#' (missing ratios are mean-imputed and the PCA is unscaled), so strong biological
#' structure can dominate it - read it for sample swaps / population separation,
#' not as evidence the sequencing run itself was good or bad. The `"balance"` plot
#' is the direct read-out of whether the assay yields interpretable dosage signal;
#' it shows every cell with a defined ratio (including below `min.depth`, marked by
#' a dotted line) and, with `ploidy`, a depth-dependent binomial expectation band.
#'
#' @param madc Path to a fixed allele ID MADC file, or an already-read data.frame.
#' @param plot.type One or more of `"pca"`, `"marker"`, `"heatmap"`, `"missing"`,
#'   `"balance"`, `"depth"`, `"circos"`. A single value returns that plot; several
#'   (excluding `"circos"`) return a multi-panel figure. `"circos"` must be
#'   requested on its own.
#' @param target.only Logical; restrict to target `|Ref`/`|Alt` alleles. Default `FALSE`.
#' @param min.depth Minimum locus x sample read depth to be considered present. Default `10`.
#' @param metadata Optional sample-metadata data.frame (see [madc_summary()]) used
#'   to color the PCA / group the missing-data boxplot.
#' @param group.col Category column in `metadata` (e.g. species/population) used
#'   to color the PCA and group the missing-data boxplot / heatmap.
#' @param shape.col Optional second `metadata` column mapped to point shape in
#'   the PCA (e.g. plate or location), so two variables can be shown at once.
#' @param palette Optional vector of colors for the categorical scales (PCA
#'   point color, missing-boxplot fill); interpolated to the number of
#'   categories. `NULL` uses the ggplot2 defaults.
#' @param ploidy Optional integer species ploidy for the `"balance"` plot. When
#'   supplied, dashed guides are drawn at the expected alt-read ratios
#'   `(0:ploidy)/ploidy` (e.g. 0, 0.5, 1 for a diploid) and a depth-dependent 95%
#'   binomial interval band is drawn around each interior expected ratio (the
#'   expected spread of `A/D` under `Binomial(D, ratio)`), so observed dispersion
#'   wider than sampling noise stands out. `NULL` (default) draws the ratio x
#'   depth density with no guides/band, which also suits mixed/unknown ploidy.
#' @param balance.density Logical; for a lone `plot.type = "balance"`, stack a
#'   marginal frequency-polygon of the alt-read ratio above the heatmap (aligned
#'   x-axis). This returns a composite grob instead of a `ggplot`, so it is
#'   ignored when `"balance"` is combined with other plot types. Default `FALSE`.
#' @param markers_info Optional marker lookup (id + `Chr` + `Pos`) for the marker/circos plot.
#' @param pc.x,pc.y Principal components to plot (PCA). Defaults `1`, `2`.
#' @param loci.miss.max,sample.miss.max Drop loci/samples with missingness above
#'   these fractions before PCA. Defaults `0.5`.
#' @param miss.sort Ordering of the category boxes in the `"missing"` plot by their
#'   median missing rate: `"none"` (default, keeps category order), `"desc"`
#'   (highest to lowest), or `"asc"` (lowest to highest). Only applies when a
#'   grouping category is supplied.
#' @param horizontal Logical; draw the `"missing"` boxplots horizontally
#'   (`coord_flip`). Default `FALSE`. With `miss.sort`, the order then reads
#'   top-to-bottom.
#' @param fill Heatmap fill metric: `"depth"` (default) or `"mhaps"`.
#' @param facet.chrom Logical; facet the heatmap by chromosome. Default `FALSE`.
#'   When `metadata`/`group.col` are supplied, the heatmap also sorts samples by
#'   category and labels each category (a facet column per group).
#' @param max.loci,max.samples Optional caps; larger heatmaps are subsampled (with a note).
#' @param density.window Window size in bp for the circos marker-density track.
#'   Default `1e6` (1 Mb).
#' @param density.col Gradient colors (low -> high) for the circos marker-density
#'   heatmap ring; interpolated across the range. Default a grey-to-blue gradient.
#' @param depth.col Length-3 vector of colors for the circos "Depth QC" ribbon:
#'   markers that are `low depth` (mean depth below `min.depth`), `normal`, or a
#'   high-depth `outlier`. Only low-depth and outlier markers are highlighted;
#'   normal markers take the neutral middle color. Default
#'   `c("#2166ac","#e0e0e0","#b2182b")`.
#' @param depth.outlier Multiplier `k` for the robust high-depth outlier cutoff
#'   `median + k * MAD` (not mean +/- SD, which a few extreme loci would drag up).
#'   Default `3`.
#' @param depth.qc.maxmiss Circos "Depth" ring: a marker is flagged problematic
#'   when its mean depth is below `min.depth` OR its missing rate exceeds this
#'   fraction - catching high-dropout markers (e.g. 50% zeros / 50% deep) that a
#'   mean-depth cut alone reads as normal. Default `0.5`.
#' @param mhap.col Fill color for the circos "mHaps" log-scaled radial bars.
#'   Default `"#35978f"`.
#' @param paralog.flag Microhaplotype count above which a circos "mHaps" locus is
#'   flagged as paralog-suspect: its bar's over-threshold excess is drawn as a red
#'   segment above a dashed reference ring. Default `10`.
#' @param tick.height,density.height,mhap.height Circos track heights (fraction of
#'   the radius) for the per-marker locus tick rug, the density heatmap ring, and
#'   the mHap bar track. Defaults `0.03`, `0.08`, and `0.18`.
#' @param label.gap Width in degrees of the label corridor opened at 12 o'clock
#'   for the circos track titles. Default `12`.
#' @param mhap.min.reads Minimum reads for a mHap to count as present. Default `1`.
#' @param output.file If not `NULL`, the figure is also saved here (`.png`/`.pdf`).
#' @param width,height,dpi Saved-figure dimensions (inches) and resolution.
#' @param verbose Logical; print progress/validation messages. Default `TRUE`.
#'
#' @return For a single `plot.type`, a `ggplot` (a composite grob when
#'   `plot.type = "balance"` and `balance.density = TRUE`); for several, an
#'   (invisible) list with `plots` (named list) and `panel` (the assembled
#'   multi-panel grob). `"circos"` draws to the active device and returns `NULL`
#'   invisibly.
#'
#' @examples
#' madc_file <- system.file("example_MADC_FixedAlleleID.csv", package = "BIGr")
#' p <- madc_plot(madc_file, plot.type = "missing")
#' \donttest{
#' madc_plot(madc_file, plot.type = c("pca", "heatmap"))
#' madc_plot(madc_file, plot.type = "balance", ploidy = 2)
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
                      shape.col      = NULL,
                      palette        = NULL,
                      ploidy         = NULL,
                      balance.density = FALSE,
                      markers_info   = NULL,
                      pc.x           = 1,
                      pc.y           = 2,
                      loci.miss.max  = 0.5,
                      sample.miss.max = 0.5,
                      miss.sort      = c("none", "desc", "asc"),
                      horizontal     = FALSE,
                      fill           = c("depth", "mhaps"),
                      facet.chrom    = FALSE,
                      max.loci       = NULL,
                      max.samples    = NULL,
                      density.window = 1e6,
                      density.col    = c("grey92", "#4292c6", "#08306b"),
                      depth.col      = c("#2166ac", "#e0e0e0", "#b2182b"),
                      depth.outlier  = 3,
                      depth.qc.maxmiss = 0.5,
                      mhap.col       = "#35978f",
                      paralog.flag   = 10,
                      tick.height    = 0.03,
                      density.height = 0.08,
                      mhap.height    = 0.18,
                      label.gap      = 12,
                      mhap.min.reads = 1,
                      output.file    = NULL,
                      width = 8, height = 6, dpi = 300,
                      verbose        = TRUE) {

  plot.type <- match.arg(plot.type, c("pca", "marker", "heatmap", "missing",
                                      "balance", "depth", "circos"),
                         several.ok = TRUE)
  fill      <- match.arg(fill)
  miss.sort <- match.arg(miss.sort)

  report <- .read_and_check_madc(madc, verbose = verbose)
  m <- .madc_metrics(report, min.depth = min.depth, target.only = target.only,
                     mhap.min.reads = mhap.min.reads, markers_info = markers_info,
                     verbose = FALSE)
  grp <- .madc_group(metadata, group.col, m$samples)
  shp <- if (!is.null(shape.col)) .madc_group(metadata, shape.col, m$samples) else NULL

  # circos is a standalone base-graphics (circlize) figure - it cannot be embedded
  # in a grid multi-panel, so it must be requested on its own.
  if ("circos" %in% plot.type) {
    if (length(plot.type) > 1)
      stop("plot.type = 'circos' must be requested on its own (it renders as a standalone circular figure).")
    save_to <- output.file
    # no output.file and no graphics device open: capture to a temp file rather
    # than letting base graphics leave a stray Rplots.pdf. When a device is already
    # open (interactive session, knitr/R Markdown, or a user-opened device) draw
    # straight to it so the circos renders in place.
    if (is.null(save_to) && grDevices::dev.cur() == 1L) {
      save_to <- tempfile(fileext = ".png")
      if (verbose) message("No graphics device / output.file; circos written to ", save_to)
    }
    if (!is.null(save_to)) {
      ext <- tolower(tools::file_ext(save_to))
      if (ext == "pdf") grDevices::pdf(save_to, width = width, height = height)
      else grDevices::png(save_to, width = width, height = height, units = "in", res = dpi)
      on.exit(grDevices::dev.off(), add = TRUE)
    }
    .madc_plot_circos(m, paralog.flag = paralog.flag, density.window = density.window,
                      density.col = density.col, depth.col = depth.col,
                      depth.outlier = depth.outlier, depth.qc.maxmiss = depth.qc.maxmiss,
                      mhap.col = mhap.col,
                      tick.height = tick.height, density.height = density.height,
                      mhap.height = mhap.height, label.gap = label.gap)
    return(invisible(save_to))
  }

  # the marginal-density composite only applies to a lone "balance" plot (it
  # returns a grob, which can't embed in the ggplot multi-panel)
  bal_density <- isTRUE(balance.density) && length(plot.type) == 1L

  build <- list(
    pca     = function() .madc_plot_pca(m, pc.x, pc.y, loci.miss.max, sample.miss.max, grp, shp, palette),
    marker  = function() .madc_plot_marker(m, facet.chrom = facet.chrom),
    heatmap = function() .madc_plot_heatmap(m, fill, facet.chrom, max.loci, max.samples, verbose, grp),
    missing = function() .madc_plot_missing(m, grp, palette, sort = miss.sort,
                                            horizontal = horizontal),
    balance = function() .madc_plot_balance(m, ploidy = ploidy, density = bal_density, palette = palette),
    depth   = function() .madc_plot_depth(m)
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
                           loci.miss.max = 0.5, sample.miss.max = 0.5,
                           grp = NULL, shp = NULL, palette = NULL) {
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
  has_grp <- !is.null(grp)
  has_shp <- !is.null(shp)
  if (has_grp) df$group <- grp[match(df$sample, m$samples)]
  if (has_shp) df$shape <- shp[match(df$sample, m$samples)]

  aes_args <- list(x = quote(PCx), y = quote(PCy))
  if (has_grp) aes_args$colour <- quote(group)
  if (has_shp) aes_args$shape  <- quote(shape)

  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = 2.6, alpha = 0.9) +
    ggplot2::labs(title = "PCA of alternative read ratios",
                  x = sprintf("PC%d (%.1f%%)", pc.x, 100 * ve[pc.x]),
                  y = sprintf("PC%d (%.1f%%)", pc.y, 100 * ve[pc.y])) +
    .madc_theme()

  if (has_grp && !is.null(palette)) {
    n <- length(unique(stats::na.omit(df$group)))
    p <- p + ggplot2::scale_colour_manual(values = grDevices::colorRampPalette(palette)(n))
  }
  if (has_shp) {
    shapes <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0, 5, 6, 7, 9, 10, 11, 12, 13, 14)
    nsh <- length(unique(stats::na.omit(df$shape)))
    p <- p + ggplot2::scale_shape_manual(values = rep(shapes, length.out = nsh))
  }
  p
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
.madc_plot_circos <- function(m, paralog.flag = 10, density.window = 1e6,
                              density.col = c("grey92", "#4292c6", "#08306b"),
                              depth.col = c("#2166ac", "#e0e0e0", "#b2182b"),
                              depth.outlier = 3, depth.qc.maxmiss = 0.5,
                              mhap.col = "#35978f",
                              tick.height = 0.03, density.height = 0.08,
                              mhap.height = 0.18, label.gap = 12) {
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("plot.type = 'circos' requires the 'circlize' package. Install it with install.packages('circlize').")

  df <- .madc_marker_df(m)
  df$chr <- factor(df$chr, levels = sort(unique(df$chr)))
  df <- df[order(df$chr, df$pos), , drop = FALSE]

  # per-window marker density (markers per `density.window` bp) per chromosome
  dens <- lapply(levels(df$chr), function(cc) {
    p  <- df$pos[df$chr == cc]
    br <- seq(0, max(p) + density.window, by = density.window)
    h  <- graphics::hist(p, breaks = br, plot = FALSE)
    data.frame(start = utils::head(br, -1), end = br[-1], count = h$counts)
  })
  names(dens) <- levels(df$chr)
  dmax <- max(1, vapply(dens, function(d) max(d$count), numeric(1)))
  col_dens <- circlize::colorRamp2(seq(0, dmax, length.out = length(density.col)), density.col)

  # mHap bar half-width scaled to marker spacing (visible but not overlapping)
  bar_hw <- max(vapply(split(df$pos, df$chr), function(p)
    if (length(p) > 1) stats::median(diff(sort(p))) * 0.35 else diff(range(df$pos)) * 0.01,
    numeric(1)))

  # Combined depth QC. The outlier cut is the robust median + k*MAD (not mean +/-
  # SD, which a few extreme loci would drag up). A marker is "problematic" if its
  # mean depth is below min.depth OR it is missing (below min.depth) in more than
  # depth.qc.maxmiss of samples - the high-dropout case (e.g. 50% zeros / 50% deep)
  # that a mean-depth cut alone reads as "normal". Colors are precomputed per
  # marker and looked up inside the track by a chr_pos key.
  depth_thr <- stats::median(df$depth, na.rm = TRUE) +
    depth.outlier * stats::mad(df$depth, na.rm = TRUE)
  if (!is.finite(depth_thr) || depth_thr <= 0) depth_thr <- max(df$depth, na.rm = TRUE)
  miss_thr <- m$min.depth
  qc_problem <- (df$depth < miss_thr) | ((1 - df$call_rate) > depth.qc.maxmiss)
  qc_col <- rep(depth.col[2], nrow(df))                        # normal
  qc_col[qc_problem] <- depth.col[1]                           # low depth / high missing
  qc_col[df$depth >= depth_thr & !qc_problem] <- depth.col[3]  # high-depth outlier
  names(qc_col) <- paste(df$chr, df$pos, sep = "_")

  op <- graphics::par(mar = c(1, 1, 2, 1)); on.exit(graphics::par(op), add = TRUE)
  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  N <- nlevels(df$chr)
  top_gap <- if (N == 1) max(label.gap, 30) else label.gap  # roomier corridor when one sector wraps
  gaps <- if (N == 1) top_gap else c(rep(3, N - 1), top_gap)
  circlize::circos.par(gap.after = gaps, start.degree = 90 - top_gap / 2,  # center gap at 12 o'clock
                       cell.padding = c(0, 0, 0, 0), track.margin = c(0.006, 0.006),
                       canvas.xlim = c(-1.35, 1.35), canvas.ylim = c(-1.35, 1.35),
                       points.overflow.warning = FALSE)
  circlize::circos.initialize(sectors = df$chr, x = df$pos)

  chr_lab_y <- 1 + 0.18 / tick.height

  # Track 1 (Loci): per-marker position ticks (kept dominant) + Mb coordinate
  # axis (longer/heavier/darker so it stays distinct) + chromosome label
  circlize::circos.track(
    sectors = df$chr, x = df$pos, ylim = c(0, 1), track.height = tick.height, bg.border = NA,
    panel.fun = function(x, y) {
      xl <- circlize::CELL_META$xlim
      circlize::circos.segments(x, rep(0, length(x)), x, rep(1, length(x)),
                                col = "grey25", lwd = 0.4)
      circlize::circos.text(circlize::CELL_META$xcenter, chr_lab_y,
                            circlize::CELL_META$sector.index,
                            facing = "bending.inside", niceFacing = TRUE, cex = 0.8)
      at <- pretty(xl, n = 5); at <- at[at >= xl[1] & at <= xl[2]]
      circlize::circos.axis(h = "top", major.at = at, labels = paste0(round(at / 1e6, 1)),
                            labels.cex = 0.45, major.tick.length = circlize::mm_y(1.5),
                            minor.ticks = 0, lwd = 0.9, col = "black")
    })

  # Track 2 (Density): marker density heatmap over a grey background
  circlize::circos.track(
    sectors = df$chr, ylim = c(0, 1), track.height = density.height,
    bg.col = "grey92", bg.border = "grey80",
    panel.fun = function(x, y) {
      d <- dens[[circlize::CELL_META$sector.index]]
      circlize::circos.rect(d$start, 0, d$end, 1, col = col_dens(d$count), border = NA)
    })

  # Track 3 (mHaps): log2 bars (real-count axis) as two-tone bars - a neutral base
  # up to the paralog threshold and a red segment for the over-threshold excess -
  # plus a dashed reference ring at the threshold
  ymax_l <- log2(max(df$n_mhaps, na.rm = TRUE) + 1)
  ref_l  <- log2(paralog.flag + 1)
  circlize::circos.track(
    sectors = df$chr, x = df$pos, y = log2(df$n_mhaps + 1), ylim = c(0, ymax_l),
    track.height = mhap.height, bg.border = "grey90",
    panel.fun = function(x, y) {
      base <- pmin(y, ref_l)
      circlize::circos.rect(x - bar_hw, 0, x + bar_hw, base, col = mhap.col, border = NA)
      hi <- y > ref_l
      if (any(hi))
        circlize::circos.rect(x[hi] - bar_hw, ref_l, x[hi] + bar_hw, y[hi],
                              col = "#b2182b", border = NA)
      if (ref_l < ymax_l)
        circlize::circos.lines(circlize::CELL_META$xlim, c(ref_l, ref_l),
                               col = "grey45", lty = 2, lwd = 0.6)
    })
  # radial scale on a sector near the bottom, clear of the 12 o'clock corridor
  mlabs <- c(0, 2, 5, 10, 25, 50, 100); mk <- log2(mlabs + 1) <= ymax_l
  ax_sector <- levels(df$chr)[1L + (N %/% 2L)]
  circlize::circos.yaxis("left", at = log2(mlabs[mk] + 1), labels = mlabs[mk],
                         sector.index = ax_sector, track.index = 3,
                         labels.cex = 0.4, lwd = 0.5)

  # Track 4 (Depth): contiguous categorical tile ribbon (gapless midpoint edges)
  circlize::circos.track(
    sectors = df$chr, x = df$pos, ylim = c(0, 1),
    track.height = 0.055, bg.border = "grey90",
    panel.fun = function(x, y) {
      xs  <- sort(x)
      mid <- (utils::head(xs, -1) + xs[-1]) / 2
      key <- paste(circlize::CELL_META$sector.index, xs, sep = "_")
      circlize::circos.rect(c(circlize::CELL_META$xlim[1], mid), 0,
                            c(mid, circlize::CELL_META$xlim[2]), 1,
                            col = qc_col[key], border = NA)
    })

  # horizontal track titles in the 12 o'clock corridor, aligned to each ring
  titles <- c("Loci", "Density", "mHaps", "Depth")
  s1 <- levels(df$chr)[1]
  for (ti in seq_along(titles)) {
    rt <- circlize::get.cell.meta.data("cell.top.radius",    s1, track.index = ti)
    rb <- circlize::get.cell.meta.data("cell.bottom.radius", s1, track.index = ti)
    graphics::text(0, (rt + rb) / 2, titles[ti], adj = c(0.5, 0.5),
                   cex = 0.5, font = 2, xpd = NA)
  }

  # compact color/scale keys in the empty center
  .circos_center_legend(col_dens, dmax, density.window, mhap.col, paralog.flag,
                        depth.col, miss_thr, depth_thr, depth.qc.maxmiss)

  graphics::title("Marker distribution", cex.main = 1.1)
  invisible(NULL)
}

# compact legend drawn in the center hole of the circos (base-graphics, circos coords)
#' @keywords internal
#' @noRd
.circos_center_legend <- function(col_dens, dmax, density.window, mhap.col, paralog.flag,
                                  depth.col, miss_thr, depth_thr, depth.qc.maxmiss = 0.5) {
  # density: small horizontal continuous colour bar
  n <- 64; bx <- seq(-0.17, 0.17, length.out = n + 1); yb <- 0.15; yt <- 0.19
  vals <- seq(0, dmax, length.out = n)
  for (i in seq_len(n))
    graphics::rect(bx[i], yb, bx[i + 1], yt, col = col_dens(vals[i]), border = NA, xpd = NA)
  graphics::rect(-0.17, yb, 0.17, yt, border = "grey55", lwd = 0.5, xpd = NA)
  at <- pretty(c(0, dmax), n = 4); at <- at[at >= 0 & at <= dmax]
  graphics::text(-0.17 + 0.34 * at / dmax, yb - 0.028, labels = at, cex = 0.42, xpd = NA)
  graphics::text(0, yt + 0.035, sprintf("Markers/%g Mb", density.window / 1e6),
                 cex = 0.5, font = 2, xpd = NA)

  # mHaps key
  graphics::text(0, 0.075, "mHaps", cex = 0.5, font = 2, xpd = NA)
  graphics::rect(-0.16, 0.03, -0.13, 0.055, col = mhap.col, border = NA, xpd = NA)
  graphics::text(-0.115, 0.0425, "count (log2)", cex = 0.42, adj = 0, xpd = NA)
  graphics::rect(0.03, 0.03, 0.06, 0.055, col = "#b2182b", border = NA, xpd = NA)
  graphics::text(0.07, 0.0425, sprintf("> %d (paralog)", paralog.flag), cex = 0.42, adj = 0, xpd = NA)

  # depth key
  graphics::text(0, -0.03, "Depth", cex = 0.5, font = 2, xpd = NA)
  dl <- c(sprintf("low <%g / miss >%g", miss_thr, depth.qc.maxmiss), "normal",
          sprintf("outlier (>= %g)", round(depth_thr)))
  yk <- -0.065 - c(0, 0.035, 0.07)
  for (j in 1:3) {
    graphics::rect(-0.16, yk[j] - 0.011, -0.13, yk[j] + 0.011,
                   col = depth.col[j], border = "grey70", xpd = NA)
    graphics::text(-0.115, yk[j], dl[j], cex = 0.42, adj = 0, xpd = NA)
  }
}

# shared marker data.frame (positions + per-marker #mHaps and mean depth)
#' @keywords internal
#' @noRd
.madc_marker_df <- function(m) {
  if (all(is.na(m$pos)))
    stop("CloneIDs are not in Chr_Pos format; supply `markers_info` with Chr/Pos for the marker/circos plot.")
  df <- data.frame(CloneID = m$markers, chr = m$chr, pos = m$pos,
                   n_mhaps = as.integer(m$n_mhaps_marker),
                   depth = rowMeans(m$depth_total),
                   call_rate = 1 - rowMeans(m$missing_mask),
                   stringsAsFactors = FALSE)
  df[!is.na(df$pos), , drop = FALSE]
}

# ---- Heatmap ---------------------------------------------------------------
#' @keywords internal
#' @noRd
.madc_plot_heatmap <- function(m, fill = "depth", facet.chrom = FALSE,
                               max.loci = NULL, max.samples = NULL, verbose = TRUE,
                               grp = NULL) {
  mat <- if (fill == "mhaps") m$n_mhaps_present else m$depth_total
  fill_lab <- if (fill == "mhaps") "# mHaps" else "Read depth"

  # order loci by chr + pos
  ord <- order(m$chr, m$pos, na.last = TRUE)
  mat <- mat[ord, , drop = FALSE]
  chr_ord <- m$chr[ord]

  # optional subsample for very large panels - deterministic evenly-spaced pick so
  # repeated calls give identical figures and genomic/sample spread is preserved
  if (!is.null(max.loci) && nrow(mat) > max.loci) {
    idx <- unique(round(seq(1, nrow(mat), length.out = max.loci)))
    mat <- mat[idx, , drop = FALSE]; chr_ord <- chr_ord[idx]
    if (verbose) message("Heatmap: subsampled to ", length(idx), " of ",
                         length(ord), " loci (evenly spaced).")
  }
  if (!is.null(max.samples) && ncol(mat) > max.samples) {
    idx <- unique(round(seq(1, ncol(mat), length.out = max.samples)))
    mat <- mat[, idx, drop = FALSE]
    if (verbose) message("Heatmap: subsampled to ", length(idx), " samples (evenly spaced).")
  }

  long <- reshape2::melt(mat, varnames = c("marker", "sample"), value.name = "value")
  long$marker <- factor(long$marker, levels = rownames(mat))
  if (facet.chrom) long$chr <- chr_ord[match(long$marker, rownames(mat))]

  # when a category is given, sort samples by category and add ONE x-axis label
  # per category (at its centre) so category labels never overlap, with white
  # separators between categories.
  has_grp <- !is.null(grp) && !all(is.na(grp))
  cat_breaks <- cat_labels <- boundaries <- NULL
  if (has_grp) {
    sgrp <- grp[match(as.character(long$sample), m$samples)]
    keep <- !is.na(sgrp)
    long <- long[keep, , drop = FALSE]; sgrp <- sgrp[keep]
    samp_order <- unique(long$sample[order(sgrp, as.character(long$sample))])
    long$sample <- factor(long$sample, levels = samp_order)
    og   <- grp[match(as.character(samp_order), m$samples)]   # group per ordered sample
    runs <- rle(og)
    ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1L
    cat_breaks <- as.character(samp_order)[round((starts + ends) / 2)]
    cat_labels <- runs$values
    boundaries <- utils::head(ends, -1) + 0.5
  }

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
    .madc_theme()
  if (nrow(mat) > 60) p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())

  if (has_grp) {
    p <- p +
      ggplot2::geom_vline(xintercept = boundaries, color = "white", linewidth = 0.4) +
      ggplot2::scale_x_discrete(breaks = cat_breaks, labels = cat_labels) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                     axis.ticks.x = ggplot2::element_blank())
  } else if (ncol(mat) > 60) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank())
  } else {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  }

  # facet by chromosome (rows), keeping the chromosome labels visible on the left
  if (facet.chrom)
    p <- p + ggplot2::facet_grid(chr ~ ., scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0),
                     strip.placement = "outside", panel.spacing.y = grid::unit(2, "pt"))
  p
}

# ---- Missing-data boxplot --------------------------------------------------
#' @keywords internal
#' @noRd
.madc_plot_missing <- function(m, grp = NULL, palette = NULL,
                               sort = "none", horizontal = FALSE) {
  df <- data.frame(sample = m$samples,
                   missing_rate = colMeans(m$missing_mask),
                   stringsAsFactors = FALSE)
  if (!is.null(grp)) df$group <- grp[match(df$sample, m$samples)]

  grouped <- !is.null(grp) && !all(is.na(df$group))
  if (grouped) {
    df <- df[!is.na(df$group), , drop = FALSE]
    # order the category boxes by their median missing rate when requested; in
    # horizontal mode reverse so the sort reads top-to-bottom (coord_flip stacks
    # the first factor level at the bottom)
    if (sort != "none") {
      med <- tapply(df$missing_rate, df$group, stats::median, na.rm = TRUE)
      ord <- names(base::sort(med, decreasing = (sort == "desc")))
      if (horizontal) ord <- rev(ord)
      df$group <- factor(df$group, levels = ord)
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = missing_rate, fill = group)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
      ggplot2::labs(x = NULL) +
      ggplot2::theme(legend.position = "none")
    # angled tick labels only help the vertical layout; horizontal reads flat
    if (!horizontal)
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    if (!is.null(palette)) {
      n <- length(unique(df$group))
      p <- p + ggplot2::scale_fill_manual(values = grDevices::colorRampPalette(palette)(n))
    }
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = "all samples", y = missing_rate)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey80") +
      ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
      ggplot2::labs(x = NULL)
  }
  p <- p + ggplot2::labs(title = "Per-sample missing rate", y = "Missing rate") +
    ggplot2::ylim(0, NA) + .madc_theme()
  if (horizontal) p <- p + ggplot2::coord_flip()
  p
}

# ---- Allele read-ratio balance --------------------------------------------
#' @keywords internal
#' @noRd
.madc_plot_balance <- function(m, ploidy = NULL, density = FALSE, palette = NULL) {
  # One point-density cell per (marker x sample) with a defined allele ratio: alt
  # read ratio on x, target depth (ref + alt) on a log y. Shows whether the
  # chemistry yields interpretable dosage signal (bands), whether balance depends
  # on depth, and - with `ploidy` - whether observed spread exceeds binomial
  # sampling. `palette` is accepted for signature symmetry but unused (the fill is
  # a continuous cell count).
  ratio <- m$alt_ratio
  size  <- m$size_matrix
  keep  <- is.finite(ratio) & is.finite(size) & size >= 1   # keep sub-min.depth cells too
  if (!any(keep))
    stop("No (marker x sample) cells have a defined allele ratio for the balance plot.")
  df <- data.frame(alt_ratio = as.numeric(ratio[keep]),
                   depth      = as.numeric(size[keep]))

  gx <- NULL
  if (!is.null(ploidy)) {
    ploidy <- as.integer(ploidy)
    if (is.na(ploidy) || ploidy < 1L) stop("`ploidy` must be a positive integer.")
    gx <- (0:ploidy) / ploidy
  }

  # shared x scale so the optional marginal panel lines up with the heatmap
  x_scale <- ggplot2::scale_x_continuous(limits = c(-0.02, 1.02),
                                         breaks = seq(0, 1, 0.25), expand = c(0, 0))
  guides  <- if (!is.null(gx))
    ggplot2::geom_vline(xintercept = gx, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4, alpha = 0.6) else NULL
  comma <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

  # depth-dependent 95% binomial interval of A/D around each interior expected
  # ratio (funnels toward the ratio as depth rises); homozygous 0/1 are degenerate
  env <- if (!is.null(gx)) .balance_binom_envelope(gx, max(df$depth)) else NULL

  # --- heatmap: count in BOTH fill and alpha (log) so sparse noise recedes ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = alt_ratio, y = depth)) +
    ggplot2::geom_bin2d(ggplot2::aes(alpha = ggplot2::after_stat(count)), bins = 60) +
    ggplot2::scale_fill_viridis_c("Observations per bin", trans = "log10",
                                  breaks = c(1, 10, 100, 1000, 10000), labels = comma) +
    ggplot2::scale_alpha_continuous(trans = "log10", range = c(0.3, 1), guide = "none") +
    x_scale + ggplot2::scale_y_log10() +
    ggplot2::geom_hline(yintercept = m$min.depth, linetype = "dotted",
                        linewidth = 0.5, colour = "grey30")
  if (!is.null(env))
    p <- p +
      ggplot2::geom_ribbon(data = env, inherit.aes = FALSE,
        ggplot2::aes(y = depth, xmin = ratio_lo, xmax = ratio_hi, group = grp),
        orientation = "y", fill = "grey30", alpha = 0.10) +
      ggplot2::geom_path(data = env, inherit.aes = FALSE,
        ggplot2::aes(x = ratio_lo, y = depth, group = grp),
        linetype = "dotted", linewidth = 0.35, colour = "grey20") +
      ggplot2::geom_path(data = env, inherit.aes = FALSE,
        ggplot2::aes(x = ratio_hi, y = depth, group = grp),
        linetype = "dotted", linewidth = 0.35, colour = "grey20")
  if (!is.null(guides)) p <- p + guides

  sub <- if (!is.null(gx))
    paste0(sprintf("Expected ratios (ploidy %d): %s", ploidy,
                   paste(round(gx, 3), collapse = ", ")),
           if (!is.null(env)) "; band = 95% binomial interval" else "") else NULL
  p <- p +
    ggplot2::labs(title = "Allele balance vs. target depth", subtitle = sub,
                  x = "ALT read ratio (ALT / [REF + ALT])",
                  y = "Target depth (REF + ALT reads, log scale)",
                  caption = sprintf("Dotted line = min.depth (%g); cells below would be dropped at that threshold",
                                    m$min.depth)) +
    .madc_theme()

  if (!isTRUE(density)) return(p)

  # --- optional marginal observed-ratio density stacked above (aligned x) ---
  # histogram (not KDE) because allele ratios are discrete at low depth. Use a
  # data-range scale + coord_cartesian (zoom, not hard limits) so the boundary
  # bins aren't clipped/warned; the -0.02..1.02 view matches the heatmap panel.
  top <- ggplot2::ggplot(df, ggplot2::aes(x = alt_ratio)) +
    ggplot2::geom_histogram(bins = 80, fill = "grey35", colour = NA) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(-0.02, 1.02)) +
    ggplot2::labs(title = "Allele balance vs. target depth", subtitle = sub, y = "Obs.") +
    .madc_theme() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  if (!is.null(guides)) top <- top + guides
  p <- p + ggplot2::labs(title = NULL, subtitle = NULL)   # one title, on the top panel

  .arrange_grobs(list(top, p), nrow = 2, ncol = 1,
                 heights = grid::unit(c(0.25, 0.75), "null"), align_widths = TRUE)
}

# depth-stratified 95% binomial interval of A/D around each interior expected ratio
#' @keywords internal
#' @noRd
.balance_binom_envelope <- function(gx, max_depth, n = 200, lo = 0.025, hi = 0.975) {
  interior <- gx[gx > 0 & gx < 1]
  if (!length(interior)) return(NULL)
  Dgrid <- unique(round(10^seq(0, log10(max(max_depth, 2)), length.out = n)))
  Dgrid <- Dgrid[Dgrid >= 1]
  do.call(rbind, lapply(seq_along(interior), function(i) {
    p <- interior[i]
    data.frame(depth    = Dgrid,
               ratio_lo = stats::qbinom(lo, Dgrid, p) / Dgrid,
               ratio_hi = stats::qbinom(hi, Dgrid, p) / Dgrid,
               grp      = paste0("p", i),
               stringsAsFactors = FALSE)
  }))
}

# ---- Marker depth (uniformity) distribution --------------------------------
#' @keywords internal
#' @noRd
.madc_plot_depth <- function(m) {
  d   <- rowMeans(m$depth_total)
  d   <- d[is.finite(d)]
  med <- stats::median(d)
  q1  <- stats::quantile(d, 0.25, names = FALSE)
  q3  <- stats::quantile(d, 0.75, names = FALSE)
  nz  <- sum(d <= 0)
  df  <- data.frame(depth = d[d > 0])
  if (!nrow(df)) stop("All markers have zero mean depth; nothing to plot.")

  sub <- sprintf("median %.1f (IQR %.1f-%.1f); dashed line = min.depth %g",
                 med, q1, q3, m$min.depth)
  if (nz > 0) sub <- paste0(sub, sprintf("; %d zero-depth marker(s) omitted", nz))

  ggplot2::ggplot(df, ggplot2::aes(x = depth)) +
    ggplot2::geom_histogram(bins = 40, fill = "#4292c6", colour = "white", linewidth = 0.2) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = m$min.depth, linetype = 2,
                        colour = "#b2182b", linewidth = 0.5) +
    ggplot2::labs(title = "Marker mean-depth distribution", subtitle = sub,
                  x = "Mean read depth per marker (log10)", y = "Markers") +
    .madc_theme()
}
