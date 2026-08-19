#' Convert rhAmpSeq Haplotype Genotypes to MADC Format
#'
#' Reads a rhAmpSeq haplotype genotype file and the matching FASTA of haplotype
#' allele sequences, then assembles a Missing Allele Discovery Counts (MADC)
#' data frame. For each locus, all alleles are globally aligned to the
#' \code{#1} reference sequence (or the allele with the smallest available
#' numeric ID if \code{#1} is absent). A synthetic \code{Alt} allele is
#' introduced at the reference position with the fewest surrounding
#' polymorphisms, ensuring that no existing allele coincides with it — making
#' all original alleles \code{RefMatch}. Because the synthetic \code{Alt} is
#' assigned ID \code{0002}, every \code{RefMatch} allele receives its original
#' rhAmpSeq numeric ID (from \code{#N} in the FASTA) incremented by one.
#'
#' Alleles with sequences identical to the reference (\code{#1}) are detected
#' and their read depths are merged into the reference row; a warning is
#' reported via \code{vmsg()}. Loci absent from the FASTA are silently skipped
#' (with a \code{vmsg()} message) and excluded from all outputs. When every
#' reference position is covered by at least one polymorphism, the synthetic
#' SNP is placed at the centre of the longest run of positions with the minimum
#' polymorphism count across alleles.
#'
#' Homozygous genotype cells with a single depth value (e.g. \code{"2/2:39"})
#' are collapsed to a single allele entry (\code{"2:39"}) before depth
#' aggregation, so the reported depth is attributed once to the unique allele.
#' Missing calls (\code{"./.:0"}) are excluded from depth aggregation; any
#' allele absent from all samples receives a depth of zero.
#'
#' @param hap_genotype_file Path to a tab-delimited haplotype genotype file.
#'   Column 1 must contain locus names; column 2 must contain haplotype
#'   frequencies in \code{"N(freq);"} format (informational, not validated);
#'   columns 3 onward must contain per-sample genotypes in
#'   \code{"allele1/allele2:depth1,depth2"} format (\code{"./.:0"} for missing
#'   calls). Homozygous calls with a single depth (e.g. \code{"2/2:39"}) are
#'   collapsed to a single allele/depth pair before aggregation.
#' @param haplotype_allele_fasta Path to a FASTA file of haplotype allele
#'   sequences. Sequence names must follow the \code{"LocusID#N"} convention
#'   (e.g. \code{"rhMAS_5GT_cons95#1"}). Loci present in
#'   \code{hap_genotype_file} but absent from this file are skipped.
#' @param n_cores Number of cores for parallel locus processing. Defaults to
#'   \code{1} (sequential). On Windows a PSOCK cluster is always used; on
#'   Unix/macOS the type is controlled by \code{parallel_type}.
#' @param parallel_type Character string controlling the cluster type on
#'   Unix/macOS when \code{n_cores > 1}. Use \code{"PSOCK"} for socket-based
#'   workers (safer, slightly more overhead) or \code{NULL} (default) for
#'   \code{"FORK"}-based workers (faster, shares memory). Ignored on Windows,
#'   which always uses PSOCK.
#' @param prefix Optional character string. When provided, three files are
#'   written: \code{<prefix>_madc.csv}, \code{<prefix>_target_positions.csv},
#'   and \code{<prefix>.fasta}. Defaults to \code{NULL} (no files written).
#' @param verbose Logical. When \code{TRUE}, progress and per-locus warning
#'   messages are printed via \code{vmsg()}. Defaults to \code{TRUE}.
#'
#' @return A named \code{list} with four elements:
#'   \describe{
#'     \item{madc}{A \code{data.frame} in MADC format with columns
#'       \code{AlleleID}, \code{CloneID}, \code{AlleleSequence}, and one
#'       column per sample containing integer read depths. Only loci found in
#'       the FASTA are included. Sample depth columns contain no \code{NA}
#'       values; unobserved alleles and missing calls are represented as
#'       \code{0}.}
#'     \item{target_positions}{A \code{data.frame} with one row per processed
#'       locus and columns \code{AlleleID_rhAmpSeq} (reference sequence name,
#'       e.g. \code{"LocusID#1"}), \code{AlleleID} (locus name),
#'       \code{target_position} (integer position of the synthetic SNP within
#'       the reference sequence), and \code{target_base} (substitution in
#'       \code{"from/to"} format, e.g. \code{"A/C"}).}
#'     \item{new_fasta}{A \code{DNAStringSet} containing all allele sequences
#'       from \code{madc}, with \code{AlleleID} values as sequence names.}
#'     \item{info}{A named \code{list} of diagnostic character vectors:
#'       \code{no_fasta_ids} (locus names absent from the FASTA),
#'       \code{fallback_ref_ids} (locus names for which \code{#1} was absent
#'       and a fallback reference was used),
#'       \code{dup_ref_ids} (full FASTA sequence IDs — including the
#'       \code{#N} suffix — of alleles found to be identical to the reference
#'       and whose depths were merged into it, e.g. \code{"LocusID#3"}), and
#'       \code{retry_failed_ids} (locus names that failed even after sequential
#'       retry).}
#'   }
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet subseq replaceLetterAt writeXStringSet
#' @importFrom pwalign pairwiseAlignment nmismatch
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply mclapply
#'
#' @export
rhampseq2madc <- function(hap_genotype_file,
                          haplotype_allele_fasta,
                          n_cores = 1,
                          prefix = NULL,
                          verbose = TRUE,
                          parallel_type = NULL) {
  vmsg("Running BIGr rhampseq2madc", verbose = verbose, level = 0, type = ">>")
  vmsg("hap_genotype_file      : %s", verbose = verbose, level = 1, type = ">>", hap_genotype_file)
  vmsg("haplotype_allele_fasta : %s", verbose = verbose, level = 1, type = ">>", haplotype_allele_fasta)
  vmsg("n_cores                : %s", verbose = verbose, level = 1, type = ">>", n_cores)
  vmsg("prefix                 : %s", verbose = verbose, level = 1, type = ">>", if (is.null(prefix)) "NULL" else prefix)

  vmsg("Reading input files", verbose = verbose, level = 0, type = ">>")
  hapgeno <- read.table(hap_genotype_file, sep = "\t", header = TRUE)
  sequences <- readDNAStringSet(haplotype_allele_fasta)
  vmsg("%s loci and %s sequences read",
       verbose = verbose, level = 1, type = ">>",
       nrow(hapgeno), length(sequences)
  )

  vmsg("Checking inputs", verbose = verbose, level = 0, type = ">>")
  # --- Input validation ---
  if (ncol(hapgeno) < 3) {
    stop("'hap_genotype_file' must have at least 3 columns: Locus, Haplotypes, and at least one sample column.")
  }

  if (!is.character(hapgeno[[1]])) {
    stop("Column 1 of 'hap_genotype_file' (Locus) must contain character locus names.")
  }

  # Sample columns: each cell must be a genotype like "12/6:101,72", "1:58", "./.:0", or NA
  geno_pattern <- "^(\\./\\.|[0-9]+(/[0-9]+)*):[0-9]+(,[0-9]+)*$"
  sample_data <- hapgeno[, 3:ncol(hapgeno), drop = FALSE]
  bad_cells <- !apply(sample_data, 1:2, function(x) is.na(x) | grepl(geno_pattern, x))
  if (any(bad_cells)) {
    bad_coords <- which(bad_cells, arr.ind = TRUE)
    stop(paste0(
      "Unexpected genotype format in sample column(s) of 'hap_genotype_file'. ",
      "First offending cell — row: ", bad_coords[1, 1],
      ", column: ", colnames(sample_data)[bad_coords[1, 2]],
      ", value: '", sample_data[bad_coords[1, 1], bad_coords[1, 2]], "'."
    ))
  }

  # FASTA sequence names must follow "LocusID#N" convention
  invalid_seq_names <- !grepl("#[0-9]+$", names(sequences))
  if (any(invalid_seq_names)) {
    stop(paste0(
      "'haplotype_allele_fasta' sequence names must end with '#N' (e.g. 'LocusID#1'). ",
      "First offending name: '", names(sequences)[which(invalid_seq_names)[1]], "'."
    ))
  }

  # Cross-platform parallel backend
  # Pre-build locus-prefix → sequence-indices map once, avoiding O(n_loci × n_seqs) grep per worker
  seq_prefix_index <- split(seq_along(sequences), sub("#[0-9]+$", "", names(sequences)))

  vmsg("Processing %s loci with %s core(s)",
       verbose = verbose, level = 0, type = ">>",
       nrow(hapgeno), n_cores
  )
  # mclapply (fork) is unavailable on Windows; use a PSOCK cluster there instead
  if (n_cores > 1) {
    if (.Platform$OS.type == "windows") {
      cl <- makeCluster(n_cores, type = "PSOCK")
    } else if (!is.null(parallel_type) && parallel_type == "PSOCK") {
      cl <- makeCluster(n_cores, type = "PSOCK")
    } else {
      cl <- makeCluster(n_cores, type = "FORK")
    }
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl,
                  varlist = c("hapgeno", "sequences", "verbose", "seq_prefix_index"),
                  envir = environment()
    )
    clusterEvalQ(cl, {
      library(Biostrings)
      library(pwalign)
    })
    par_fun <- function(X, FUN) {
      safe <- function(x) tryCatch(FUN(x), error = function(e) {
        structure(list(message = conditionMessage(e)), class = "try-error")
      })
      parLapply(cl, X, safe)
    }
  } else {
    par_fun <- function(X, FUN) lapply(X, FUN)
  }

  .locus_worker <- function(t) {
    ## Debug code
    locus <- unique(hapgeno$Locus)
    t <- which(hapgeno$Locus == locus[20])
    ###

    .info <- list(fallback_ref = FALSE, dup_ref = FALSE)
    onetag <- hapgeno[t, ]

    # Extract read depths from sample columns (col 3 to last)
    sample_cols <- as.data.frame(onetag[, 3:ncol(onetag), drop = FALSE])
    sample_names <- colnames(sample_cols)

    # Vectorised: bulk-split all sample cells, build one data.frame at the end
    vals <- unname(vapply(sample_cols[1, ], as.character, character(1)))
    has_colon <- grepl(":", vals, fixed = TRUE) & !is.na(vals)
    vals_v <- vals[has_colon]
    samps_v <- sample_names[has_colon]

    if (length(vals_v) == 0L) {
      long_df <- NULL
    } else {
      parts <- strsplit(vals_v, ":", fixed = TRUE)
      allele_parts <- vapply(parts, "[[", character(1), 1L)
      depth_parts <- vapply(parts, "[[", character(1), 2L)

      ok <- !grepl(".", allele_parts, fixed = TRUE)
      allele_parts <- allele_parts[ok]
      depth_parts <- depth_parts[ok]
      samps_v <- samps_v[ok]

      if (length(allele_parts) == 0L) {
        long_df <- NULL
      } else {
        alleles_list <- strsplit(allele_parts, "/", fixed = TRUE)
        depths_list <- strsplit(depth_parts, ",", fixed = TRUE)
        # Collapse homozygous calls with a single depth (e.g. "2/2:39" → allele 2, depth 39).
        keep <- mapply(function(a, d) {
          if (length(d) < length(a) && length(unique(a)) == 1L) seq_len(1L) else seq_along(a)
        }, alleles_list, depths_list, SIMPLIFY = FALSE)
        alleles_list <- mapply(function(a, k) a[k], alleles_list, keep, SIMPLIFY = FALSE)
        depths_list  <- mapply(function(d, k) d[k], depths_list,  keep, SIMPLIFY = FALSE)
        lens <- lengths(alleles_list)
        long_df <- data.frame(
          allele_id = as.integer(unlist(alleles_list, use.names = FALSE)),
          depth = as.integer(unlist(depths_list, use.names = FALSE)),
          sample = rep(samps_v, lens),
          stringsAsFactors = FALSE
        )
      }
    }

    # Build depth matrix directly — avoids slow reshape()
    if (is.null(long_df) || nrow(long_df) == 0L) {
      depth_wide <- data.frame(allele_id = integer(0), check.names = FALSE)
      for (s in sample_names) depth_wide[[s]] <- integer(0)
    } else {
      all_allele_ids <- sort(unique(long_df$allele_id))
      depth_mat <- matrix(0L,
                          nrow = length(all_allele_ids), ncol = length(sample_names),
                          dimnames = list(NULL, sample_names)
      )
      depth_mat[cbind(
        match(long_df$allele_id, all_allele_ids),
        match(long_df$sample, sample_names)
      )] <- long_df$depth
      depth_wide <- data.frame(allele_id = all_allele_ids, depth_mat, check.names = FALSE)

      # Samples absent from depth_wide (all calls missing) get an all-zero column
      missing_samples <- setdiff(sample_names, colnames(depth_wide))
      for (s in missing_samples) depth_wide[[s]] <- 0L
    }

    # Define IDs
    cloneID <- as.character(onetag[[1]])
    AlleleIDs_idx <- seq_prefix_index[[cloneID]]
    if (is.null(AlleleIDs_idx)) AlleleIDs_idx <- integer(0)

    if (length(AlleleIDs_idx) == 0) {
      # Return a skip marker — vmsg inside a worker goes to worker stdout, not the console.
      # The main process will collect and report all skip reasons.
      return(list(.skipped = TRUE, reason = "no sequences found in FASTA"))
    }

    reftag_n <- paste0(cloneID, "#1")
    reftag_idx <- which(names(sequences) == reftag_n)

    if (length(reftag_idx) == 0) {
      # #1 absent — fall back to the allele with the smallest available ID
      all_allele_nums <- as.integer(sapply(strsplit(names(sequences)[AlleleIDs_idx], "#"), "[[", 2))
      fallback_num <- min(all_allele_nums)
      reftag_n <- paste0(cloneID, "#", fallback_num)
      reftag_idx <- which(names(sequences) == reftag_n)
      vmsg(
        "Locus '%s': reference sequence '#1' not found in FASTA; using '#%s' as reference instead.",
        verbose = verbose, level = 1, type = ">>",
        cloneID, fallback_num
      )
      .info$fallback_ref <- TRUE
    }
    reftag <- sequences[[reftag_idx]]

    # Align and count mismatches + total indel bases
    results <- lapply(AlleleIDs_idx, function(i) {
      aln <- pairwiseAlignment(sequences[[i]], reftag, type = "global")
      indels <- rbind(as.data.frame(aln@subject@indel@unlistData), as.data.frame(aln@pattern@indel@unlistData))
      list(
        data.frame(
          mismatches = nmismatch(aln),
          indels = nrow(indels)
        ),
        # Use subject (reftag) coordinates so positions are always within reference bounds
        pos_mismatches = aln@subject@mismatch@unlistData,
        pos_indels = indels
      )
    })

    # Transpose for easier reading
    results_df <- t(sapply(results, "[[", 1))

    # positions
    pos_mismatches <- sapply(results, "[[", 2)

    # Indels return a table with start, end and width
    pos_indels <- lapply(results, "[[", 3)

    rownames(results_df) <- names(pos_indels) <- names(pos_mismatches) <- names(sequences)[AlleleIDs_idx]

    # Check for sequences identical to the reference (#1);
    # when duplicates exist, merge their depths into the reference row.
    refs_tag_idx <- which(results_df[, 1] == 0 & results_df[, 2] == 0)

    if (length(refs_tag_idx) > 1) {
      dup_ref_seqnames <- rownames(results_df)[refs_tag_idx]
      extra_ref_names <- dup_ref_seqnames[dup_ref_seqnames != reftag_n]
      extra_ids <- as.integer(sapply(strsplit(extra_ref_names, "#"), "[[", 2))

      vmsg(
        "Locus '%s': %d sequence(s) identical to reference (#1) detected: %s — their read depths will be summed into the reference.",
        verbose = verbose, level = 1, type = ">>",
        cloneID, length(extra_ref_names), paste(extra_ref_names, collapse = ", ")
      )
      .info$dup_ref_seqnames <- extra_ref_names

      # Merge extra-ref depths into the reference row of depth_wide
      ref_row_idx <- match(1L, depth_wide$allele_id)
      extra_row_idx <- match(extra_ids, depth_wide$allele_id)
      extra_row_idx <- extra_row_idx[!is.na(extra_row_idx)]

      if (!is.na(ref_row_idx) && length(extra_row_idx) > 0) {
        depth_wide[ref_row_idx, -1] <- depth_wide[ref_row_idx, -1] +
          colSums(depth_wide[extra_row_idx, -1, drop = FALSE])
      } else if (is.na(ref_row_idx) && length(extra_row_idx) > 0) {
        # Reference allele absent from genotype data; re-label first duplicate as allele 1
        depth_wide[extra_row_idx[1], "allele_id"] <- 1L
        extra_row_idx <- extra_row_idx[-1]
      }
      if (length(extra_row_idx) > 0) {
        depth_wide <- depth_wide[-extra_row_idx, , drop = FALSE]
      }

      # Exclude merged alleles from non-ref processing
      pos_mismatches <- pos_mismatches[!names(pos_mismatches) %in% extra_ref_names]
      pos_indels <- pos_indels[!names(pos_indels) %in% extra_ref_names]
    }

    # Make a fake one target SNP in a position that polymorphis still doesn't exist
    # this way all tags will be RefMatch simplifying the code
    # the position of the fake one will depend on the distribution of existing polymorphisms
    # the fake mutation should be created on regions with less existing true polymorphisms

    # Defining target position
    # Collect all polymorphic positions (mismatches + indels) across all alleles
    non_ref_names <- names(pos_mismatches)[names(pos_mismatches) != reftag_n]

    all_poly_pos <- unique(c(
      unlist(pos_mismatches[non_ref_names]),
      unlist(lapply(pos_indels[non_ref_names], function(ir) {
        if (length(ir) == 0) {
          return(integer(0))
        }
        unlist(Map(seq.int, ir$start, ir$start + ir$width - 1L))
      }))
    ))

    # Clip to valid reference positions (guard against alignment-coordinate overflow)
    seq_len <- length(reftag)
    all_poly_pos <- all_poly_pos[all_poly_pos >= 1L & all_poly_pos <= seq_len]

    poly_sorted <- sort(all_poly_pos)
    boundaries <- c(0L, poly_sorted, seq_len + 1L)
    gap_sizes <- diff(boundaries) - 1L
    max_gap_idx <- which.max(gap_sizes)
    gap_start <- boundaries[max_gap_idx] + 1L
    gap_end <- boundaries[max_gap_idx + 1L] - 1L

    if (gap_start <= gap_end) {
      # Normal case: at least one clean (polymorphism-free) position exists
      mid_pos <- as.integer(round((gap_start + gap_end) / 2))
    } else {
      # Every position is covered by at least one polymorphism.
      # Build a per-position count of how many alleles are polymorphic there,
      # then pick the centre of the longest run at the minimum count.
      vmsg(
        "Locus '%s': all reference positions are polymorphic; selecting the least-affected position for the synthetic SNP.",
        verbose = verbose, level = 1, type = ">>",
        cloneID
      )
      poly_cov <- integer(seq_len)
      for (nm in non_ref_names) {
        mm <- pos_mismatches[[nm]]
        if (length(mm) > 0) {
          mm <- mm[mm >= 1L & mm <= seq_len]
          poly_cov[mm] <- poly_cov[mm] + 1L
        }
        ir <- pos_indels[[nm]]
        if (length(ir) > 0) {
          ip <- unlist(Map(seq.int, ir$start, ir$start + ir$width - 1L))
          ip <- ip[ip >= 1L & ip <= seq_len]
          poly_cov[ip] <- poly_cov[ip] + 1L
        }
      }
      min_pos <- which(poly_cov == min(poly_cov))
      run_breaks <- which(diff(min_pos) > 1L)
      run_starts <- c(1L, run_breaks + 1L)
      run_ends <- c(run_breaks, length(min_pos))
      best_run_idx <- which.max(run_ends - run_starts)
      best_run <- min_pos[run_starts[best_run_idx]:run_ends[best_run_idx]]
      mid_pos <- as.integer(best_run[ceiling(length(best_run) / 2)])
    }
    mid_base <- as.character(subseq(reftag, mid_pos, mid_pos))

    ## Has a specific change dynamic
    new_base <- switch(mid_base,
                       "A" = "C",
                       "C" = "A",
                       "G" = "T",
                       "T" = "G"
    )

    alttag <- replaceLetterAt(reftag, mid_pos, new_base)

    idx_non_ref_names <- match(non_ref_names, names(sequences))
    ids_real <- sapply(strsplit(non_ref_names, "#"), "[[", 2)
    # Alt sequence is fake and must have id 0002, so all the others gain a +1
    ids_new <- as.numeric(ids_real) + 1
    # padding with zeros to match the four digits IDs
    ids_new <- sprintf("%04d", as.integer(ids_new))

    if(length(ids_new) >0){
      AlleleID <- c(
        paste0(cloneID, "|Ref_0001"),
        paste0(cloneID, "|Alt_0002"),
        paste0(cloneID, "|RefMatch_", ids_new)
      )
    } else {
      AlleleID <- c(
        paste0(cloneID, "|Ref_0001"),
        paste0(cloneID, "|Alt_0002")
      )
    }

    AlleleSequence <- c(
      as.character(reftag),
      as.character(alttag),
      sapply(sequences[idx_non_ref_names], as.character)
    )
    names(AlleleSequence) <- NULL

    madc13_one <- data.frame(
      AlleleID = AlleleID,
      CloneID = cloneID,
      AlleleSequence = unlist(AlleleSequence)
    )

    tomerge_idx <- match(c(1, "alt", ids_real), depth_wide$allele_id)

    madc_one <- cbind(madc13_one, depth_wide[tomerge_idx, -1])

    # Replace NAs in all sample columns (unobserved alleles and missing ./.:0 calls)
    madc_one[, -(1:3)][is.na(madc_one[, -(1:3)])] <- 0L

    list(
      madc = madc_one,
      target_pos = data.frame(
        AlleleID_rhAmpSeq = reftag_n,
        AlleleID          = cloneID,
        target_position   = mid_pos,
        target_base       = paste0(mid_base, "/", new_base),
        stringsAsFactors  = FALSE
      ),
      info = .info
    )
  }
  madc_list <- par_fun(seq_len(nrow(hapgeno)), .locus_worker)

  # Loci that errored inside the parallel worker are returned as try-error objects.
  # Retry them sequentially — they usually succeed outside the parallel context.
  retry_failed_ids <- character(0)
  failed_idx <- which(vapply(madc_list, inherits, logical(1), what = "try-error"))
  if (length(failed_idx) > 0) {
    vmsg("%s locus/loci failed in parallel; retrying sequentially.",
         verbose = verbose, level = 1, type = ">>", length(failed_idx))
    for (t in failed_idx) {
      madc_list[[t]] <- tryCatch(
        .locus_worker(t),
        error = function(e) {
          retry_failed_ids <<- c(retry_failed_ids, as.character(hapgeno[[1]][t]))
          vmsg("Locus '%s': sequential retry also failed \u2014 %s",
               verbose = verbose, level = 1, type = ">>",
               hapgeno[[1]][t], conditionMessage(e))
          NULL
        }
      )
    }
  }

  # Report loci intentionally skipped inside workers and convert to NULL
  # (vmsg inside workers writes to worker stdout, not the console)
  no_fasta_ids <- character(0)
  skipped_mask <- vapply(madc_list, function(x) is.list(x) && isTRUE(x$.skipped), logical(1))
  if (any(skipped_mask)) {
    skip_idx <- which(skipped_mask)
    no_fasta_ids <- as.character(hapgeno[[1]][skip_idx])
    for (i in skip_idx) {
      vmsg("Locus '%s': %s — skipped.",
           verbose = verbose, level = 1, type = ">>",
           hapgeno[[1]][i], madc_list[[i]]$reason)
    }
    madc_list[skip_idx] <- list(NULL)
  }

  fallback_ref_ids <- as.character(hapgeno[[1]][vapply(madc_list, function(x) is.list(x) && isTRUE(x$info$fallback_ref), logical(1))])
  dup_ref_ids      <- unlist(lapply(madc_list, function(x) {
    if (is.list(x) && length(x$info$dup_ref_seqnames) > 0) x$info$dup_ref_seqnames else character(0)
  }), use.names = FALSE)

  # Drop all NULL entries (skipped + sequential-retry failures)
  madc_list <- Filter(Negate(is.null), madc_list)

  vmsg("Assembling results from %s processed loci",
       verbose = verbose, level = 0, type = ">>",
       length(madc_list)
  )
  madc_final <- do.call(rbind, lapply(madc_list, "[[", "madc"))
  target_positions <- do.call(rbind, lapply(madc_list, "[[", "target_pos"))
  rownames(target_positions) <- NULL

  new_fasta <- DNAStringSet(madc_final$AlleleSequence)
  names(new_fasta) <- madc_final$AlleleID

  if (!is.null(prefix)) {
    vmsg("Writing output files with prefix '%s'", verbose = verbose, level = 0, type = ">>", prefix)
    write.csv(madc_final, paste0(prefix, "_madc.csv"), row.names = FALSE)
    write.csv(target_positions, paste0(prefix, "_target_positions.csv"), row.names = FALSE)
    writeXStringSet(new_fasta, paste0(prefix, ".fasta"))
  }

  vmsg("Done!", verbose = verbose, level = 0, type = ">>")
  list(
    madc             = madc_final,
    target_positions = target_positions,
    new_fasta        = new_fasta,
    info             = list(
      no_fasta_ids     = no_fasta_ids,
      fallback_ref_ids = fallback_ref_ids,
      dup_ref_ids      = dup_ref_ids,
      retry_failed_ids = retry_failed_ids
    )
  )
}
