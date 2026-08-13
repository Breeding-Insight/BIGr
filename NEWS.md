# BIGr 0.10.0

* New `madc_summary()`: read-only summary tables for a fixed allele ID MADC - a per-marker microhaplotype frequency distribution, per-sample and per-marker total read depth (`depth`) and on-target `|Ref`/`|Alt` depth (`target_depth`) alongside missing rate, number of mHaps, and off-target read fraction, a missingness threshold sweep, an overall summary, and optional aggregation by a sample-metadata category (e.g. species).
* New `madc_plot()`: visualize a fixed allele ID MADC via a `plot.type` selector - a PCA of alternative read ratios (color by a metadata category, optional second variable by point `shape.col`, custom `palette`), a linear marker-distribution genome plot, a markers x samples read-depth/#mHaps heatmap (failed cells black; depth on a colorblind-safe diverging scale; samples sorted/labeled by category and optionally faceted by chromosome), a per-sample missing-data boxplot (optionally ordered by median missing rate via `miss.sort` and drawn horizontally via `horizontal`), and a publication-oriented `circos` circular plot. The `circos` figure carries four concentric tracks - dominant per-marker locus ticks with a distinct Mb coordinate axis ("Loci"), a marker-density heatmap ring ("Density", window set by `density.window`), log-scaled #mHaps two-tone bars - a neutral base with a red over-threshold segment above a dashed reference ring - flagging paralog-suspect loci ("mHaps", threshold `paralog.flag`), and a gapless depth-QC ribbon flagging low-depth and high-depth-outlier loci ("Depth QC") - identified by horizontal track titles in a 12 o'clock corridor (`label.gap`) with a single compact legend in the center; track colors are customizable (`mhap.col`, `depth.col`, `density.col`). Requesting several (non-circos) types returns a multi-panel figure; figures can be saved via `output.file`. Plotting uses `ggplot2` + base `grid`; the `circos` plot uses the suggested `circlize` package.
* `filterMADC` updates:
    * Now only accepts fixed allele ID MADC files (processed through HapApp). The input is validated with `check_madc_sanity` and raw DArT MADC files (or any file without fixed AlleleIDs) are rejected with an error. This also fixes a crash on raw MADC files caused by non-unique AlleleIDs.
    * Removed the `n.summary.columns` argument (breaking change). Fixed allele ID MADC has a fixed column layout, so summary columns no longer need to be specified or detected.
    * Replaced the per-mhap `min.mean.reads`/`max.mean.reads` arguments with a per-locus depth window `min.locus.depth`/`max.locus.depth` (breaking change). Depth is the mean reads per sample at a locus (sum of all mhap reads at the locus divided by the number of samples), and loci outside the window are removed whole (Ref and Alt together) - e.g. `max.locus.depth` to exclude paralogous over-amplification.
    * `max.mhaps.per.loci` and `target.only` now retain strictly the target `|Ref`/`|Alt` alleles at affected loci, also removing `|Other` alleles (previously only `|RefMatch`/`|AltMatch` were removed).
    * The `min.ind.with.reads` presence filter now protects the target `|Ref`/`|Alt` alleles, only pruning non-target (`|RefMatch`/`|AltMatch`/`|Other`) mhaps, so a locus never loses its Ref/Alt pair.


# BIGr 0.9.0

* `check_madc_sanity` updates: distinguish presence of IUPAC codes on REF/ALT (return logical variable IUPACcodes) from RefMatch/AltMatch/Others (returned logical variable IUPACcodes_MatchAlleles) and from Identical IUPAC code in identical positions in REF/ALT (returned logical variable IUPACcodes_IdenticalRefAlt)
* Two new arguments to the `madc2vcf_all` function:
    * `others_min_dist` (default 5bp)
    * `others_max_close_snps` (default 3 SNPs)

* By default, `Other` tags will be discarded if they have more than 3 SNPs with less than 5bp distance between them   
* Adapt `madc2vcf_all` code to let pass identical IUPAC code in identical position in REF/ALT but ignore polymorphims in Match or Other alleles at the same position
* `madc2vcf_targets` and `madc2vcf_multi` will only throw an error if different IUPAC are present in REF/ALT sequences (IUPACcodes = TRUE)


# BIGr 0.8.1

- Fixed `filterVCF()` discarding SNPs whose INFO values are written in scientific notation. The `OD`, `BIAS`, and `PMC` filters read INFO values with a pattern that accepted only digits and decimal points, so a value such as `PMC=4.78e-07` was read as `4.78` and failed a `filter.PMC = 0.05` threshold that it actually passed by five orders of magnitude. `updog2vcf()` writes these values with `paste0()`, so R's default formatting produces scientific notation for small numbers, and updog's `od` and `prop_mis` are routinely small for good markers - meaning the SNPs most likely to be dropped were the cleanest ones. INFO values are now parsed with `as.numeric()` on the whole field, accepting any valid numeric representation. Field names are also matched against a complete INFO entry rather than as a substring, so a field such as `XOD` is no longer mistaken for `OD`.
- Fixed `filterVCF()` writing corrupt all-`NA` variants when a filter value could not be read. A logical index containing `NA` does not drop a row in R, it inserts a row of `NA`s, so any variant with an unreadable `OD`, `BIAS` or `PMC` value, or an `NA` minor allele frequency, was written into the output with an empty `CHROM`, `POS`, `REF` and `ALT`. The `NA` frequency case was reachable through normal use, because `vcfR::maf()` returns `NA` for a variant left with no called genotypes and the `filter.DP` and `filter.MPP` masking can produce exactly that. These variants are now removed, and the number removed is reported in a warning so that they are never discarded silently.
- `filterVCF()` now removes variants when none of them have a readable value for a requested filter, instead of returning the data unfiltered. Previously the `OD`, `BIAS` and `PMC` filters skipped filtering entirely in that situation and reported only "No valid values found", so a run that filtered nothing could look like a run that succeeded. All of the filters now behave the same way, and the warning names the likely cause so it can be diagnosed.
- Fixed `filterVCF()` skipping a requested `OD`, `BIAS` or `PMC` filter when the first variant in the file did not carry that INFO field. The available INFO fields were read from the first record only, so a field that was absent from that one record disabled its filter for the entire file, with no warning and no indication in the output that filtering had not happened. A field may legally be absent from any individual record, so requested filters are now applied to every record, and records without a usable value are removed and reported. This also removes the last case where the order of records in the file could change the result of filtering.
- Corrected and expanded the documentation for `filterVCF()`. `filter.DP` was described as a total read depth filter applied to each SNP; it is applied to each genotype call using the `DP` value in the FORMAT column, setting calls below the threshold to missing rather than removing variants. `filter.MPP` works the same way and was previously undescribed. The remaining filters now state which field they read and which direction the comparison runs, `filter.SAMPLE.miss` and `filter.SNP.miss` state that they take a proportion rather than a percentage, and `filter.BIAS.min` and `filter.BIAS.max` state that neither has any effect unless both are supplied. The return value is documented as a vcfR object when `output.file` is omitted, rather than always a gzipped file.
- Updated madc2vcf_all and madc2vcf_targets. Before, it was possible for POS to be exported as scientific notation instead of integers, and for negative POS values to be present for off target SNPs. POS are corrected to be integers, and SNPs with a negative POS value are removed.
- Remove BIGpopA functions - now it is a independent package: https://github.com/Breeding-Insight/BIGpopA 
- Fixed `madc2vcf_all()` error "invalid substring arguments" that occurred with `add_others = TRUE` when an off-target ("Other") allele aligned to the reference with no mismatch positions remaining after the target SNP position was removed. The reference/alternate base lookups for off-target alleles are now guarded by the existing non-empty check, matching how the off-target Match alleles are already handled.

# BIGr 0.7.2

- Fixed manual text errors

# BIGr 0.7.1

- Updated `check_ped()` to return corrected pedigree data in the result list instead of assigning objects to the global environment
- Skipped long remote `madc2vcf_all` integration tests on CRAN while keeping them enabled in GitHub Actions

# BIGr 0.7.0

## Updates on `dosage2vcf`

- Added support for DArT SNP/INDEL 1-row and 2-row report formats
- `dosage2vcf` now validates marker and sample sets between report and counts files, then aligns counts to the report order before writing VCF genotypes
- VCF `CHROM` and `POS` are derived from `Chrom`/`ChromPos` when present, otherwise from `MarkerName`; `MarkerName` is retained in the VCF `ID` field
- Missing SNP/INDEL genotype calls (`-`/`NA`) are written as diploid missing genotypes (`./.`)

## New function `madc2vcf_multi`

- New function `madc2vcf_multi` to convert a DArTag MADC file to a VCF using the polyRAD pipeline for multiallelic genotyping
- Runs `check_madc_sanity` before loading the data and stops with informative errors if:
    - Required columns are missing
    - IUPAC (non-ATCG) codes are present in AlleleSequence
    - Ref/Alt sequences are unpaired (`RefAltSeqs = FALSE`)
    - Allele IDs have not been fixed by HapApp (`FixAlleleIDs = FALSE`)
    - CloneIDs do not follow `Chr_Pos` format and no `markers_info` is provided
- New argument `markers_info`: optional path or URL to a CSV with `CloneID`/`BI_markerID`, `Chr`, and `Pos` columns; required when CloneIDs do not follow the `Chr_Pos` format
- Runs `check_botloci` to validate and reconcile CloneIDs between the MADC and botloci file, automatically fixing padding mismatches
- A corrected temp file is written and passed to `readDArTag` only when needed (all-NA rows/columns detected, CloneIDs remapped by `check_botloci`, or botloci IDs remapped)
- Accepts paths or URLs for `madc_file`, `botloci_file`, and `markers_info`
- Estimates overdispersion with `polyRAD::TestOverdispersion`, iterates priors with `polyRAD::IterateHWE`, and exports the result with `polyRAD::RADdata2VCF`
- `polyRAD` is a soft dependency (listed under `Suggests`); an informative error is raised if it is not installed

# BIGr 0.6.6

## Updates on `madc2vcf_all`

- New arguments for controlling processing of `Other` alleles:
    - `add_others`: if `TRUE` (default), alleles labeled "Other" in the MADC are included in off-target SNP extraction
    - `others_max_snps`: discards Other alleles with more than this many SNP differences relative to the Ref sequence (default: 5)
    - `others_rm_with_indels`: discards Other alleles containing insertions or deletions relative to the Ref sequence (default: `TRUE`)
- Others alleles that carry a different base at the target SNP position are now reported as a 3rd allele in the VCF instead of being silently dropped
- Target position is now correctly removed from Others alignments, preventing duplicate VCF positions and marker IDs
- Fixed a bug where Others alleles with "Ref_" or "Alt_" in their AlleleID would corrupt the target SNP REF/ALT fields and read depth counts in `merge_counts`
- Improved verbose messages throughout: counts of Other alleles found, kept, and discarded (by indel filter and by max SNP filter) are now reported; multiallelic target SNPs with a 3rd allele from Others are counted and reported
- Debug-level message (level 3) listing each Other allele added and its genomic position

# BIGr 0.6.5

## Updates on madc2vcf functions
Details:

- both functions targets and all (targets + off-targets) markers now have `check_madc_sanity` function implemented. It tests:
    - [Columns] If MADC has the expected columns
    - [allNArow | allNAcol] Presence of columns and rows with all NA (happens often when people open the MADC in excel before loading in R)
    - [IUPACcodes] Presence of IUPAC codes on AlleleSequence
    - [LowerCase] Presence of lower case bases on AlleleSequence
    - [Indels] Presence of Indels
    - [ChromPos] If CloneID follows the format Chr_Pos
    - [RefAltSeqs] If all Ref Allele has corresponding Alt and vice-versa
    - [OtherAlleles] If "Other" exists in the MADC AlleleID

- Better messages if `verbose = TRUE` in `madc2vcf_all`
- `madc2vcf_all` support for Indels - markers_info with Indels position is required; only the target indel is extracted, off-targets are ignored for the tag
- `madc2vcf_targets` doesn’t run if: 
    - MADC Column names are not correct
    - Ignore Other alleles - but inform the user if they exist or not and direct them to `madc2vcf_all` in case they want to extract them as well
- See the table for madc2vcf_targets requirements accordingly to MADC content:

  | check status | get_REF_ALT | Requires
-- | -- | -- | --
IUPAC | TRUE | TRUE | markers_info REF/ALT
  | TRUE | FALSE | -
  | FALSE | TRUE | botloci or markers_info REF/ALT
  | FALSE | FALSE | -
Indels | TRUE | TRUE | markers_info REF/ALT
  | TRUE | FALSE | -
  | FALSE | TRUE | botloci or markers_info REF/ALT
  | FALSE | FALSE | -
ChromPos | TRUE | TRUE | botloci or markers_info REF/ALT
  | TRUE | FALSE | -
  | FALSE | TRUE | markers_info CHR/POS/REF/ALT or markers_info CHR/POS/ + botloci
  | FALSE | FALSE | markers_info CHR/POS
FixAlleleIDs | TRUE | TRUE | botloci or markers_info REF/ALT
  | TRUE | FALSE | -
  | FALSE | TRUE | markers_info REF/ALT
  | FALSE | FALSE | -

# BIGr 0.6.4

- Add function `vmsg` to organize messages printed on the console
- Add metadata to VCF header from madc2vcf_targets
- Add argument `madc_object` to `get_countsMADC` to avoid reading the MADC file twice and to get directly the MADC fixed padding output from `check_botloci`
- Organize messages from `madc2vcf_targets` checks
- Add argument `collapse_matches_counts` and `verbose` to `madc2vcf_targets` function

# BIGr 0.6.3

- New function to check MADC files: `check_madc_sanity`. Currently, it checks for the presence of required columns, whether fixed allele IDs were assigned, the presence of IUPAC codes, lowercase sequence bases, indels, and chromosome and position information.
- Added new argument `markers_info`, which allows users to provide a CSV file with marker information such as CHROM, POS, marker type, and position of indels. For BI species, this information is available from [PanelHub](https://github.com/Breeding-Insight/BIGapp-PanelHub).
- Checked inputs for `madc2vcf_all`.
- Updated affiliation in `DESCRIPTION`.

# BIGr 0.6.2

- Fixed the doi and name list in the CITATION file

# BIGr 0.6.1

- Added new functions for filtering MADC files and converting to relationship matrices
- Added function thinSNPs to thin SNPs based on physical distance
- Added bug fixes and improvements to existing functions

# BIGr 0.5.5

- Updated DESCRIPTION
- Added return value for merge_MADCs
- Added optional seed for check_ped
- Added verbose option

# BIGr 0.5.4

-   Updated dosage2vcf example

# BIGr 0.5.3

-   Updated madc2vcf_all example

# BIGr 0.5.2

-   madc2vcf function changed to madc2vcf_targets
-   get_OffTargets function changed to madc2vcf_all
-   Updates to testthat tests and function examples

# BIGr 0.5.1

-   Improvements of testthat tests
-   Add check_replicates and check_homozygous_trios for pedigree relationship quality check

# BIGr 0.5.0

-   Add imputation_concordance function to estimate accuracy of imputed and original dataset
-   Add get_OffTargets function to extract target and off-target SNPs from a MADC file
-   Add merge_MADCs function to merge two or more MADC files together
-   Improved documentation and examples for all functions
-   Add tests for all functions

# BIGr 0.3.3

-   Adapt updog2vcf to model f1, f1pp, s1 and s1pp

# BIGr 0.3.2

-   updog2vcf function option to output compressed VCF (.vcf.gz) - set as default
-   remove need for defining ploidy
-   add metadata at the VCF header
