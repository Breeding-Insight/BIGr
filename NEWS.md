# BIGr 0.9.0

* `check_madc_sanity` updates: distinguish presence of IUPAC codes on REF/ALT (return logical variable IUPACcodes) from RefMatch/AltMatch/Others (returned logical variable IUPACcodes_MatchAlleles) and from Identical IUPAC code in identical positions in REF/ALT (returned logical variable IUPACcodes_IdenticalRefAlt)
* Two new arguments to the `madc2vcf_all` function:
    * `others_min_dist` (default 5bp)
    * `others_max_close_snps` (default 3 SNPs)

By default, “Other“ tags will be discarded if they have more than 3 SNPs with less than 5bp distance between them   
* Adapt `madc2vcf_all` code to let pass identical IUPAC code in identical position in REF/ALT but ignore polymorphims in Match or Other alleles at the same position
* `madc2vcf_targets` and `madc2vcf_multi` will only throw an error if different IUPAC are present in REF/ALT sequences (IUPACcodes = TRUE)


# BIGr 0.8.1

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
