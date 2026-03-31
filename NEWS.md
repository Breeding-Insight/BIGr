# BIGr 0.6.5

# Updates on madc2vcf functions
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

