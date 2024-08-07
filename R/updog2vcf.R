#' Export Updog Results as VCF
#'
#' This function will convert an Updog output to a VCF file
#'
#' When performing dosage calling for multiple SNPs using Updog, the output file contains information for all loci and all samples.
#' This function will convert the updog output file to a VCF file, while retaining the information for the values that are commonly
#' used to filter low quality and low confident dosage calls.
#'
#' @param multidog.object updog output object with class "multidog" from dosage calling
#' @param output.file output file name and path
#' @param updog_version character defining updog package version used to generate the multidog object
#' @param compress logical. If TRUE returns a vcf.gz file
#' @return A vcf file
#' @import dplyr
#' @import tidyr
#' @importFrom Rdpack reprompt
#' @references
#' Updog R package
#' @export
updog2vcf <- function(multidog.object, output.file, updog_version = NULL, compress = TRUE) {
  #Making the VCF (This is highly dependent on snps being in a format where the SNP IDs are the CHR_POS)

  mout <- multidog.object
  ploidy <- as.numeric(unique(multidog.object$snpdf$ploidy))
  if(!grepl(".vcf", output.file)) output.file <- paste0(output.file,".vcf") #Make sure it ends with the vcf file extension
  model_selected <- unique(multidog.object$snpdf$model)

  updog_meta <- paste0('##UpdogCommandLine.multidog=<ID=Multidog,Version="',
                       updog_version, '",CommandLine="> multidog(refmat = matrices$ref_matrix, sizemat = matrices$size_matrix, ploidy = ',ploidy,
                       ', model = "',model_selected,'")>')

  bigr_meta <- paste0('##BIGrCommandLine.updog2vcf=<ID=updog2vcf,Version"',
                      packageVersion("BIGr"), '",Data="',
                      Sys.time(),'", CommandLine="> updog2vcf(',deparse(substitute(multidog.object)),',',
                      deparse(substitute(output.file)), ',',
                      deparse(substitute(updog_version)),')>')

  #Make a header separate from the dataframe
  vcf_header <- c(
    "##fileformat=VCFv4.3",
    "##reference=NA",
    "##contig=<ID=NA,length=NA>",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##INFO=<ID=ADS,Number=R,Type=Integer,Description="Depths for the ref and each alt allele in the order listed">',
    '##INFO=<ID=BIAS,Number=1,Type=Float,Description="The estimated allele bias of the SNP from updog">',
    '##INFO=<ID=OD,Number=1,Type=Float,Description="The estimated overdispersion parameter of the SNP from updog">',
    '##INFO=<ID=PMC,Number=1,Type=Float,Description="The estimated proportion of individuals misclassified in the SNP from updog">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, where 1 is the count of alternate alleles">',
    '##FORMAT=<ID=UD,Number=1,Type=Integer,Description="Dosage count of reference alleles from updog, where 0 = homozygous alternate">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=RA,Number=1,Type=Integer,Description="Reference allele read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    '##FORMAT=<ID=MPP,Number=1,Type=Float,Description="Maximum posterior probability for that dosage call from updog">',
    updog_meta,
    bigr_meta
  )

  #Get the total depth and total ref and total alt depths for each SNP across all samples for the VCF
  depth_df <- mout$inddf %>%
    group_by(snp) %>%
    summarize(
      total_ref = sum(ref),
      total_size = sum(size),
      total_alt = sum(size)-sum(ref)
    )

  #Make sure it is sorted
  depth_df <- depth_df %>%
    arrange(match(snp, mout$snpdf$snp))

  #CHROM and POS from SNP ID
  new_df <- mout$snpdf %>%
    separate(snp, into = c("CHROM", "POS"), sep = "_") %>%
    select(CHROM, POS)
  new_df$POS <- sub("^0+", "", new_df$POS)

  #Make the VCF df
  vcf_df <- data.frame(
    CHROM = new_df$CHROM,
    POS = new_df$POS,
    ID = mout$snpdf$snp,
    REF = "A",
    ALT = "B",
    QUAL = ".",
    FILTER = ".",
    INFO = NA,
    FORMAT = NA
  )

  #Add the INFO column for each SNP
  vcf_df$INFO <- paste0("DP=",depth_df$total_size,";",
                        "ADS=",depth_df$total_ref,",",depth_df$total_alt,";",
                        "BIAS=",mout$snpdf$bias,";",
                        "OD=",mout$snpdf$od,";",
                        "PMC=",mout$snpdf$prop_mis)

  #Add the FORMAT label for each SNP
  vcf_df$FORMAT <- paste("GT","UD","DP","RA","AD","MPP",sep=":")

  #Convert genotypes from dosage to gt
  convert_dosage_to_genotype <- function(dosage, ploidy) {
    # Handle missing values
    if (is.na(dosage)) {
      return("./.")
    }

    # Generate the genotype string based on the dosage and ploidy
    # Updog uses the ref counts, which is not typical, so this corrects it
    ref_count <- dosage
    alt_count <- ploidy - dosage
    genotype <- paste(c(rep("0", ref_count), rep("1", alt_count)), collapse = "/")

    return(genotype)
  }

  #creating a temporary df for the converted genotypes
  geno_df <- mout$inddf[,c("snp","geno")] %>%
    mutate(genotype = sapply(geno, convert_dosage_to_genotype, ploidy = as.numeric(ploidy)))

  #Format the mout FORMAT info
  format_df <- data.frame(
    snp = mout$inddf$snp,
    ind = mout$inddf$ind,
    format = paste0(geno_df$genotype,":",
                    mout$inddf$geno,":",
                    mout$inddf$size,":",
                    mout$inddf$ref,":",
                    mout$inddf$ref,",",(mout$inddf$size-mout$inddf$ref),":",
                    mout$inddf$maxpostprob)
  )
  #Make sample the columns
  format_wide <- format_df %>%
    pivot_wider(names_from = ind, values_from = format)

  #Add samples to vcf df
  vcf_df <- cbind(vcf_df,format_wide[,-1])

  # Add # to the CHROM column name
  colnames(vcf_df)[1] <- "#CHROM"

  if(!compress){
    # Write the header to the file
    writeLines(vcf_header, con = output.file)

    # Append the dataframe to the file in tab-separated format
    suppressWarnings(
      write.table(vcf_df, file = output.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
    )
  } else {

    gz <- gzfile(paste0(output.file, ".gz"), "w")
    write(vcf_header, gz)
    header <- colnames(vcf_df)
    header <- paste(header, collapse="\t")
    write(header, gz)
    close(gz)
    vcfR:::.write_vcf_body(fix = as.matrix(vcf_df[,1:8]), gt = as.matrix(vcf_df[,9:ncol(vcf_df)]), filename = paste0(output.file, ".gz"))
  }
}
