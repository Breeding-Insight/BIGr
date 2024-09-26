library(vcfR)
library(Biostrings)

report_csv <- read.csv("DSp23-8876_MADC_rename_LSU_filter_miss.csv")
result_csv <- read.csv("DSp23-8876_MADC_rename_LSU_filter_miss_snps.csv")
result_vcf <- read.vcfR("DSp23-8876_MADC_rename_LSU_filter_miss_snps.vcf")

# Comparing when botloci
head(report_csv[,1:4],7)
head(result_csv[,1:7],13)
head(result_vcf@fix,7)
head(result_vcf@gt[,1:2],7)

# Comparing when not botloci
report_csv[grep("Chr01_000239479", report_csv$AlleleID),1:4]
result_csv[grep("Chr01_000239479", result_csv$AlleleID),1:7]
result_vcf@fix[grep("Chr01_000239479", result_vcf@fix[,3]),]
result_vcf@gt[grep("Chr01_000239479", result_vcf@fix[,3]),1:2]


botloci <- read.csv("sweetpotato_20K_SNPset_f180bp_forDArT_3K_f180bp.botloci", header = F)
botloci <- botloci$V1

had_seq <- read.table("sweetpotato_allele_db_v001.fa", header = F)

report <- read.csv("DSp23-8876_MADC_rename_LSU_filter_miss.csv")

had_seq <- get_ref_alt_hap_seq(had_seq)

loop_though_dartag_report <- function(report, botloci, hap_seq, verbose = TRUE){}

nsamples <- ncol(report) - 3
new.file <- data.frame("AlleleID" = report[,1],
                       "Chromosome" = NA, "SNP_position_in_Genome" = NA,
                       "Ref" = NA, "Alt" =NA,
                       "CloneID" = report[,2], report[,3:ncol(report)])

by_cloneID <- split.data.frame(new.file, new.file$CloneID)

ref_t <- alt_t <- vector()
for(i in 1:length(by_cloneID)){

  # Add ref and alt with absent
  cloneID <- by_cloneID[[i]]$CloneID[1]
  ref <- paste0(cloneID, "|Ref_0001")
  alt <- paste0(cloneID, "|Alt_0002")

  ref_t <- c(ref_t,!ref %in% by_cloneID[[i]]$AlleleID)
  if(!ref %in% by_cloneID[[i]]$AlleleID){
    empty_allele <- c(ref, NA, NA, NA,NA, cloneID, had_seq[which(had_seq$AlleleID == ref),2], rep(0, nsamples))
    by_cloneID[[i]] <- rbind(by_cloneID[[i]], empty_allele)
  }

  alt_t <- c(alt_t,!alt %in% by_cloneID[[i]]$AlleleID)
  if(!alt %in% by_cloneID[[i]]$AlleleID){
    empty_allele <- c(alt, NA, NA, NA,NA, cloneID, had_seq[which(had_seq$AlleleID == alt),2], rep(0, nsamples))
    by_cloneID[[i]] <- rbind(by_cloneID[[i]], empty_allele)
  }

  # Get SNP positions based on the reference
  compare <- function(){}
  by_cloneID[[i]]$AlleleID

  # Target SNP at the position in the ID
  ref_seq <- by_cloneID[[i]]$AlleleSequence[grep("Ref_0001$",by_cloneID[[i]]$AlleleID)]
  alt_seq <- by_cloneID[[i]]$AlleleSequence[grep("Alt_0002$",by_cloneID[[i]]$AlleleID)]
  pos_target <- as.numeric(sapply(strsplit(by_cloneID[[i]]$CloneID[1], "_"), function(x) x[length(x)]))

  sigma <- nucleotideSubstitutionMatrix(match = 1, mismatch = -0.5, baseOnly = FALSE)
  align <- pairwiseAlignment(ref_seq,
                             alt_seq,
                             substitutionMatrix = sigma,gapOpening=-1.5, gapExtension=-0.1, type = "global")

  pos_target_idx <- align@pattern@mismatch@unlistData
  pos_target_alt_idx <- align@subject@mismatch@unlistData
  ref_base <- substring(ref_seq, align@pattern@mismatch@unlistData, align@pattern@mismatch@unlistData)
  alt_base <- substring(alt_seq, align@subject@mismatch@unlistData, align@subject@mismatch@unlistData)

  pos_start <- pos_target - pos_target_idx

  altMatch_seq <- by_cloneID[[i]][grep("AltMatch",by_cloneID[[i]]$AlleleID),]
  for(j in 1:nrow(altMatch_seq)){
    align <- pairwiseAlignment(ref_seq,
                               altMatch_seq[j,]$AlleleSequence,
                               substitutionMatrix = sigma,gapOpening=-1.5, gapExtension=-0.1, type = "global")
    pos_idx <- align@pattern@mismatch@unlistData
    pos_idx <- pos_idx[-which(pos_idx == pos_target_idx)]
    ref_base <- substring(ref_seq, pos_idx, pos_idx)
    pos_idx <- align@subject@mismatch@unlistData
    pos_idx <- pos_idx[-which(pos_idx == pos_target_alt_idx)]
    alt_base <- substring(altMatch_seq[j,]$AlleleSequence, pos_idx, pos_idx)
  }

  refMatch_seq <- by_cloneID[[i]][grep("RefMatch",by_cloneID[[i]]$AlleleID),]
  for(j in 1:nrow(refMatch_seq)){
    align <- pairwiseAlignment(ref_seq,
                               refMatch_seq[j,]$AlleleSequence,
                               substitutionMatrix = sigma,gapOpening=-1.5, gapExtension=-0.1, type = "global")
    pos_idx <- align@pattern@mismatch@unlistData
    ref_base <- substring(ref_seq, pos_idx, pos_idx)
    pos_idx <- align@subject@mismatch@unlistData
    alt_base <- substring(altMatch_seq[j,]$AlleleSequence, pos_idx, pos_idx)
  }

  # Stopped here - store the results in the table and keep going

}

if(verbose){
  cat("Reference alleles sequence added:", sum(ref_t))
  cat("Alternative alleles sequence added:", sum(alt_t))
}

with_ref_alt <- do.call(rbind, by_cloneID)


# Converts the fasta to a data.frame with first column the AlleleID and and second the AlleleSequence
# The function will work even if the sequence is split in multiple lines
get_ref_alt_hap_seq <- function(hap_seq){
  headers <- had_seq$V1[grep(">",had_seq$V1)]
  headers <- gsub(">", "", headers)

  seqs <- had_seq$V1[-grep(">",had_seq$V1)]

  times <- diff(grep(">",had_seq$V1))-1
  last <- length(had_seq$V1) - grep(">",had_seq$V1)[length(grep(">",had_seq$V1))]
  times <- c(times, last)
  index <- vector()
  for(i in 1:length(headers)){
    index <- c(index, rep(i, times[i]))
  }

  seqs <- split(seqs, index)
  seqs <- sapply(seqs, function(x) paste0(x, collapse = ""))

  hap_seq <- data.frame(AlleleID = headers, AlleleSequence = seqs)
  return(hap_seq)
}



