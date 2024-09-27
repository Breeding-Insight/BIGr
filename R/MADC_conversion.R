library(vcfR)
library(pwalign)
library(Biostrings)
library(parallel)

# report_csv <- read.csv("DSp23-8876_MADC_rename_LSU_filter_miss.csv")
# result_csv <- read.csv("DSp23-8876_MADC_rename_LSU_filter_miss_snps.csv")
# result_vcf <- read.vcfR("DSp23-8876_MADC_rename_LSU_filter_miss_snps.vcf")
#
# # Comparing when botloci
# head(report_csv[,1:4],7)
# head(result_csv[,1:7],13)
# head(result_vcf@fix,7)
# head(result_vcf@gt[,1:2],7)
#
# # Comparing when not botloci
# report_csv[grep("Chr01_000239479", report_csv$AlleleID),1:4]
# result_csv[grep("Chr01_000239479", result_csv$AlleleID),1:7]
# result_vcf@fix[grep("Chr01_000239479", result_vcf@fix[,3]),]
# result_vcf@gt[grep("Chr01_000239479", result_vcf@fix[,3]),1:2]
#
# # Inputs
# botloci <- read.csv("sweetpotato_20K_SNPset_f180bp_forDArT_3K_f180bp.botloci", header = F)
# botloci <- botloci$V1
# had_seq <- read.table("sweetpotato_allele_db_v001.fa", header = F)
# report <- read.csv("DSp23-8876_MADC_rename_LSU_filter_miss.csv")
# n.cores <- 5
#
# my_results_csv <- loop_though_dartag_report(report, botloci, hap_seq, n.cores=5, verbose = TRUE)
#
# # Compare with python result
# dim(result_csv)
# result_csv <- result_csv[order(result_csv$Chromosome, result_csv$SNP_position_in_Genome),]
# my_results_csv$SNP_position_in_Genome <- as.numeric(my_results_csv$SNP_position_in_Genome)
# my_results_csv <- my_results_csv[order(my_results_csv$Chromosome, my_results_csv$SNP_position_in_Genome),]
#
# sum(!result_csv$AlleleID %in% my_results_csv$AlleleID)
# dim(my_results_csv)[1] - dim(result_csv)[1]


#' Include SNP_position_in_Genome, Ref, and Alt information
#'
loop_though_dartag_report <- function(report, botloci, hap_seq, n.cores=1, verbose = TRUE){

  had_seq <- get_ref_alt_hap_seq(had_seq)

  nsamples <- ncol(report) - 3
  new.file <- data.frame("AlleleID" = report[,1],
                         "Chromosome" = NA, "SNP_position_in_Genome" = NA,
                         "Ref" = NA, "Alt" =NA,
                         "CloneID" = report[,2], report[,3:ncol(report)])

  by_cloneID <- split.data.frame(new.file, new.file$CloneID)

  clust <- makeCluster(n.cores)
  clusterExport(clust, c("had_seq","add_ref_alt"))
  add_ref_alt_results <- parLapply(clust, by_cloneID, function(x) add_ref_alt(x, had_seq))
  stopCluster(clust)

  updated_by_cloneID <- lapply(add_ref_alt_results, "[[",1)
  ref_index <- sapply(add_ref_alt_results, "[[",2)
  alt_index <- sapply(add_ref_alt_results, "[[",3)

  if(verbose){
    cat("The Ref_0001 sequence had to be added for:", sum(ref_index),"tags\n")
    cat("The Alt_0002 sequence had to be added for:", sum(alt_index),"tags\n")
  }

  clust <- makeCluster(n.cores)
  #clusterExport(clust, c("botloci","compare"))
  clusterExport(clust, c("botloci","compare", "DNAString", "reverseComplement", "pairwiseAlignment", "nucleotideSubstitutionMatrix"))
  compare_results <- parLapply(clust, updated_by_cloneID, function(x) compare(x, botloci))
  stopCluster(clust)

  my_results_csv <- lapply(compare_results, "[[", 1)
  my_results_csv <- do.call(rbind, my_results_csv)

  rm_score <- sapply(compare_results, "[[", 2)
  rm_score <- unlist(rm_score)
  rm_N <- sapply(compare_results, "[[", 3)
  rm_N <- unlist(rm_N)

  if(verbose){
    cat("Number of tags removed because of low alignment score:", length(rm_score),"tags\n")
    cat("Number of tags removed because of N in the alternative sequence:", length(rm_N),"tags\n")
  }

  rownames(my_results_csv) <- NULL
  return(my_results_csv)
}

#' Check if Ref_0001 and Alt_0002 tags are present, if not, add them from the had_seq input
#'
#' @param one_tag
#' @param had_seq
#'
add_ref_alt <- function(one_tag, had_seq) {

  # Add ref and alt
  cloneID <- one_tag$CloneID[1]
  ref <- paste0(cloneID, "|Ref_0001")
  alt <- paste0(cloneID, "|Alt_0002")

  # Precompute ref and alt indices
  ref_index <- alt_index <- 0

  # Use match() for faster lookup in had_seq
  ref_row <- match(ref, had_seq$AlleleID)
  alt_row <- match(alt, had_seq$AlleleID)

  # Initialize a list to store potential new rows to add
  new_rows <- list()

  if (!ref %in% one_tag$AlleleID) {
    ref_index <- 1
    # Only extract the relevant sequence data once from had_seq
    ref_seq <- if (!is.na(ref_row)) had_seq[ref_row, 2] else NA
    empty_allele <- c(ref, rep(NA, 4), cloneID, ref_seq, rep(0, nsamples))
    new_rows[[length(new_rows) + 1]] <- empty_allele
  }

  if (!alt %in% one_tag$AlleleID) {
    alt_index <- 1
    alt_seq <- if (!is.na(alt_row)) had_seq[alt_row, 2] else NA
    empty_allele <- c(alt, rep(NA, 4), cloneID, alt_seq, rep(0, nsamples))
    new_rows[[length(new_rows) + 1]] <- empty_allele
  }

  # Only rbind once if new rows were added
  if (length(new_rows) > 0) {
    one_tag <- rbind(one_tag, do.call(rbind, new_rows))
  }

  return(list(one_tag, ref_index, alt_index))
}


#' Get SNP positions, reference and alternative alleles based on the reference
#' Align alternatives to reference and discard low score alignment tags
#' Discard tags if alternative in the target locus is N
#' Do the complement reverse if cloneID present in the botloci vector
#'
#' @param one_tag
#' @param botloci
#'
#' @importFrom Biostrings DNAString reverseComplement
#' @importFrom pwalign pairwiseAlignment nucleotideSubstitutionMatrix
#'
compare <- function(one_tag, botloci){

  cloneID <- one_tag$CloneID[1]
  isBotLoci <- cloneID %in% botloci
  # If marker is present in the botloc list, get the reverse complement sequence
  if(isBotLoci) one_tag$AlleleSequence <- sapply(one_tag$AlleleSequence, function(x) as.character(reverseComplement(DNAString(x))))

  chr <- sapply(strsplit(cloneID, "_"), function(x) x[-length(x)])
  # Target SNP at the position in the ID
  ref_seq <- one_tag$AlleleSequence[grep("Ref_0001$",one_tag$AlleleID)]
  alt_seq <- one_tag$AlleleSequence[grep("Alt_0002$",one_tag$AlleleID)]
  pos_target <- as.numeric(tail(strsplit(cloneID, "_")[[1]], 1))

  sigma <- nucleotideSubstitutionMatrix(match = 1, mismatch = -0.5, baseOnly = FALSE) # baseOnly = FALSE avoid breaking when N is present (they will be filtered after))
  align <- pairwiseAlignment(ref_seq,
                             alt_seq,
                             substitutionMatrix = sigma,gapOpening=-1.4, gapExtension=-0.1, type = "global")

  # The score is a bit different from the python script despite same weights
  if(align@score > 40){ # if score for the target sequence is smaller than 40, the tag will be discarted
    pos_target_idx <- align@pattern@mismatch@unlistData
    pos_target_alt_idx <- align@subject@mismatch@unlistData
    ref_base <- substring(ref_seq, align@pattern@mismatch@unlistData, align@pattern@mismatch@unlistData)
    alt_base <- substring(alt_seq, align@subject@mismatch@unlistData, align@subject@mismatch@unlistData)

    # If target alternative have N, discard whole tag
    if(all(alt_base %in% c("A", "T", "C", "G"))) {

      # Update with new info - always only one polymorphism, if there are more than one, not sure which is the target
      update_tag <- one_tag[grep("Ref_0001$",one_tag$AlleleID),]
      update_tag[,2:5] <- c(chr, pos_target, ref_base, alt_base)
      update_tag_temp <- one_tag[grep("Alt_0002$",one_tag$AlleleID),]
      update_tag_temp[,2:5] <- c(chr, pos_target, ref_base, alt_base)
      update_tag <- rbind(update_tag, update_tag_temp)

      Match_seq <- one_tag[grep("Match",one_tag$AlleleID),]
      if(nrow(Match_seq) >0){
        for(j in 1:nrow(Match_seq)){
          align <- pairwiseAlignment(ref_seq,
                                     Match_seq[j,]$AlleleSequence,
                                     substitutionMatrix = sigma,gapOpening=-1.4, gapExtension=-0.1, type = "global")
          pos_ref_idx <- align@pattern@mismatch@unlistData
          rm_target <- which(pos_ref_idx == pos_target_idx)                 # remove target position when is AltMatch
          if(length(rm_target) >0) pos_ref_idx <- pos_ref_idx[-rm_target]
          # Cases found where the AltMatch is another alternative for the target SNP - they are discarted
          if(length(pos_ref_idx) >0){
            ref_base <- substring(ref_seq, pos_ref_idx, pos_ref_idx)
            pos_alt_idx <- align@subject@mismatch@unlistData                 # If there are indels, the position in the alternative is not the same as the reference
            if(length(rm_target) >0) pos_alt_idx <- pos_alt_idx[-rm_target]   # remove target position when is AltMatch - but the order in the sequence is the same
            alt_base <- substring(Match_seq[j,]$AlleleSequence, pos_alt_idx, pos_alt_idx)

            # If Match sequences have N, do not consider as polymorphism
            if(any(!alt_base %in% c("A", "T", "C", "G"))) {
              alt_base <- alt_base[-which(!alt_base %in% c("A", "T", "C", "G"))]
              ref_base <- ref_base[-which(!alt_base %in% c("A", "T", "C", "G"))]
            }

            if(length(alt_base) >0){ # If the N is the only polymorphis found, the Match tag will be discarted
              # The reported position is always on reference
              # If present in the botloc, the position will be before target, if not, after
              pos <- pos_target - (pos_target_idx - pos_ref_idx)
              # if(isBotLoci) pos <- pos_target - (pos_target_idx - pos_ref_idx) else
              #   pos <- pos_target + (pos_target_idx - pos_ref_idx)

              # Sometimes there are more than one polymorphism in the sequence, we need to add rows to the table
              update_tag_temp <- one_tag[grep("Match",one_tag$AlleleID)[j],][rep(1, length(alt_base)), ]

              update_tag_temp$Chromosome <- chr
              update_tag_temp$SNP_position_in_Genome <- pos
              update_tag_temp$Ref <- ref_base
              update_tag_temp$Alt <- alt_base

              update_tag <- rbind(update_tag, update_tag_temp)
            }
          }
        }
      }
      return(list(update_tag = update_tag, # updated data.frame, NULL if discarted
                  rm_score = NULL,         # cloneID if removed because of low alignment score, NULL if kept
                  rm_N = NULL))            # cloneID if removed because of N in the target alternative, NULL if kept
    } else {
      return(list(update_tag = NULL,
                  rm_score = NULL,
                  rm_N = cloneID))
    }
  } else{
    return(list(update_tag = NULL,
                rm_score = cloneID,
                rm_N = NULL))
  }
}

#' Converts the fasta to a data.frame with first column the AlleleID and and second the AlleleSequence
#' The function will work even if the sequence is split in multiple lines
#'
#' @param hap_seq
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



