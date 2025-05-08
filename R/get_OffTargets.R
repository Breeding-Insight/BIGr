
#' Converts MADC file to VCF recovering target and off-target SNPs
#' @param madc path to MADC file
#' @param botloci path to file containing the target IDs that were designed in the bottom strand
#' @param hap_seq path to haplotype DB FASTA file
#' @param rm_multiallelic_SNP logical. If TRUE, SNP with more than one alternative base will be removed. If FALSE, check `multiallelic_SNP_dp_thr` specs
#' @param multiallelic_SNP_dp_thr numerical. If `rm_multiallelic_SNP` is FALSE, set a minimum depth by tag threshold `multiallelic_SNP_dp_thr` combined
#' with minimum number of samples `multiallelic_SNP_sample_thr` to eliminate low frequency SNP allele. If the threshold does not eliminate the multiallelic
#' aspect of the marker, the marker is discarded. This is likely to happen to paralogous sites.
#' @param multiallelic_SNP_sample_thr numerical. If `rm_multiallelic_SNP` is FALSE, set a minimum depth by tag threshold combined with minimum number of
#' samples `multiallelic_SNP_sample_thr` to eliminate low frequency SNP allele. If the threshold does not eliminate the multiallelic aspect of the marker,
#' the marker is discarded. This is likely to happen to paralogous sites.
#' @param n.cores number of cores to be used in the parallelization
#' @param out_vcf output VCF file name
#' @param verbose print metrics on the console
#'
#' @importFrom utils packageVersion read.csv
#' @import vcfR
#'
#' @export
get_OffTargets <- function(madc = NULL,
                           botloci = NULL,
                           hap_seq = NULL,
                           n.cores = 5,
                           rm_multiallelic_SNP = FALSE,
                           multiallelic_SNP_dp_thr = 0,
                           multiallelic_SNP_sample_thr = 0,
                           out_vcf = NULL,
                           verbose = TRUE){

  bigr_meta <- paste0('##BIGrCommandLine.getOffTargets=<ID=getOffTargets,Version="',
                      packageVersion("BIGr"), '",Data="',
                      Sys.time(),'", CommandLine="> getOffTargets(',deparse(substitute(madc)),', ',
                      "botloci= ", botloci, ', ',
                      "hap_seq= ", hap_seq, ', ',
                      "n.cores= ", n.cores, ', ',
                      "rm_multiallelic_SNP= ", rm_multiallelic_SNP, ', ',
                      "multiallelic_SNP_dp_thr= ", multiallelic_SNP_dp_thr, ', ',
                      "multiallelic_SNP_sample_thr= ", multiallelic_SNP_sample_thr, ', ',
                      "out_vcf= ", out_vcf, ', ',
                      "verbose= ", verbose,')">')

  report <- read.csv(madc)
  botloci <- read.csv(botloci, header = F)
  hap_seq <- read.table(hap_seq, header = F)

  my_results_csv <- loop_though_dartag_report(report,
                                              botloci,
                                              hap_seq,
                                              n.cores=n.cores,
                                              verbose = verbose)

  vcf_body <- create_VCF_body(csv = my_results_csv,
                              n.cores = n.cores,
                              rm_multiallelic_SNP = rm_multiallelic_SNP,
                              multiallelic_SNP_dp_thr = multiallelic_SNP_dp_thr,
                              multiallelic_SNP_sample_thr = multiallelic_SNP_sample_thr,
                              verbose = verbose)

  #Make a header separate from the dataframe
  vcf_header <- c(
    "##fileformat=VCFv4.3",
    "##reference=NA",
    "##contig=<ID=NA,length=NA>",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##INFO=<ID=ADS,Number=R,Type=Integer,Description="Depths for the ref and each alt allele in the order listed">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=RA,Number=1,Type=Integer,Description="Reference allele read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    bigr_meta
  )

  vcf_term <- sapply(strsplit(out_vcf, "[.]"), function(x) x[length(x)])
  if(length(vcf_term) != 0) if(vcf_term != "vcf") out_vcf <- paste0(out_vcf,".vcf")

  writeLines(vcf_header, con = out_vcf)
  suppressWarnings(
    write.table(vcf_body, file = out_vcf, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
  )
}

#' Include SNP_position_in_Genome, Ref, and Alt information
#'
#' @param report MADC file
#' @param botloci file containing the target IDs that were designed in the bottom strand
#' @param hap_seq haplotype DB fasta file
#' @param n.cores number of cores to be used in the parallelization
#' @param verbose print metrics on the console
#'
#' @import parallel
#'
#' @noRd
loop_though_dartag_report <- function(report, botloci, hap_seq, n.cores=1, verbose = TRUE){

  hap_seq <- get_ref_alt_hap_seq(hap_seq)

  nsamples <- ncol(report) - 3
  new.file <- data.frame("AlleleID" = report[,1],
                         "Chromosome" = NA, "SNP_position_in_Genome" = NA,
                         "Ref" = NA, "Alt" =NA,
                         "CloneID" = report[,2], report[,3:ncol(report)])

  by_cloneID <- split.data.frame(new.file, new.file$CloneID)

  clust <- makeCluster(n.cores)
  #clusterExport(clust, c("hap_seq","add_ref_alt"))
  add_ref_alt_results <- parLapply(clust, by_cloneID, function(x) add_ref_alt(x, hap_seq, nsamples))
  stopCluster(clust)

  updated_by_cloneID <- lapply(add_ref_alt_results, "[[",1)
  ref_index <- sapply(add_ref_alt_results, "[[",2)
  alt_index <- sapply(add_ref_alt_results, "[[",3)

  if(verbose){
    cat("The Ref_0001 sequence had to be added for:", sum(ref_index),"tags\n")
    cat("The Alt_0002 sequence had to be added for:", sum(alt_index),"tags\n")
  }

  clust <- makeCluster(n.cores)
  #clusterExport(clust, c("botloci", "compare", "nucleotideSubstitutionMatrix", "pairwiseAlignment", "DNAString", "reverseComplement"))
  #clusterExport(clust, c("botloci"))
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

#' Check if Ref_0001 and Alt_0002 tags are present, if not, add them from the hap_seq input. Function made for parallelization.
#'
#' @param one_tag madc file split by tag
#' @param hap_seq haplotype DB
#' @param nsamples number of samples
#'
#' @noRd
add_ref_alt <- function(one_tag, hap_seq, nsamples) {

  # Add ref and alt
  cloneID <- one_tag$CloneID[1]
  ref <- paste0(cloneID, "|Ref_0001")
  alt <- paste0(cloneID, "|Alt_0002")

  # Precompute ref and alt indices
  ref_index <- alt_index <- 0

  # Use match() for faster lookup in hap_seq
  ref_row <- match(ref, hap_seq$AlleleID)
  alt_row <- match(alt, hap_seq$AlleleID)

  # Initialize a list to store potential new rows to add
  new_rows <- list()
  if (!ref %in% one_tag$AlleleID) {
    ref_index <- 1
    # Only extract the relevant sequence data once from hap_seq
    ref_seq <- if (!is.na(ref_row)) hap_seq[ref_row, 2] else NA
    empty_allele <- c(ref, rep(NA, 4), cloneID, ref_seq, rep(0, nsamples))
    new_rows[[length(new_rows) + 1]] <- empty_allele
  }

  if (!alt %in% one_tag$AlleleID) {
    alt_index <- 1
    alt_seq <- if (!is.na(alt_row)) hap_seq[alt_row, 2] else NA
    empty_allele <- c(alt, rep(NA, 4), cloneID, alt_seq, rep(0, nsamples))
    new_rows[[length(new_rows) + 1]] <- empty_allele
  }

  # Only rbind once if new rows were added
  if (length(new_rows) > 0) {
    if(length(new_rows) == 1) one_tag <- rbind(one_tag, new_rows[[1]]) else  {
      new_matrix <- do.call(rbind, new_rows)
      colnames(new_matrix) <- colnames(one_tag)
      one_tag <- rbind(one_tag, new_matrix)
    }
  }

  return(list(one_tag, ref_index, alt_index))
}


#' Get SNP positions, reference and alternative alleles based on the reference
#' Align alternatives to reference and discard low score alignment tags
#' Discard tags if alternative in the target locus is N
#' Do the complement reverse if cloneID present in the botloci vector
#'
#' @param one_tag madc file split by tag
#' @param botloci file containing the target IDs that were designed in the bottom strand
#'
#' @importFrom Biostrings DNAString reverseComplement
#' @importFrom pwalign pairwiseAlignment nucleotideSubstitutionMatrix
#'
#' @noRd
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
#' @param hap_seq haplotype db
#'
#' @noRd
get_ref_alt_hap_seq <- function(hap_seq){
  headers <- hap_seq$V1[grep(">",hap_seq$V1)]
  headers <- gsub(">", "", headers)

  seqs <- hap_seq$V1[-grep(">",hap_seq$V1)]

  times <- diff(grep(">",hap_seq$V1))-1
  last <- length(hap_seq$V1) - grep(">",hap_seq$V1)[length(grep(">",hap_seq$V1))]
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


#' Creates VCF body from CSV generated by loop_though_dartag_report
#'
#' @param csv CSV file generated by loop_though_dartag_report
#' @param rm_multiallelic_SNP logical. If TRUE, SNP with more than one alternative base will be removed. If FALSE, check `multiallelic_SNP_dp_thr` specs
#' @param multiallelic_SNP_dp_thr numerical. If `rm_multiallelic_SNP` is FALSE, set a minimum
#' depth by tag threshold combined with minimum number of samples (`multiallelic_SNP_dp_thr` + `multiallelic_SNP_sample_thr`)
#' to eliminate low frequency SNP allele. If the threshold does not eliminate the multiallelic aspect of the marker, the marker
#' is discarded. This is likely to happen to paralogous sites.
#' @param multiallelic_SNP_sample_thr numerical. If `rm_multiallelic_SNP` is FALSE, set a minimum depth by tag threshold `multiallelic_SNP_dp_thr` combined
#' with minimum number of samples `multiallelic_SNP_sample_thr` to eliminate low frequency SNP allele. If the threshold does not eliminate the multiallelic
#' aspect of the marker, the marker is discarded. This is likely to happen to paralogous sites.
#' @param n.cores number of cores to be used in the parallelization
#' @param verbose print metrics on the console
#'
#' @import parallel
#'
#' @noRd
create_VCF_body <- function(csv,
                            rm_multiallelic_SNP = TRUE,
                            multiallelic_SNP_dp_thr = 2,
                            multiallelic_SNP_sample_thr = 10,
                            n.cores = 1,
                            verbose = TRUE){

  # Make sure counts are numeric
  csv[,-c(1:7)] <- apply(csv[,-c(1:7)], 2, as.numeric)

  #csv <- my_results_csv
  cloneID <- split.data.frame(csv, csv$CloneID)

  clust <- makeCluster(n.cores)
  #clusterExport(clust, c("merge_counts","rm_multiallelic_SNP", "multiallelic_SNP_dp_thr", "multiallelic_SNP_sample_thr"))
  vcf_tag_list <- parLapply(clust, cloneID, function(x) merge_counts(x, rm_multiallelic_SNP, multiallelic_SNP_dp_thr, multiallelic_SNP_sample_thr))
  stopCluster(clust)

  vcf_tag_list1 <- lapply(vcf_tag_list, "[[", 1)
  rm_mks <- sapply(vcf_tag_list, "[[" ,2)

  if(verbose) print(paste("SNP removed because presented more than one allele:", sum(rm_mks)))

  for(i in 1:length(vcf_tag_list1)) {
    if(is.vector(vcf_tag_list1[[i]])) {
      vcf_tag_list1[[i]] <- c(target = names(vcf_tag_list1)[i],vcf_tag_list1[[i]])
    } else  vcf_tag_list1[[i]] <- cbind(target = names(vcf_tag_list1)[i],vcf_tag_list1[[i]])
  }

  vcf_body <- do.call("rbind",vcf_tag_list1)
  vcf_body <- as.data.frame(vcf_body)
  vcf_body$V3 <- as.numeric(vcf_body$V3)
  rownames(vcf_body) <- NULL

  vcf_body$target <- paste0(sapply(strsplit(vcf_body$target, "_"),"[[", 1), "_",as.numeric(sapply(strsplit(vcf_body$target, "_"),"[[", 2)))

  # Dealing with repeated positions
  # discard the ones that are not the target and keep only the first if all are off-targets
  if(length(which(duplicated(vcf_body[,3]))) > 0){
    repeated <- vcf_body[which(duplicated(vcf_body[,3])), 4]

    if(verbose) print(paste("Different primers pair capture same SNP positions in", length(repeated), "locations. The repeated were discarded."))

    repeated_tab <- vcf_body[which(vcf_body[,4] %in% repeated),]
    vcf_body_new <- vcf_body[-which(vcf_body[,4] %in% repeated),]

    # discard the off-targets that overlap with the target position
    repeated_tab_stay1 <- repeated_tab[which(repeated_tab$target == repeated_tab$V4),]
    repeated_tab <- repeated_tab[-which(repeated_tab$V4 %in% repeated_tab_stay1$V4),]

    # Keep only the first appearance of the repeated off-targets
    repeated_tab <- repeated_tab[order(repeated_tab$target, repeated_tab$V4),]
    repeated_tab_stay2 <- repeated_tab[-which(duplicated(repeated_tab$V4)),]

    repeated_tab_stay <- rbind(repeated_tab_stay1, repeated_tab_stay2)

    vcf_body_new <- rbind(vcf_body_new, repeated_tab_stay)
  } else vcf_body_new <- vcf_body

  vcf_body_new <- vcf_body_new[,-1]

  colnames(vcf_body_new) <- c("#CHROM", "POS", "ID", "REF", "ALT","QUAL", "FILTER", "INFO","FORMAT", colnames(csv)[-c(1:7)])
  vcf_body_new <- vcf_body_new[order(vcf_body_new[,1], vcf_body_new[,2]),]

  return(vcf_body_new)
}


#' Function made for parallelization of create_VCF_body function
#'
#' @param cloneID_unit one item of csv file split by cloneID
#' @param rm_multiallelic_SNP logical. If TRUE, SNP with more than one alternative base will be removed. If FALSE, check `multiallelic_SNP_dp_thr` specs
#' @param multiallelic_SNP_dp_thr numerical. If `rm_multiallelic_SNP` is FALSE, set a minimum depth by tag threshold combined with minimum number of samples
#' `multiallelic_SNP_sample_thr` to eliminate low frequency SNP allele. If the threshold does not eliminate the multiallelic aspect of the marker, the marker
#' is discarded. This is likely to happen to paralogous sites.
#' @param multiallelic_SNP_sample_thr numerical. If `rm_multiallelic_SNP` is FALSE, set a minimum depth by tag threshold `multiallelic_SNP_dp_thr` combined
#' with minimum number of samples `multiallelic_SNP_sample_thr` to eliminate low frequency SNP allele. If the threshold does not eliminate the multiallelic
#' aspect of the marker, the marker is discarded. This is likely to happen to paralogous sites.
#'
#' @noRd
merge_counts <- function(cloneID_unit, rm_multiallelic_SNP = FALSE, multiallelic_SNP_dp_thr = 0,  multiallelic_SNP_sample_thr = 0){

  #Get counts for target SNP
  rm <- 0
  RefTag <- apply(cloneID_unit[which(grepl("Ref", cloneID_unit$AlleleID) & !duplicated(cloneID_unit$AlleleID)),-c(1:7)], 2, sum)
  AltTag <- apply(cloneID_unit[which(grepl("Alt", cloneID_unit$AlleleID) & !duplicated(cloneID_unit$AlleleID)),-c(1:7)], 2, sum)
  tab_counts <- paste0(RefTag + AltTag, ":", RefTag, ":", RefTag, ",", AltTag)

  info <- cloneID_unit[grep("Ref_", cloneID_unit$AlleleID),]
  info <- c(info$Chromosome,
            info$SNP_position_in_Genome,
            paste0(info$Chromosome, "_", info$SNP_position_in_Genome),
            info$Ref, info$Alt, ".", ".", paste0("DP=", sum(c(RefTag, AltTag)),";",
                                                 "ADS=",sum(RefTag),",",sum(AltTag)), "DP:RA:AD")

  vcf_tag <- c(info, tab_counts)

  # Check if there are more than one alternative allele by loci
  off_tag <- cloneID_unit[-which(grepl("Ref_", cloneID_unit$AlleleID) | grepl("Alt_", cloneID_unit$AlleleID)),]
  if(nrow(off_tag)){ # If there are off target SNP
    by_pos <- split.data.frame(off_tag, off_tag$SNP_position_in_Genome)

    ## Possibility to filter alt alleles with low read
    # The filtering is done after the counts are collapsed so all reads can be counted
    for(i in 1:length(by_pos)){
      alleles <- unique(by_pos[[i]]$AlleleID)

      if(length(unique(by_pos[[i]]$Alt)) > 1){ # If SNP is multiallelic

        if(rm_multiallelic_SNP){ # option to remove multiallelics
          rm <- rm + 1
          next()
        } else if(multiallelic_SNP_dp_thr > 0 & multiallelic_SNP_sample_thr > 0){ # If not removed, user can set threshold to remove low frequency alleles
          rm.idx <- which(apply(by_pos[[i]][,-c(1:7)], 1, function(x) sum(x > multiallelic_SNP_dp_thr) < multiallelic_SNP_sample_thr))
          if(length(rm.idx))  up_by_pos <- by_pos[[i]][-rm.idx,] else up_by_pos <- by_pos[[i]]

          if(length(unique(up_by_pos$Alt)) == 0) { # If after applied filter all tags are gone
            rm <- rm + 1
            next()

          } else if (length(unique(up_by_pos$Alt)) > 1 ){ # If after applied filter the SNP remains multiallelic
            by_alt <- split.data.frame(up_by_pos, up_by_pos$Alt)
            by_alt_counts <- lapply(by_alt, function(x) apply(x[,-c(1:7)], 2, sum))
            total_counts <- sapply(by_alt_counts, sum)
            by_alt_counts <- by_alt_counts[order(total_counts, decreasing = T)] # Most frequency base will come first
            total_alt <- 0
            for(j in 1:length(by_alt_counts)){
              total_alt <- total_alt + by_alt_counts[[j]]
              if(j == 1) alt <- by_alt_counts[[j]] else  alt <- paste0(alt, ",",by_alt_counts[[j]])
            }
            info <- up_by_pos[,2:5]
            info$Alt <- paste0(names(by_alt_counts), collapse = ",")
            info <- unique.data.frame(info)

          } else { # If after applied filter, only one alternative remains
            alt <- apply(up_by_pos[,-c(1:7)], 2, sum)
            total_alt <- alt
            info <- unique.data.frame(up_by_pos[,2:5])
          }

        } else { # If rm_multiallelic_SNP set to FALSE and threshold is 0, keep all multiallelics
          by_alt <- split.data.frame(by_pos[[i]], by_pos[[i]]$Alt)
          by_alt_counts <- lapply(by_alt, function(x) apply(x[,-c(1:7)], 2, sum))
          total_counts <- sapply(by_alt_counts, sum)
          by_alt_counts <- by_alt_counts[order(total_counts, decreasing = T)] # Most frequency base will come first
          total_alt <- 0
          for(j in 1:length(by_alt_counts)){
            total_alt <- total_alt + by_alt_counts[[j]]
            if(j == 1) alt <- by_alt_counts[[j]] else  alt <- paste0(alt, ",",by_alt_counts[[j]])
          }

          info <- by_pos[[i]][,2:5]
          info$Alt <- paste0(names(by_alt_counts), collapse = ",")
          info <- unique.data.frame(info)
        }
      } else { # If SNP is not multiallelic
        up_by_pos <- by_pos[[i]]
        alt <- apply(up_by_pos[,-c(1:7)], 2, sum)
        total_alt <- alt
        info <- unique.data.frame(by_pos[[i]][,2:5])
      }

      ref <- apply(cloneID_unit[-which(cloneID_unit$AlleleID %in% alleles | duplicated(cloneID_unit$AlleleID)),-c(1:7)], 2, sum)
      tab_counts <- paste0(ref + total_alt, ":", ref, ":", ref, ",", alt)

      info <- c(info$Chromosome,
                info$SNP_position_in_Genome,
                paste0(info$Chromosome, "_", info$SNP_position_in_Genome),
                info$Ref, info$Alt, ".", ".", paste0("DP=", sum(c(ref, total_alt)),";",
                                                     "ADS=",sum(ref),",",sum(total_alt)), "DP:RA:AD")

      vcf_off_tag <- c(info, tab_counts)

      vcf_tag <- rbind(vcf_tag, vcf_off_tag)
    }
  }

  return(list(vcf_tag, rm))
}
