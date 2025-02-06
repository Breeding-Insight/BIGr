
##' Merge MADC files
##'
##' If duplicated samples exist in different files, a suffix will
##' be added at the end of the sample name. If run_ids is defined,
##' they are used as suffix, if not, files will be identified from
##' 1 to number of files, considering the order that was defined
##' in the function.
##'
##' @param ... one or more MADC files path
##' @param madc_list list containing path to MADC files to be merged
##' @param out_madc output merged MADC file path
##' @param run_ids vector of character defining the run ID for each file.
##' This ID will be added as a suffix in repeated sample ID in case they
##' exist in different files.
##'
##' @import dplyr
##'
##' @export
merge_MADCs <- function(..., madc_list=NULL, out_madc=NULL, run_ids=NULL){

  if(is.null(out_madc)) stop("Define output file name")

  if(!is.null(madc_list)) files <- lapply(madc_list, read.csv) else files <- lapply(list(...), read.csv)

  if(is.null(run_ids)) run_ids <- paste0("file",1:length(files)) else if(length(run_ids) != length(files))
    stop("run_ids vector should have some number of elements as the number of files listed as input.")

  # Check if there are same samples in different files
  # if there are, _file:1 and _file:2 will be added at the end of the sample name
  col_names <- unlist(sapply(files, function(x) colnames(x)[-c(1:3)]))
  repeated <- col_names[which(duplicated(col_names))]

  if(length(repeated) >0){
    for(i in 1:length(files)){
      repeated1 <- colnames(files[[i]])[which(colnames(files[[i]]) %in% repeated)]
      colnames(files[[i]])[match(repeated1, colnames(files[[i]]))] <- paste0(repeated1,"_", run_ids[i])
    }
  }
  for(i in 1:(length(files)-1)){
    if(i==1) new <- full_join(files[[i]], files[[i+1]], by = c("AlleleID", "CloneID", "AlleleSequence"), )
    else new <- full_join(new, files[[i+1]], by = c("AlleleID", "CloneID", "AlleleSequence"))
  }

  new_zeros <- new %>% mutate_if(is.numeric,coalesce,0)

  write.csv(new_zeros, out_madc, row.names = FALSE, quote = FALSE)
}
