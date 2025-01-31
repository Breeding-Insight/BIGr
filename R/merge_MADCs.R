
##' Merge MADC files
##'
##' If duplicated samples exist in different files, a extension (e.g."_file:1") will
##' be added at the end of the sample name
##'
##' @param ... one or more MADC files path
##'
##' @import dplyr
##'
##' @export
merge_MADCs <- function(..., out_madc=NULL){

  if(is.null(out_madc)) stop("Define output file name")

  files <- lapply(list(...), read.csv)

  # Check if there are same samples in different files
  # if there are, _rep:1 and _rep:2 will be added at the end of the sample name

  col_names <- unlist(sapply(files, function(x) colnames(x)[-c(1:3)]))
  repeated <- col_names[which(duplicated(col_names))]

  if(length(repeated) >0){
    for(i in 1:length(files)){
      repeated1 <- colnames(files[[i]])[which(colnames(files[[i]]) %in% repeated)]
      colnames(files[[i]])[match(repeated1, colnames(files[[i]]))] <- paste0(repeated1,"_file:",i)
    }
  }
  for(i in 1:(length(files)-1)){
    if(i==1) new <- full_join(files[[i]], files[[i+1]], by = c("AlleleID", "CloneID", "AlleleSequence"), )
    else new <- full_join(new, files[[i+1]], by = c("AlleleID", "CloneID", "AlleleSequence"))
  }

  new_zeros <- new %>% mutate_if(is.numeric,coalesce,0)

  write.csv(new_zeros, out_madc, row.names = FALSE, quote = FALSE)
}
