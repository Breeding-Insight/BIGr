
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
##' @examples
##' # First generating example MADC files
##' temp_dir <- tempdir()
##' file1_path <- file.path(temp_dir, "madc1.csv")
##' file2_path <- file.path(temp_dir, "madc2.csv")
##' out_path <- file.path(temp_dir, "merged_madc.csv")
##'
##' # Data for file 1: Has SampleA and SampleB
##' df1 <- data.frame(
##'   AlleleID = c("chr1.1_0001|Alt_0002", "chr1.1_0001|Ref_0001", "chr1.1_0001|AltMatch_0001"),
##'   CloneID = c("chr1.1_0001", "chr1.1_0001", "chr1.1_0001"),
##'   AlleleSequence = c("GGG", "AAA", "TTT"),
##'   SampleA = c(10, 8, 0),
##'   SampleB = c(5, 4, 9),
##'   stringsAsFactors = FALSE,
##'   check.names = FALSE
##' )
##' write.csv(df1, file1_path, row.names = FALSE, quote = FALSE)
##'
##' # Data for file 2: Has SampleA (duplicate name) and SampleC, different rows
##' df2 <- data.frame(
##'   AlleleID = c("chr1.1_0001|Alt_0002", "chr1.1_0001|Ref_0001", "chr1.1_0001|AltMatch_0001"),
##'   CloneID = c("chr1.1_0001", "chr1.1_0001", "chr1.1_0001"),
##'   AlleleSequence = c("GGG", "AAA", "TTT"),
##'   SampleA = c(11, 7, 20),
##'   SampleC = c(1, 2, 6),
##'   stringsAsFactors = FALSE,
##'   check.names = FALSE
##' )
##' write.csv(df2, file2_path, row.names = FALSE, quote = FALSE)
##'
##' # --- 2. Run the merge function ---
##' # Use default suffixes (.x, .y) for the duplicated "SampleA"
##' merge_MADCs(madc_list = list(file1_path, file2_path),
##'             out_madc = out_path)
##'
##' # --- 3. Display the results (optional) ---
##' # Check if output file exists and print its content
##' if (file.exists(out_path)) {
##'   cat("--- Content of merged MADC file: ---\n")
##'   merged_data <- read.csv(out_path, check.names = FALSE)
##'   print(merged_data)
##'   # Expected output shows:
##'   # - 4 rows (union of alleles)
##'   # - Columns: ID cols, SampleA.x, SampleB, SampleA.y, SampleC
##'   # - NA values generated during merge replaced with 0
##' }
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
