#' Read a set of Fluidigm ssqPCR data
#' 
#' Read a TSV file derived from Fluidigm output. Return as matrix formatted as
#' we want.

f2mat <- function(file_name) {
  dat <- read.delim(file_name, check.names = F)
  dat_m <- dat
  dat_m <- dat_m[,-1]
  dat_m <- as.matrix(dat_m)
  rownames(dat_m) <- dat[,1]
  colnames(dat_m) <- colnames(dat)[2:ncol(dat)]
  return(dat_m)
}



#' Read a set of Fluidigm ssqPCR data
#'
#' depreciated -- please use the R6 class interface
#'
#' @param data_dir Path to your data files
#' @param ct_file  TSV file containing raw \eqn{C_{q}} values
#' @param qc_file  TSV file containing pass/fail calls from the machine. 
#' @param exclude_regex A vector of sample name prefixes to exclude from the 
#' returned data set -- e.x. \code{'^NTC*'} will exclude all rows with a sample
#' name starting with "NTC".
#'
#' @export


read_experiment <- function(xl_file = "", sheet = 1, exclude_regex = NULL) {

    
NULL  
  
  
  
  
}
