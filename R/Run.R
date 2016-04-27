#' @title Run class
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @description A class representation of a Bio-Rad plate.
#'
#' @format \code{\link{Plate}} class generator
#'
#' @usage \code{run <- Run$new(xl_file = '~/Google Drive/spectrum/stats/tarnawski/yafa/data/Heat map results YaFa Explant 1.xlsx')}
#'
#' @keywords data
#'


###
# notes:
#     - can you steal some stuff from here: http://www.biomedcentral.com/content/pdf/1471-2164-13-296.pdf
###



Run <- R6::R6Class(
  "Run",
  public = list(

    raw_data = NA,
    quality_data = NA, 
    normalized_data = NA,
    normalization_method = NA, 
    metadata = NA,

    
    # public methods
    initialize = function(xl_file, sheet = 1) {

      if (!missing(xl_file)) {
        dat <- readxl::read_excel(xl_file, sheet)
        
        # find where data is located
        ct_idx <- which(dat[,1] == 'FAM-MGB Ct')
        qual_idx <- which(dat[,1] == 'Quality Results')
        
        # parse out the ct_values
        g_start <- ct_idx + 2
        gene_names <- borgmisc::trim(setNames(unlist(dat[g_start, -c(1,2)]), NULL))
        sample_names <- borgmisc::trim(setNames(unlist(dat[(g_start + 1):(qual_idx - 1), 2]), NULL))
        ct_values <- as.matrix(dat[(g_start + 1):(qual_idx - 1), 3:ncol(dat)])
        ct_values <- apply(ct_values, 2, as.numeric)
        colnames(ct_values) <- gene_names
        rownames(ct_values) <- sample_names
        
        # parse out the quality data      
        g_start <- qual_idx + 2
        s_end <- ((g_start + 1):nrow(dat))[min(which(is.na(dat[(g_start + 1):nrow(dat), 2])))]
        gene_names <- borgmisc::trim(setNames(unlist(dat[g_start, -c(1,2)]), NULL))
        sample_names <- borgmisc::trim(setNames(unlist(dat[(g_start + 1):(s_end - 1), 2]), NULL))
        q_values <- as.matrix(dat[(g_start + 1):(s_end - 1), 3:ncol(dat)])
        colnames(q_values) <- gene_names
        rownames(q_values) <- sample_names

        # sanity check
        all_clear <- all(
          rownames(ct_values) == rownames(ct_values), 
          colnames(ct_values) == colnames(ct_values)
          )
        
        if (all_clear) {
          self$raw_data <- ct_values
          self$quality_data <- q_values
        } else {
          stop('There is a problem with the input file... Exiting...')
        }
      }
      
    },
    
    # getting around known issue: https://github.com/wch/R6/issues/51
    say_hi = function(x) {
      print('hello')
    }
  )
)



Run$set("public", "threshold_data", function(LOD = 25) {

  # set up our masks
  lod_mask <- self$raw_data > LOD
  fail_mask <- self$quality_data == 'Fail'
  na_mask <- !lod_mask & fail_mask # anything that was detected but was marked as FAIL needs to be NA
  
  ct <- self$raw_data
  ct[lod_mask] <- LOD
  ct[na_mask] <- NA
  
  self$normalized_data <- LOD - ct

})



Run$set("public", "merge_tech_reps", function(na_rm = T) {
  
  # create new matrix -- 1 col per gene
  all_genes <- unique(colnames(self$normalized_data))
  new_dat <- matrix(NA, ncol = length(all_genes), nrow = nrow(self$normalized_data))
  rownames(new_dat) <- rownames(self$normalized_data)
  colnames(new_dat) <- all_genes
  
  # average across plate replicates
  for (i in seq_along(all_genes)) {
    sel <- which(colnames(self$normalized_data) == all_genes[i])
    
    if (length(sel) > 1) {
      to_ave <- self$normalized_data[,sel]
      meaners <- apply(to_ave, 1, mean, na.rm = na_rm)
    } else {
      meaners <- self$normalized_data[,sel]
    }
    new_dat[ ,all_genes[i]] <- meaners
  }
  
  self$normalized_data <- new_dat
  
})


