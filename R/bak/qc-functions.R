#' QC Thresholding
#' 
#' Use \code{threshold_data} to combine  This implementation follows
#' the rules outlined in the \emph{How to Handle Missing Data} for the Fluidigm
#' provided \emph{SINGuLAR Analysis Toolset 3.0 User Guide}. "For expression data retrieved from a BioMark qPCR experiment, a missing
#' value is defined as data that does not have a Ct value of 999 and is a FAIL." Returns
#' \code{Log2Ex} values as defined in the tech note.
#' 
#' @param ct A matrix containing \eqn{C_{t}} values for the experiment
#' @param qc A matrix containing \code{PASS/FAIL} values for the experiment
#' @param LOD The limit of detection for the experiment. Defaults to 25.
#' 

threshold_data <- function(ct, qc, LOD = 25) {
  
  # confirm that all the rownames and column names match
  stopifnot(
    all(colnames(ct) == colnames(qc)), 
    all(rownames(ct) == rownames(qc))
  )
  
  # set up our masks
  lod_mask <- ct > LOD
  fail_mask <- qc == 'Fail'
  na_mask <- !lod_mask & fail_mask # anything that was detected but was marked as FAIL needs to be NA
  
  # set our new values
  ct[lod_mask] <- LOD
  ct[na_mask] <- NA
  
  return(LOD - ct)
  
}


#' Merge Technical Replicates
#' 
#' Use \code{merge_tech_reps} to merge technical replicates from a Fluidigm
#' experiment. Currently just takes the mean of the two columns that are run
#' in duplicate. 
#' 
#' @param ct A matrix of ct values
#' @param na_rm Should NA values be ignored when merging values?
#' 

merge_tech_reps <- function(ct, na_rm = T) {
  
  # create new matrix -- 1 col per gene
  all_genes <- unique(colnames(ct))
  new_dat <- matrix(NA, ncol = length(all_genes), nrow = nrow(ct))
  rownames(new_dat) <- rownames(ct)
  colnames(new_dat) <- all_genes
  
  # average across plate replicates
  for (i in seq_along(all_genes)) {
    sel <- which(colnames(ct) == all_genes[i])
    
    if (length(sel) > 1) {
      to_ave <- ct[,sel]
      meaners <- apply(to_ave, 1, mean, na.rm = na_rm)
    } else {
      meaners <- ct[,sel]
    }
    new_dat[ ,all_genes[i]] <- meaners
  }
  
  return(new_dat)
}



#' Make pass/fail QC plot from QC data frame
#' 
#' Use \code{pass_fail_qc} to perform graphical PASS/FAIL qualtiy assessment
#' for a set of Fluidigm data. \code{qc} should be a \code{data.frame} containing 
#' the PASS/FAIL assessment for each of the 
#'
#' @param qc A vector of sample name prefixes to exclude from the 
#' @param merge_genes Logical determining whether or not to merge technical 
#' replicate genes prior to assessment. If true, 
#' @param exclude_regex A vector of sample name prefixes to exclude from the 
#' returned data set -- e.x. \code{'^NTC*'} will exclude all rows with a sample
#'
#' @export

pass_fail_qc <- function(qc = NULL, merge_genes = F, exclude_regex = NULL) {
  
  
  if (is.null(qc)) stop('Please provide QC input.')
  
  
  # clean any samples to exclude..
  if (!is.null(exclude_regex)) {
    
    sapply()
    
    bs <- grep("^NT", qc[['Sample_ID']])
    
    if (length(bs)) qc <- qc[-bs, ]
    
    
  }
  
  
  qcm <- qc
  qcm <- qc[,-1]
  qcm <- as.matrix(qcm)
  qcm[qcm == 'Fail'] <- 0
  qcm[qcm == 'Pass'] <- 1
  qcm <- apply(qcm, 2, as.numeric)
  rownames(qcm) <- qc[,1]
  
  borgmisc::heat_misc(qcm, z_score = F, row_clust = T, col_clust = T, rang = c(0,1), 
                      axis_scale = .5)
  
  
  miss <- list()
  
  ## quality plots
  mm <- melt(qc, id.vars = 'Sample_ID')
  
  p1 <- ggplot(mm, aes(x = variable, y = Sample_ID, fill = value)) + 
    geom_tile() + scale_fill_brewer(palette = 'Set1') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  mm_sum <- mm %>%
    group_by(variable) %>%
    summarise(
      pct_fail = sum(value == 'Fail') / length(value),
      pct_pass = sum(value == 'Pass') / length(value)
    ) %>%
    ungroup() %>%
    arrange(pct_fail)
  
  names(mm_sum)[1] <- 'gene'
  mm_sum[['gene']] <- factor(mm_sum$gene, levels = mm_sum$gene)
  p2 <- mm_sum %>%
    melt(id.vars = 'gene') %>%
    ggplot(aes(x = gene, y = value, fill = variable)) + 
    geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_brewer(palette = 'Set1')
  
  miss[[1]] <- list()
  miss[[1]][[1]] <- mm_sum
  
  mm_sum <- mm %>%
    group_by(Sample_ID) %>%
    summarise(
      pct_fail = sum(value == 'Fail') / length(value),
      pct_pass = sum(value == 'Pass') / length(value)
    ) %>%
    ungroup() %>%
    arrange(pct_fail)
  
  names(mm_sum)[1] <- 'id'
  mm_sum[['id']] <- factor(mm_sum$id, levels = mm_sum$id)
  p3 <- mm_sum %>%
    melt(id.vars = 'id') %>%
    ggplot(aes(x = id, y = value, fill = variable)) + 
    geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_brewer(palette = 'Set1')
  
  miss[[1]][[2]] <- mm_sum
  
  borgmisc::multiplot(p1, p2, p3, cols = 2)
  
}

