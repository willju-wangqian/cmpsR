#' Compute the CMPS score
#'
#' Compute the CMPS score from a list of positions of (consistent) correlation peaks.
#' @param input.ccp a list of positions for (consistent) correlation peaks
#' @param Tx integer, the tolerance zone is `+/- Tx`
#' @param order boolean, whether or not to order peak positions based on their CMPS score
#'
#' @return a list of three components:
#' * `CMPS.score`: computed CMPS score
#' * `rec.position`: the recommended position that results in the CMPS score
#' * `pos.df`: a dataframe that includes all positions and their corresponding CMPS score 
#' @export
#' @importFrom assertthat assert_that
#'
#' @examples
get_CMPS <- function(input.ccp, Tx = 25, order = T) {
  
  assert_that(is.numeric(Tx), is.logical(order))

  tt <- unlist(input.ccp)
  
  if(is.null(tt)) { 
    return(list(CMPS.score = 0, rec.position = NULL, pos.df = NULL)) 
  }
  
  pos.df <- data.frame(position = seq(min(tt), max(tt)))
  
  pos.df$cmps <- sapply(1:nrow(pos.df), function(i) {
    ccp.count <- sapply(1:length(input.ccp), function(j) {
      # if any one of the peaks belongs to the interval 
      # centered at this position, count 1
      any(abs(pos.df$position[i] - input.ccp[[j]]) <= Tx)
    })
    
    sum(ccp.count)
  })
  
  # obtain the CMPS score
  CMPS <- pos.df[order(pos.df$cmps, decreasing = T)[1], ]$cmps  
  
  if(order) {
    pos.df <- pos.df[order(pos.df$cmps, decreasing = T), ]
  }
  
  # recommended position for counting the CMPS score
  rec.position <- floor(median(pos.df[pos.df$cmps == CMPS, ]$position))
  
  return(list(CMPS.score = CMPS, rec.position = rec.position, pos.df = pos.df))
}

#' Computes the CMPS score of a comparison between two bullet profiles/signatures
#' 
#' Compute the Congruent Matching Profile Segments (CMPS) score based on two bullet profiles/signatures.
#' The reference profile will be divided into consecutive, non-overlapping, basis segments of the same length.
#' Then the number of segments that are congruent matching will be found as the CMPS score. 
#' By default, `extract_feature_cmps` implements the algorithm with multi-peak insepction at three 
#' different segment scales. By setting `seg_scale_max = 1`, users can switch to the algorithm 
#' with multi-peak inspection at the basis scale only.
#'
#' @param x a numeric vector, vector of the reference bullet signature/profile that will be divided into basis segments
#' @param y a numeric vector, vector of the comparison bullet signature/profile
#' @param seg_length integer, integer, the length of a basis segment
#' @param seg_scale_max integer, the maximum scale of the segment.
#' * If setting `seg_scale_max = 1`, the algorithm uses multi-peak inspection only at the basis scale;
#' * If setting `seg_scale_max` to be an integer greater than 1, the algorithm uses multi-peak inspection at 
#'    different segment scales. 
#' * By default, `seg_scale_max = 3`. Increasing `seg_scale_max` will reduce the number of false positive results
#' @param Tx integer, the tolerance zone is `+/- Tx`
#' @param npeaks.set a numeric vector, specify the number of peaks to be found for each different scale of segment. The 
#' number of peaks for the highest segment scale must be 1, i.e., the very last element of `npeaks.set` must be 1.
#' @param full_result boolean, whether or not to return the full CMPS result
#'
#' @return integer or a list
#' * if `full_result = FALSE`, return the CMPS number only
#' * if `full_result = TRUE`, return a list of four elements:
#'     + `CMPS.score`: computed CMPS score
#'     + `rec.position`: the recommended position that results in the CMPS score
#'     + `pos.df`: a dataframe that includes all positions and their corresponding CMPS score 
#'     + `nseg`: the number of basis segments obtained from the reference profile
#' 
#' @export
#' @importFrom assertthat assert_that
#'
#' @examples
#' library(bulletxtrctr)
#' library(tidyverse)
#' library(x3ptools)
#' library(CMPS)
#' 
#' data("bullets")
#' land2_3 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
#' land1_2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]
#' 
#' # compute cmps
#' 
#' # algorithm with multi-peak insepction at three different segment scales
#' cmps_with_multi_scale <- extract_feature_cmps(land2_3$sig, land1_2$sig, full_result = TRUE)
#' 
#' # algorithm with multi-peak inspection at the basis scale only
#' cmps_without_multi_scale <- extract_feature_cmps(land2_3$sig, land1_2$sig, seg_scale_max = 1, 
#'                                                  npeaks.set = 5, full_result = TRUE)
#' 
#' # full comparison of the two bullets
#' lands <- unique(bullets$bulletland)
#' 
#' comparisons <- data.frame(expand.grid(land1 = lands[1:6], land2 = lands[7:12]), 
#'                           stringsAsFactors = FALSE)
#' 
#' comparisons <- comparisons %>% mutate(
#'   aligned = purrr::map2(.x = land1, .y = land2, 
#'                         .f = function(xx, yy) {
#'                           land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
#'                           land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
#'                           land1$bullet <- "first-land"
#'                           land2$bullet <- "second-land"
#'                           
#'                           sig_align(land1$sig, land2$sig)
#'                         }))
#' 
#' comparisons <- comparisons %>% 
#'   mutate(cmps = aligned %>% purrr::map(.f = function(a) {
#'     extract_feature_cmps(a$lands$sig1, a$lands$sig2, full_result = TRUE)
#'   }))
#' 
#' comparisons <- comparisons %>% 
#'   mutate(
#'     cmps_score = sapply(comparisons$cmps, function(x) x$CMPS.score),
#'     cmps_nseg = sapply(comparisons$cmps, function(x) x$nseg)
#'   )
#' 
#' cp1 <- comparisons %>% select(land1, land2, cmps_score, cmps_nseg)
#' cp1
extract_feature_cmps <- function(x, y, seg_length = 50, seg_scale_max = 3, Tx = 25, npeaks.set = c(5, 3, 1),
                                 full_result = FALSE) {
  
  assert_that(
    is.numeric(x), is.numeric(y), is.numeric(seg_length),
    is.numeric(seg_scale_max), is.numeric(Tx), is.numeric(npeaks.set),
    is.logical(full_result)
  )

  if (length(npeaks.set) != seg_scale_max) { 
    stop("seg_scale_max doesn't match the length of npeaks.set.")
  }
  
  if (seg_scale_max > 1 & npeaks.set[length(npeaks.set)] != 1) {
    stop("the length of the highest level must be 1, i.e., the last element of npeaks.set must be 1")
  }
  
  segments <- get_segs(x, seg_length)
  nseg <- length(segments$segs)
  if (seg_scale_max == 1) {
    ccp.list <- lapply(1:nseg, function(nseg) {
      ccr <- get_ccr_peaks(y, segments, seg_scale = seg_scale_max, 
                           nseg = nseg, npeaks = npeaks.set[seg_scale_max])
      ccr$peaks.pos
      
    })
  } else if(seg_scale_max > 1) {
    ccp.list <- lapply(1:nseg, function(nseg) {
      ccr.list <- lapply(1:seg_scale_max, function(seg_scale) {
        get_ccr_peaks(y, segments, seg_scale = seg_scale, nseg = nseg, npeaks = npeaks.set[seg_scale])
      })
      
      get_ccp(ccr.list, Tx = Tx)
    })
  } else {
    stop("seg_scale_max is invalid. Please use a positive integer instead.")
  }
  
  cmps <- get_CMPS(ccp.list, Tx = Tx)
  cmps$nseg <- nseg
  
  if(full_result) { return(cmps) } 
  else { return(cmps$CMPS.score) }
}