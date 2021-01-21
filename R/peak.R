#' Identify peaks of a cross correlation curve
#' 
#' Given a comparison profile and a segment, `get_ccr_peaks` computes the 
#' cross correlation curve and finds peaks of the curve.
#' @param comp a nueric vector, vector of the bullet comparison profile
#' @param segments list with basis segments and their corresponding indices in the original profile, obtianed by `get_segs()`
#' @param seg_scale integer. the scale of the segment.
#' * `seg_scale` = 1: the length remains the same;
#' * `seg_scale` = 2: the length will be doubled. 
#' * `seg_scale` = 3: the length will be doubled twice. 
#' @param nseg integer. `nseg` = 3: the third segment in `segments`
#' @param npeaks integer. the number of peaks to be identified.
#'
#' @return a list consisting of:
#' * `ccr`: the cross correlation curve
#' * `adj.pos`: indices of the curve
#' * `peaks.pos`: position of the identified peaks
#' * `peaks.heights`: the cross correlation value (height of the curve) of the peaks
#' @export
#' 
#' @importFrom zoo rollapply
#' @importFrom assertthat assert_that
#'
#' @examples
#' data("bullets")
#' land2_3 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
#' land1_2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]
#' x <- land2_3$sig
#' y <- land1_2$sig
#' 
#' segments <- get_segs(x, len = 50)
#' 
#' # compute ccf based on y and segment 7 with scale 1, then identify 5 highest peaks
#' ccrpeaks <- get_ccr_peaks(y, segments = segments, nseg = 7, npeaks = 5, seg_scale = 1)
get_ccr_peaks <- function(comp, segments, seg_scale, nseg = 1, npeaks = 5){

  assert_that(
    is.numeric(comp), is.numeric(seg_scale), is.numeric(nseg), is.numeric(npeaks)
  )
  
  # tic("time before finding the peaks")
  if(seg_scale == 1) {
    # compute for the basis segment
    # tic("time for get_ccf4 1")
    seg <- segments$segs[[nseg]]
    
    min.overlap <- min(length(seg[!is.na(seg)])*0.9, round(length(comp)*0.1))
    
    ccr <- get_ccf4(comp, seg, min.overlap = min.overlap)
    # toc()
    tmp_pos <- segments$index[[nseg]][1]
  } else{
    # find the increased segment, then compute
    # tic("time for get scale")
    tt <- get_seg_scale(segments, nseg, scale = seg_scale)
    # toc()
    
    # if(length(tt$aug_seg) == 1676) {browser()}
    
    min.overlap <- min(length(tt$aug_seg[!is.na(tt$aug_seg)])*0.9, round(length(comp)*0.1))
    
    # tic("time for get_ccf4 2")
    ccr <- get_ccf4(comp, tt$aug_seg, min.overlap = min.overlap)
    # toc()
    tmp_pos <- tt$aug_idx[1]
  }
  # toc()
  
  # tic("find peaks")
  # get all peaks and peak.heights
  find_maxs <- rollapply(ccr$ccf, 3, function(x) max(x) == x[2], 
                         fill = list(NA, NA, NA))
  peaks <- which(find_maxs)
  
  # new stuff
  od <- order(ccr$ccf[peaks], decreasing = TRUE)[1:npeaks]
  # adjust the position
  adj_pos <- ccr$lag - tmp_pos + 1
  
  peaks.heights <- ccr$ccf[peaks][od]
  peaks.pos <- adj_pos[peaks][od]
  # toc()
  
  # tic("sort")
  # get the position of peaks
  # peak_pos <- sort(peaks[order(peaks.heights, decreasing = T)][1:npeaks] - tmp_pos)
  # peak_height <- ccr$ccf[peak_pos + tmp_pos]
  # toc()
  
  return(list(ccr = ccr, adj.pos = adj_pos,
              peaks.pos = peaks.pos, peaks.heights = peaks.heights))
  
  
  # return(list(ccr = ccr, adj.pos = adj_pos,
  #             peaks.pos = peak_pos, peaks.heights = peak_height))
}

#' Identify at most one consistent correlation peak (ccp) 
#'
#' If multi segment lengths strategy is being used, at most one consistent correlation
#' peak (ccp) will be found for the corresponding basis segment. If the ccp cannot be identified, 
#' return `NULL`
#' @param ccr.list list, obtained by `get_ccr_peaks`
#' @param Tx integer, the tolerance zone is `+/- Tx`
#'
#' @return integer, the position of the ccp if it is identified; `NULL` otherwise. 
#' @export
#' @importFrom assertthat assert_that
#'
#' @examples
#' data("bullets")
#' land2_3 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
#' land1_2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]
#' x <- land2_3$sig
#' y <- land1_2$sig
#' 
#' segments <- get_segs(x, len = 50)
#' 
#' # identify the consistent correlation peak when ccf curves are computed
#' # based on y and segment 7 in 3 different scales;
#' # the number of peaks identified in each scale are 5, 3, and 1, respectively.
#' seg_scale_max <- 3
#' npeaks.set <- c(5,3,1)
#' 
#' ccr.list <- lapply(1:seg_scale_max, function(seg_scale) {
#'   get_ccr_peaks(y, segments, seg_scale = seg_scale, nseg = 7, npeaks = npeaks.set[seg_scale])
#' })
#' 
#' get_ccp(ccr.list, Tx = 25)
get_ccp <- function(ccr.list, Tx = 25){
  
  assert_that(is.numeric(Tx))
  
  # the number of different scales we have
  seg_level <- length(ccr.list)
  
  # the highest level has only one position, 
  # set it as a basis
  basis <- ccr.list[[seg_level]]$peaks.pos
  
  # if the entire basis segment is NA 
  if(length(basis) == 0 & all(is.na(ccr.list[[seg_level]]$ccr$ccf))) { return(NULL) }
  
  # for the purpose of debugging
  # if(length(basis) != 1) {stop("the length of the highest level must be 1, i.e., 
  #                              the last element of npeaks.set must be 1")}
  
  ccp <- lapply(1:(seg_level-1), function(level) {
    ccr.list[[level]]$peaks.pos[abs(ccr.list[[level]]$peaks.pos - basis) <= Tx]
  })
  ccp <- unlist(c(ccp, basis))
  if(length(ccp) >= seg_level) {return(basis)} else {return(NULL)}
}

