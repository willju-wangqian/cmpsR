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
#' data("bullets")
#' land2_3 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
#' land1_2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]
#' x <- land2_3$sig
#' y <- land1_2$sig
#' 
#' segments <- get_segs(x, len = 50)
#' nseg <- length(segments$segs)
#' seg_scale_max <- 3
#' npeaks.set <- c(5,3,1)
#' 
#' ccp.list <- lapply(1:nseg, function(nseg) {
#'  ccr.list <- lapply(1:seg_scale_max, function(seg_scale) {
#'    get_ccr_peaks(y, segments, seg_scale = seg_scale, nseg = nseg, npeaks = npeaks.set[seg_scale])
#'  })
#'  
#'  get_ccp(ccr.list, Tx = 25)
#' })
#' cmps <- get_CMPS(ccp.list, Tx = 25)
get_CMPS <- function(input.ccp, Tx = 25, order = T) {
  
  assert_that(is.numeric(Tx), is.logical(order))

  tt <- unlist(input.ccp)
  tt <- sort(tt)
  
  if(is.null(tt)) { 
    
    return(list(CMPS.score = 0, nseg = length(input.ccp), 
                congruent.pos = NA, congruent.seg = NA, 
                congruent.seg.idx = NA, pos.df = NA ))
    
    # return(list(CMPS.score = 0, rec.position = NULL, pos.df = NULL)) 
  }
  
  # pos.df <- data.frame(position = seq(min(tt), max(tt)))
  pos.df <- data.frame(position = tt)
  
  current.max <- 0
  # current.voter <- c()
  current.voter.list <- list()
  
  pos.df$cmps <- sapply(1:nrow(pos.df), function(i) {
    ccp.count <- sapply(1:length(input.ccp), function(j) {
      # if any one of the peaks belongs to the interval 
      # centered at this position, count 1
      any(abs(pos.df$position[i] - input.ccp[[j]]) <= Tx)
    })
    
    tmp.sum <- sum(ccp.count)
    
    if(tmp.sum >= current.max) {
      current.max <<- tmp.sum
      # current.voter <<- ccp.count
    }
    
    current.voter.list[[i]] <<- ccp.count
    
    return(tmp.sum)
    
  })
  
  # obtain the CMPS score
  ood <- order(pos.df$cmps, decreasing = T)
  ordered.pos.df <- pos.df[ood, ]
  ordered.current.voter <- current.voter.list[ood]
  
  # CMPS <- pos.df[order(pos.df$cmps, decreasing = T)[1], ]$cmps  
  CMPS <- ordered.pos.df$cmps[1]
  if(CMPS != current.max) {
    stop("unexpected: current.max didn't find the max CMPS score")
  }
  
  if(order) {
    pos.df <- ordered.pos.df
  }
  
  # recommended position for counting the CMPS score
  # rec.position <- floor(median(pos.df[pos.df$cmps == CMPS, ]$position))
  
  # congruent.pos <- round(median(ordered.pos.df$position[ordered.pos.df$cmps == CMPS]))
  congruent.pos.idx <- ceiling(sum(ordered.pos.df$cmps == CMPS)/2)
  
  if(congruent.pos.idx <= 0) {
    congruent.pos <- NA
    current.voter <- rep(FALSE, length(input.ccp))
  } else {
    congruent.pos <- ordered.pos.df$position[congruent.pos.idx]
    current.voter <- ordered.current.voter[[congruent.pos.idx]]
  }
  
  
  
  # browser()
  
  
  
  return(list(CMPS.score = CMPS, nseg = length(input.ccp), 
              congruent.pos = congruent.pos, congruent.seg = current.voter, 
              congruent.seg.idx = (1:length(input.ccp))[current.voter], pos.df = pos.df ))
}

#' Computes the CMPS score of a comparison between two bullet profiles/signatures
#' 
#' Compute the Congruent Matching Profile Segments (CMPS) score based on two bullet profiles/signatures.
#' The reference profile will be divided into consecutive, non-overlapping, basis segments of the same length.
#' Then the number of segments that are congruent matching will be found as the CMPS score. 
#' By default, `extract_feature_cmps` implements the algorithm with multi-peak insepction at three 
#' different segment scales. By setting `npeaks.set` as a single-length vector, users can switch to the algorithm 
#' with multi-peak inspection at the basis scale only.
#'
#' @param x a numeric vector, vector of the reference bullet signature/profile that will be divided into basis segments
#' @param y a numeric vector, vector of the comparison bullet signature/profile
#' @param seg_length integer, the length of a basis segment
#' @param Tx integer, the tolerance zone is `+/- Tx`
#' @param npeaks.set a numeric vector, specify the number of peaks to be found for each different scale of segment. The 
#' number of peaks for the highest segment scale must be 1, i.e., the very last element of `npeaks.set` must be 1.
#' * If `length(npeaks.set) == 1`, the algorithm uses multi-peak inspection only at the basis scale;
#' * If `length(npeaks.set) > 1`, the algorithm uses multi-peak inspection at 
#'    different segment scales. 
#' * By default, `npeaks.set = c(5,3,1)`. Including more segment scales will reduce the number of false positive results
#' @param include boolean, whether or not to return the full CMPS result (todo: update this; full_result 
#' has been discarded)
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
#' cmps_with_multi_scale <- extract_feature_cmps(land2_3$sig, land1_2$sig, include = "full_result" )
#' 
#' # algorithm with multi-peak inspection at the basis scale only
#' cmps_without_multi_scale <- extract_feature_cmps(land2_3$sig, land1_2$sig, 
#'                                                  npeaks.set = 5, include = "full_result" )
#' \dontrun{
#' library(tidyverse)
#' library(bulletxtrctr)
#' 
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
#'     extract_feature_cmps(a$lands$sig1, a$lands$sig2, include = "nseg")
#'   }))
#' 
#' # comparisons.cmps <- comparisons %>% 
#' #   mutate(cmps = aligned %>% purrr::map_dbl(.f = function(a) {
#' #     extract_feature_cmps(a$lands$sig1, a$lands$sig2, include = NULL)
#' #   }))
#' # comparisons.cmps %>% select(land1, land2, cmps) 
#' 
#' comparisons <- comparisons %>% 
#'   mutate(
#'     cmps_score = sapply(comparisons$cmps, function(x) x$CMPS.score),
#'     cmps_nseg = sapply(comparisons$cmps, function(x) x$nseg)
#'   )
#' 
#' cp1 <- comparisons %>% select(land1, land2, cmps_score, cmps_nseg)
#' cp1
#' }
#' @references 
#' Chen, Zhe, Wei Chu, Johannes A Soons, Robert M Thompson, John Song, 
#' and Xuezeng Zhao. 2019. “Fired Bullet Signature Correlation Using the 
#' Congruent Matching Profile Segments (CMPS) Method.” Forensic Science 
#' International, December, #109964. https://doi.org/10.1016/j.forsciint.2019.109964.

# extract_feature_cmps <- function(x, y, seg_length = 50, seg_scale_max = 3, Tx = 25, npeaks.set = c(5, 3, 1),
#                                  full_result = FALSE) {
extract_feature_cmps <- function(x, y, seg_length = 50, Tx = 25, npeaks.set = c(5, 3, 1),
                                 include = NULL) {
  assert_that(
    is.numeric(x), is.numeric(y), is.numeric(seg_length), # is.numeric(seg_scale_max),
    is.numeric(Tx), is.numeric(npeaks.set),
    (is.null(include) | is.character(include))
  )
  
  seg_scale_max <- length(npeaks.set)

  if (length(npeaks.set) != seg_scale_max) { 
    stop("seg_scale_max doesn't match the length of npeaks.set.")
  }
  
  # if (seg_scale_max > 1 & npeaks.set[length(npeaks.set)] != 1) {
  if (seg_scale_max > 1 & npeaks.set[seg_scale_max] != 1) {  
    stop("the length of the highest level must be 1, i.e., the last element of npeaks.set must be 1")
  }
  
  segments <- get_segs(x, seg_length)
  nseg <- length(segments$segs)
  
  # tic("Total")
  
  if (seg_scale_max == 1) {
    ccp.list <- lapply(1:nseg, function(nseg) {
      ccr <- get_ccr_peaks(y, segments, seg_scale = seg_scale_max, 
                           nseg = nseg, npeaks = npeaks.set[seg_scale_max])
      ccr$peaks.pos
      
    })
  } else if(seg_scale_max > 1) {
    
    
    ccp.list <- lapply(1:nseg, function(nseg) {
      
      # tic("time for one segment")
      ccr.list <- lapply(1:seg_scale_max, function(seg_scale) {
        
        # tic("time for get_ccr_peaks for one scale")
        get_ccr_peaks(y, segments, seg_scale = seg_scale, nseg = nseg, npeaks = npeaks.set[seg_scale])
        # toc()
      })
      
      rr <- get_ccp(ccr.list, Tx = Tx)
      
      # if(nseg == 2) {browser()}
      
      # toc()
      rr
      
    })
  } else {
    stop("seg_scale_max is invalid. Please use a positive integer instead.")
  }
  
  # browser()
  
  # tic("time for get_cmps")
  cmps <- get_CMPS(ccp.list, Tx = Tx)
  # toc()
  # cmps$nseg <- nseg
  if(cmps$nseg != nseg) {
    stop("unexpected: number of obs of ccp.list is not equal to nseg")
  }
  
  parameters <- list(x=x, y=y, seg_length=seg_length, Tx=Tx,
                     npeaks.set=npeaks.set, include=include)
  
  cmps$ccp.list <- ccp.list
  cmps$segments <- segments
  cmps$parameters <- parameters
  
  
  # toc()
  if(is.null(include)) {
    return(cmps$CMPS.score)
  } else {
    mm <- match.arg(include, c("nseg", "congruent.pos", "congruent.seg", "congruent.seg.idx", 
                               "pos.df", "ccp.list", "segments", "parameters", "full_result"),
                    several.ok = TRUE)
    if("full_result" %in% mm) {
      return(cmps)
    }
    
    mm <- c("CMPS.score", mm)
    return(cmps[mm])
  }
  
  # if(full_result) { return(cmps) } 
  # else { return(cmps$CMPS.score) }
}