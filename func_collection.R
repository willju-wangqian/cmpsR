get_ccf2 <- function(x, y, idx){
  # compute the cross correlation
  # x - a short segment
  # y - the comparison signature
  # idx - the index of x
  
  x <- as.vector(unlist(x))
  y <- as.vector(unlist(y))
  
  min.overlap = length(x[!is.na(x)])
  
  nx <- length(x)
  ny <- length(y)
  assert_that(is.numeric(x), is.numeric(y))
  assert_that(nx > 0, ny > 0, nx <= ny)
  xx <- c(rep(NA, ny - min.overlap), x, rep(NA, ny - min.overlap))
  yy <- c(y, rep(NA, length(xx) - ny))
  lag.max <- length(yy) - length(y)
  lags <- lag.max:0
  # lags <- 0:lag.max
  cors <- sapply(lags, function(lag) {
    cor(xx, lag(yy, lag), use = "pairwise.complete")
  })
  ns <- sapply(lags, function(lag) {
    dim(na.omit(cbind(xx, lag(yy, lag))))[1]
  })
  cors[ns < min.overlap] <- NA
  # lag <- lags - (ny - min.overlap)
  lag <- (ny - min.overlap) - lags - (idx[1] - 1)
  return(list(lag = lag, ccf = cors))
}

get_ccf3 <- function (x, y, min.overlap = round(0.1 * max(length(x), length(y)))) 
{
  # requires x to be the longer signature
  x <- as.vector(unlist(x))
  y <- as.vector(unlist(y))
  nx <- length(x)
  ny <- length(y)
  assert_that(is.numeric(x), is.numeric(y))
  assert_that(nx > 0, ny > 0, nx >= ny) # this is the only change
  xx <- c(rep(NA, ny - min.overlap), x, rep(NA, ny - min.overlap))
  yy <- c(y, rep(NA, length(xx) - ny))
  lag.max <- length(yy) - length(y)
  lags <- 0:lag.max
  cors <- sapply(lags, function(lag) {
    cor(xx, lag(yy, lag), use = "pairwise.complete")
  })
  ns <- sapply(lags, function(lag) {
    dim(na.omit(cbind(xx, lag(yy, lag))))[1]
  })
  cors[ns < min.overlap] <- NA
  lag <- lags - (ny - min.overlap)
  return(list(lag = lag, ccf = cors))
}

get_segs <- function(x, n){
  # divide a signature into segments
  # x - the signature to be divided
  # n - the desired nubemr of segment 
  
  segs <- split(x, ceiling(seq_along(x)/(length(x)/n))) 
  index <- split(1:length(x), ceiling(seq_along(x)/(length(x)/n)))
  return(list(segs = segs, index = index, x = x))
}

get_seg_scale <- function(segments, nseg, scale = 1){
  # obtain a longer segment, centered at the current segment
  # cannot exceed 1 or the max index of x
  # segments - the collection of all segments, obtained by get_segs()
  # nseg - the index of the segment to be augmented
  # scale - the augment scale
  
  x <- segments$x
  idx <- segments$index[[nseg]]
  
  ct <- floor(median(idx))
  unitt <- max(idx) - ct
  ############# ???
  # cut at the max or min
  min.idx <- max(ct - unitt * scale, 1)
  max.idx <- min(ct + unitt * scale, length(x))
  
  return(list(aug_seg=x[min.idx:max.idx], aug_idx=min.idx:max.idx))
}

get_ccr_peaks <- function(comp, segments, seg_scale, nseg = 1, npeaks = 5){
  # obtain the position of peaks of the cross correlation curve between 
  # the chosen segment and the comparison profile
  
  # comp - the comparison profile
  # segments - the collection of all basis segments of the reference profile;
  # generated from get_segs();
  # seg_scale - integer; the length or scale of a segment will be increased
  # to this number (seg_scale) times the length of the basis segment; 
  # nseg - which segment will be investigated;
  # npeaks - the number of peaks to be identified
  
  if(seg_scale == 1) {
    # compute for the basis segment
    ccr <- get_ccf3(comp, segments$segs[[nseg]], min.overlap = length(segments$segs[[nseg]]))
    tmp_pos <- segments$index[[nseg]][1]
  } else{
    # find the increased segment, then compute
    tt <- get_seg_scale(segments, nseg, scale = seg_scale)
    ccr <- get_ccf3(comp, tt$aug_seg, min.overlap = length(tt$aug_seg))
    tmp_pos <- tt$aug_idx[1]
  }
  
  # get all peaks
  rr <- sig_get_peaks(ccr$ccf, smoothfactor = 1, window = F, striae = F)
  
  # adjust the position
  adj_pos <- ccr$lag - tmp_pos + 1
  
  # get the position of peaks
  peak_pos <- sort(rr$peaks[order(rr$peaks.heights,decreasing = T)][1:npeaks] - tmp_pos)
  peak_height <- rr$dframe$smoothed[peak_pos + tmp_pos]
  return(list(ccr = ccr, ccrpeaks = rr, adj.pos = adj_pos, 
              peaks.pos = peak_pos, peaks.heights = peak_height))
}









