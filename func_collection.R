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

get_seg_level <- function(segments, nseg, level = 1){
  # obtain a longer segment, centered at the current segment
  # cannot exceed 1 or the max index of x
  # segments - the collection of all segments, obtained by get_segs()
  # nseg - the index of the segment to be augmented
  # level - the augment level
  
  x <- segments$x
  idx <- segments$index[[nseg]]
  
  ct <- floor(median(idx))
  unitt <- max(idx) - ct
  ############# ???
  # cut at the max or min
  min.idx <- max(ct - unitt * level, 1)
  max.idx <- min(ct + unitt * level, length(x))
  
  return(list(aug_seg=x[min.idx:max.idx], aug_idx=min.idx:max.idx))
}

