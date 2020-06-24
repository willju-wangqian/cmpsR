comparisons <- data.frame(expand.grid(land1 = lands[1:6], land2 = lands[7:12]), stringsAsFactors = FALSE)

comparisons <- comparisons %>% mutate(aligned = purrr::map2(.x = land1, .y = land2, 
                                                                      .f = function(xx, yy) {
                                                                        land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
                                                                        land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
                                                                        land1$bullet <- "first-land"
                                                                        land2$bullet <- "second-land"
                                                                        
                                                                        sig_align(land1$sig, land2$sig)
                                                                      }))
aligned <- comparisons$aligned[[15]]

x <- aligned$lands$sig1
y <- aligned$lands$sig2

get_segs2 <- function(x, len = 50){
  # divide a signature into segments
  # x - the signature to be divided
  # n - the desired nubemr of segment 
  ############################################################################
  
  idx <- 1:length(x)
  idx[is.na(x)] <- NA
  
  xx <- na.trim(x)
  idx <- na.trim(idx)
  segs <- split(xx, ceiling(seq_along(xx)/len)) 
  index <- split(idx, ceiling(seq_along(xx)/len))
  
  check_length <- sapply(segs, length) == len
  
  return(list(segs = segs[check_length], index = index[check_length], x = x))
}
 
segments <- get_segs2(y, len = 50)
# segments <- get_segs(x, 25)

# check_length <- sapply(segments$segs, length) != 45

str(segments$segs)
str(segments$index)


xx <- comparisons$aligned[[15]]$lands$sig2
iidx <- 1:length(xx)
iidx[is.na(xx)] <- NA

tt <- na.trim(xx)
iiddx <- na.trim(iidx)

which(xx == tt)


summary(na.trim(comparisons$aligned[[15]]$lands$sig1))
length(na.trim(comparisons$aligned[[15]]$lands$sig1))

extract_feature_cmps <- function(x, y, seg_length = 50, seg_scale_max = 3, Tx = 25, npeaks.set = c(5, 3, 1),
                                 full_result = FALSE) {
  # compute CMPS score
  
  # x - the reference signature/profile 
  # y - the comparison signature/profile
  # nseg - the number of basis segments used for the algorithm
  # seg_scale_max - the number of scales (different lengths) 
  # used for multi segment lengths
  # Tx - size of the tolerance zone
  # npeaks.set - vector that contains the number of peaks for each segment scale
  # full_result - whether or not to return the registration position used to 
  # find the CMPS
  #################################################################
  if (length(npeaks.set) != seg_scale_max) { 
    print("Need to specify the number of peaks for each segment scale.")
    return(NULL)
  }
  
  segments <- get_segs2(x, seg_length)
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
    print("seg_scale_max is invalid. Please use a positive integer instead.")
    return(NULL)
  }
  
  cmps <- get_CMPS(ccp.list, Tx = Tx)
  cmps$nseg <- nseg
  
  if(full_result) { return(cmps) } 
  else { return(cmps$CMPS.score) }
}

system.time({
  comparisons <- comparisons %>% 
    mutate(cmps = aligned %>% purrr::map(.f = function(a) {
      extract_feature_cmps(a$lands$sig1, a$lands$sig2, full_result = T)
    }))
})

comparisons <- comparisons %>% 
  mutate(
    cmps_score = sapply(comparisons$cmps, function(x) x$CMPS.score),
    cmps_nseg = sapply(comparisons$cmps, function(x) x$nseg)
  )

cp1 <- comparisons %>% select(land1, land2, cmps_score, cmps_nseg)

cp2 <- comparisons %>% select(land1, land2, cmps_score, cmps_nseg)
cp2

identical(cp1, cp2)


