i <- 33

x <- comparisons.cmps$aligned[[i]]$lands$sig1
y <- comparisons.cmps$aligned[[i]]$lands$sig2

nseg <- 25

segments <- get_segs(x, nseg)

rbind(data.frame(value = x, sig = "sig1", x = 1:length(x)), data.frame(value = y, sig = "sig2", x = 1:length(y))) %>%
  ggplot(aes(x = x, y = value, color = sig)) + geom_line()

seg_scale_max <- 3
npeaks.set <- c(5, 3, 1)
Tx <- 25

nseg <- 2

ccr.list <- lapply(1:seg_scale_max, function(seg_scale) {
  get_ccr_peaks(y, segments, seg_scale = seg_scale, nseg = nseg, npeaks = npeaks.set[seg_scale])
})

get_ccp(ccr.list, Tx = Tx)
  
length(ccr.list)

basis <- ccr.list[[3]]$peaks.pos

ccr <- ccr.list[[3]]$ccr

all(is.na(ccr$ccf))

nseg <- 25

ccp.list <- lapply(1:nseg, function(nseg) {
  ccr.list <- lapply(1:seg_scale_max, function(seg_scale) {
    get_ccr_peaks(y, segments, seg_scale = seg_scale, nseg = nseg, npeaks = npeaks.set[seg_scale])
  })
  
  get_ccp(ccr.list, Tx = Tx)
})


ccr_peak <- get_ccr_peaks(y, segments, seg_scale = 3, nseg = nseg, npeaks = 1)

# expected peaks
ccr_peak$ccrpeaks$peaks

tt <- get_seg_scale(segments, nseg, scale = 3)
ccr <- get_ccf3(y, tt$aug_seg, 
                min.overlap = length(tt$aug_seg[!is.na(tt$aug_seg)]))

plot(ccr$lag - tmp_pos + 1, ccr$ccf, type = 'l')

system.time({
  ccr <- get_ccf3(y, tt$aug_seg, 
                  min.overlap = length(tt$aug_seg[!is.na(tt$aug_seg)]))
})
system.time({
  ccr2 <- get_ccf4(y, tt$aug_seg, 
                  min.overlap = length(tt$aug_seg[!is.na(tt$aug_seg)]))
})

system.time({
  tt <- get_seg_scale(segments, nseg, scale = 3)
  
  find_maxs <- rollapply(ccr$ccf, 3, function(x) max(x) == x[2], 
                         fill = list(NA, NA, NA))
  
  peaks <- which(find_maxs)
  peaks.heights <- ccr$ccf[peaks]
  
  npeaks <- 5
  tmp_pos <- tt$aug_idx[1]
  
  peak_pos <- sort(peaks[order(peaks.heights, decreasing = T)][1:npeaks] - tmp_pos)
  peak_height <- ccr$ccf[peak_pos + tmp_pos]
  
})


system.time({
  t2 <- get_ccr_peaks2(y, segments, seg_scale = 3, nseg = 14, npeaks = 5)
})

system.time({
  t1 <- get_ccr_peaks(y, segments, seg_scale = 3, nseg = 14, npeaks = 5)
})

t1$peaks.pos
t2$peaks.pos

t1$peaks.heights
t2$peaks.heights






















