land1 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
land2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]
land1$bullet <- "first-land"
land2$bullet <- "second-land"
aligned <- sig_align(land1$sig, land2$sig)

nseg <- 25

segments <- get_segs(aligned$lands$sig1, nseg)
y <- aligned$lands$sig2

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






















