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

plot(ccr$lag, ccr$ccf, type = 'l')


find_maxs <- rollapply(ccr$ccf, 3, function(x) max(x) == x[2], 
                       fill = list(NA, NA, NA))
peaks <- which(find_maxs)
peaks.heights <- ccr$ccf[peaks]
