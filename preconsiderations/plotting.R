library(CMPS)
library(ggplot2)
library(tidyverse)
theme_set(theme_bw())

data("bullets")
land2_3 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
land1_2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]

x <- land2_3$sig
y <- land1_2$sig

cmps <- extract_feature_cmps(x, y, full_result = T)
cmps$CMPS.score

rbind(data.frame(value = x, index = 1:length(x), sig = "x"),
      data.frame(value = y, index = 1:length(y), sig = "y")) %>% 
  ggplot() + 
  geom_line(aes(index, value, color = sig)) +
  ylab("signature") +
  xlab("index") +
  ggtitle("plot of x and y")
ggsave("image/step0.png")





segments <- get_segs(x, len = 50)

nseg <- length(segments$segs)

df <- data.frame(value = unlist(segments$segs),
                 segs = rep(1:nseg, each = 50),
                 index = unlist(segments$index))

df.part <- data.frame(value = unlist(segments$segs),
                      segs = 0,
                      index = unlist(segments$index))
df <- rbind(df, df.part)

cutt <- sapply(segments$index, function(idx) {idx[1]})

df %>% filter(segs <= 4) %>% 
  ggplot() +
  geom_line(aes(x = index, y = value)) +
  geom_vline(xintercept = cutt, color = "red") +
  facet_grid(segs ~ .) +
  xlab("position") +
  ylab("signature")
ggsave("image/step1_1.png")

seg3 <- df %>% filter(segs == 7)
seg3$segs <- "segment 7"
df.y <- data.frame(value = y, segs = "y", index = 1:length(y))
df.y <- rbind(df.y, seg3)
df.y %>% filter(!is.na(value)) %>% 
  ggplot() +
  geom_line(aes(x = index, y = value)) +
  facet_grid(segs ~ .) +
  xlab("position") +
  ylab("signature") +
  ggtitle("Comparison Signature y and The Third Segment")
ggsave("image/step2_1.png")



ccrpeaks <- get_ccr_peaks(df.y$value, segments = segments, nseg = 7, npeaks = 1, seg_scale = 1)
df.ccf <- data.frame(value = ccrpeaks$ccr$ccf, index = ccrpeaks$adj.pos)
# tt <- get_ccf4(df.y$value, seg3$value, 50)
# tail(tt$ccf)

df.ccf %>% ggplot() +
  geom_line(aes(index, value)) + 
  geom_vline(xintercept = ccrpeaks$peaks.pos, color = "red") +
  xlim(-400, 800) +
  xlab("position") +
  ylab("ccf") +
  ggtitle("CCF of y and the 7th segment")
ggsave("image/step2_2.png")

##################################################
# step 3-1
comp <- x
npeaks <- 1
seg_scale <- 1

ccf.df <- do.call(rbind, lapply(6:12, function(nss) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, nseg = nss, npeaks = npeaks, seg_scale = seg_scale)
  df.ccf <- data.frame(value = ccrpeaks$ccr$ccf, index = ccrpeaks$adj.pos)
  df.ccf$segs <- nss
  df.ccf
}))
peak.df <- do.call(rbind, lapply(6:12, function(nss) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, nseg = nss, npeaks = npeaks, seg_scale = seg_scale)
  df.ccf <- data.frame(value = ccrpeaks$peaks.heights, index = ccrpeaks$peaks.pos)
  df.ccf$segs <- nss
  df.ccf
}))

ccf.df %>% 
  ggplot() +
  geom_line(aes(x = index, y = value)) +
  geom_point(data = peak.df, aes(x = index, y = value),
             color = "blue") +
  geom_vline(xintercept = 0, color = "red") +
  facet_grid(segs ~ .) +
  # ggtitle("Real Case: x compares to y") +
  ggtitle("Ideal Case: x compares to itself") +
  xlab("position") +
  ylab("ccf")
ggsave("image/step3_1.png")

###########################################################
# step 3-2
comp <- y
npeaks <- 5
seg_scale <- 1

ccf.df <- do.call(rbind, lapply(6:12, function(nss) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, nseg = nss, npeaks = npeaks, seg_scale = seg_scale)
  df.ccf <- data.frame(value = ccrpeaks$ccr$ccf, index = ccrpeaks$adj.pos)
  df.ccf$segs <- nss
  df.ccf
}))
peak.df <- do.call(rbind, lapply(6:12, function(nss) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, nseg = nss, npeaks = npeaks, seg_scale = seg_scale)
  df.ccf <- data.frame(value = ccrpeaks$peaks.heights, index = ccrpeaks$peaks.pos)
  df.ccf$segs <- nss
  df.ccf
}))

ccf.df %>% 
  ggplot() +
  geom_line(aes(x = index, y = value)) +
  geom_point(data = peak.df, aes(x = index, y = value),
             color = "blue") +
  geom_vline(xintercept = -6, color = "red") +
  facet_grid(segs ~ .) +
  ggtitle("Real Case: x compares to y") +
  # ggtitle("Ideal Case: x compares to itself") +
  xlab("position") +
  ylab("ccf")
ggsave("image/step3_2.png")


###############################################
# step 5
comp <- y
nseg <- 7
npeaks.set <- c(5, 3, 1)

multi.seg <- do.call(rbind, lapply(1:3, function(scale) {
  tt <- get_seg_scale(segments, nseg, scale = scale)
  df.tmp <- data.frame(value = tt$aug_seg, index=tt$aug_idx, scale=paste("scale", scale))
}))

rbind(data.frame(value = x, index = 1:length(x), scale = "x"), multi.seg) %>% filter(!is.na(value)) %>% 
  ggplot() +
  geom_line(aes(index, value)) +
  facet_grid(scale ~ .) + 
  ylab("signature") +
  xlab("position") +
  ggtitle("x and segment 7 in 3 different scales")
ggsave("image/step5_1.png")

multi.df <- do.call(rbind, lapply(1:3, function(scale) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, nseg = nseg, npeaks = npeaks.set[scale], seg_scale = scale)
  df.tmp <- data.frame(value = ccrpeaks$ccr$ccf, index = ccrpeaks$adj.pos)
  df.tmp$scale <- paste("scale", scale)
  df.tmp
}))
multi.peak.df <- do.call(rbind, lapply(1:3, function(scale) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, nseg = nseg, npeaks = npeaks.set[scale], seg_scale = scale)
  df.tmp <- data.frame(value = ccrpeaks$peaks.heights, index = ccrpeaks$peaks.pos)
  df.tmp$scale <- paste("scale", scale)
  df.tmp
}))
multi.df %>% ggplot() +
  geom_line(aes(x=index, y=value)) +
  geom_point(data = multi.peak.df, aes(x = index, y = value),
           color = "blue") +
  geom_vline(xintercept = -7, color = "red") +
  facet_grid(scale ~ .) +
  scale_x_continuous(breaks=seq(-400,800,100)) +
  ylab("ccf") +
  xlab("position") +
  ggtitle("ccf of segment 7 and y in 3 different scales")
ggsave("image/step5_2.png")


