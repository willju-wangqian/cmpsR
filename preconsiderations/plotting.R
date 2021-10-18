library(cmpsR)
library(ggplot2)
library(tidyverse)
theme_set(theme_bw())

data("bullets")
land2_3 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
land1_2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]

x <- land2_3$sig
y <- land1_2$sig

cmps <- extract_feature_cmps(x, y, include = "full")
cmps$CMPS.score

rbind(data.frame(value = x, index = 1:length(x), sig = "x"),
      data.frame(value = y, index = 1:length(y), sig = "y")) %>% 
  ggplot() + 
  geom_line(aes(index, value, color = sig)) +
  ylab("signature") +
  xlab("index") +
  ggtitle("plot of x and y")
# ggsave("man/figures/step0.png")





segments <- get_segs(x, len = 50)

nseg <- length(segments$segs)

df <- data.frame(value = unlist(segments$segs),
                 segs = rep(1:nseg, each = 50),
                 index = unlist(segments$index))
df$segs_tag <- paste("seg", df$segs)

df.partx <- data.frame(value = unlist(segments$segs),
                       segs = 0,
                       index = unlist(segments$index),
                       segs_tag = "x")
df.party <- data.frame(value = y,
                       segs = 0,
                       index = 1:length(y),
                       segs_tag = "y")


cutt <- sapply(segments$index, function(idx) {idx[1]})

rbind(df.partx, df) %>% filter(segs <= 4) %>% 
  ggplot() +
  geom_line(aes(x = index, y = value)) +
  geom_vline(xintercept = cutt, color = "red") +
  facet_grid(segs_tag ~ .) +
  xlab("position") +
  ylab("signature") +
  ggtitle("Cut x into basis segments")
ggsave("man/figures/step1_1.png")

rbind(df.party, df.partx, df) %>% filter(segs %in% c(0, 1,2,6,7) & !is.na(value)) %>% 
  ggplot() +
  geom_line(aes(x = index, y = value)) +
  geom_vline(xintercept = cutt, color = "red") +
  facet_grid(segs_tag ~ .) +
  xlab("position") +
  ylab("signature") +
  ggtitle("Compare y to basis segments of x")
ggsave("man/figures/step1_2.png")

seg3 <- df %>% filter(segs == 7)
seg3$segs <- "segment 7"
df.y <- data.frame(value = y, segs = "y", index = 1:length(y))
df.y <- rbind(df.y, seg3[,-4])
df.y$segs <- factor(df.y$segs, levels = c("y", "segment 7"))
df.y %>% filter(!is.na(value)) %>% 
  ggplot() +
  geom_line(aes(x = index, y = value)) +
  facet_grid(segs ~ .) +
  xlab("position") +
  ylab("signature") +
  ggtitle("Reference Signature y and the 7th Segment")
# ggsave("man/figures/step2_1.png")



ccrpeaks <- get_ccr_peaks(y, segments = segments, 50, nseg = 7, npeaks = 1)
df.ccf <- data.frame(value = ccrpeaks$ccr$ccf, index = ccrpeaks$adj.pos)
# tt <- get_ccf4(df.y$value, seg3$value, 50)
# tail(tt$ccf)

df.ccf %>% ggplot() +
  geom_line(aes(index, value)) + 
  geom_vline(xintercept = ccrpeaks$peaks.pos, color = "red") +
  xlab("position") +
  ylab("ccf") +
  ggtitle("CCF of y and the 7th segment")
# ggsave("man/figures/step2_2.png")

##################################################
# step 3-1
comp <- x
npeaks <- 1
seg_outlength <- 50

ccf.df <- do.call(rbind, lapply(6:12, function(nss) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, seg_outlength = seg_outlength, nseg = nss, npeaks = npeaks)
  df.ccf <- data.frame(value = ccrpeaks$ccr$ccf, index = ccrpeaks$adj.pos)
  df.ccf$segs <- nss
  df.ccf
}))
peak.df <- do.call(rbind, lapply(6:12, function(nss) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, seg_outlength = seg_outlength, nseg = nss, npeaks = npeaks)
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
  ggtitle("Ideal Case: x Compares to Itself") +
  xlab("position") +
  ylab("ccf")
# ggsave("man/figures/step3_1.png")

###########################################################
# step 3-2
comp <- y
npeaks <- 5
seg_outlength <- 50

ccf.df <- do.call(rbind, lapply(6:12, function(nss) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, seg_outlength = seg_outlength, nseg = nss, npeaks = npeaks)
  df.ccf <- data.frame(value = ccrpeaks$ccr$ccf, index = ccrpeaks$adj.pos)
  df.ccf$segs <- nss
  df.ccf
}))
peak.df <- do.call(rbind, lapply(6:12, function(nss) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, seg_outlength = seg_outlength, nseg = nss, npeaks = npeaks)
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
# ggsave("man/figures/step3_2.png")


###############################################
# step 5
comp <- y
nseg <- 7
npeaks.set <- c(5, 3, 1)

multi.seg <- do.call(rbind, lapply(1:3, function(scale) {
  tt <- get_seg_scale(segments, nseg, scale = scale)
  df.tmp <- data.frame(value = tt$aug_seg, index=tt$aug_idx, scale=paste("scale", scale))
}))

# rbind(data.frame(value = x, index = 1:length(x), scale = "x"), multi.seg) 
multi.seg %>% filter(!is.na(value)) %>% 
  ggplot() +
  geom_line(aes(index, value)) +
  facet_grid(scale ~ .) + 
  ylab("signature") +
  xlab("position") +
  ggtitle("x and segment 7 in 3 different scales")
# ggsave("man/figures/step5_1.png")

out_length <- c(50, 100, 200)

multi.df <- do.call(rbind, lapply(1:3, function(scale) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments,
                            seg_outlength = out_length[scale],
                            nseg = nseg, npeaks = npeaks.set[scale])
  df.tmp <- data.frame(value = ccrpeaks$ccr$ccf, index = ccrpeaks$adj.pos)
  df.tmp$scale <- paste("scale level", scale)
  df.tmp
}))
multi.peak.df <- do.call(rbind, lapply(1:3, function(scale) {
  ccrpeaks <- get_ccr_peaks(comp, segments = segments, seg_outlength = out_length[scale],
                            nseg = nseg, npeaks = npeaks.set[scale])
  df.tmp <- data.frame(value = ccrpeaks$peaks.heights, index = ccrpeaks$peaks.pos)
  df.tmp$scale <- paste("scale level", scale)
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
# ggsave("man/figures/step5_2.png")


