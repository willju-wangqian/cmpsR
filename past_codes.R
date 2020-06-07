###############################
# Jun. 6th
###############################

# This is a demo to see how we can modify and utilize get_ccf

library(tidyverse)
library(x3ptools)
library(randomForest)
library(bulletxtrctr)
library(assertthat)
library(pracma)

# Following codes are obtained from the case study section of Chapter 3
# codes for plotting are all commented out

# Case Study
# bullet 1
urllist1 <- c("https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/cd204983-465b-4ec3-9da8-cba515a779ff", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/0e72228c-5e39-4a42-8c4e-3da41a11f32c", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/b9d6e187-2de7-44e8-9b88-c83c29a8129d", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/fda92f6a-71ba-4735-ade0-02942d14d1e9", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/8fa798b4-c5bb-40e2-acf4-d9296865e8d4", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/81e817e5-15d8-409f-b5bd-d67c525941fe")
# bullet 2
urllist2 <- c("https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/288341e0-0fdf-4b0c-bd26-b31ac8c43f72", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/c97ada55-3a35-44fd-adf3-ac27dd202522", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/8a1805d9-9d01-4427-8873-aef4a0bd323a", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/a116e448-18e1-4500-859c-38a5f5cc38fd", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/0b7182d3-1275-456e-a9b4-ae378105e4af", 
              "https://tsapps.nist.gov/NRBTD/Studies/BulletMeasurement/DownloadMeasurement/86934fcd-7317-4c74-86ae-f167dbc2f434")

b1 <- read_bullet(urllist = urllist1)
b2 <- read_bullet(urllist = urllist2)

b1$bullet <- 1
b2$bullet <- 2
b1$land <- 1:6
b2$land <- 1:6
bullets <- rbind(b1, b2)

bullets <- bullets %>% mutate(x3p = x3p %>% purrr::map(.f = x3p_m_to_mum)) %>% 
  mutate(x3p = x3p %>% purrr::map(.f = function(x) x %>% rotate_x3p(angle = -90) %>% 
                                    y_flip_x3p()))

bullets <- bullets %>% mutate(crosscut = x3p %>% purrr::map_dbl(.f = x3p_crosscut_optimize))

bullets <- bullets %>% mutate(ccdata = purrr::map2(.x = x3p, .y = crosscut, 
                                                   .f = x3p_crosscut))

# crosscuts <- bullets %>% tidyr::unnest(ccdata)

# ggplot(data = crosscuts, aes(x = x, y = value)) + 
#   geom_line() + 
#   facet_grid(bullet ~ land, labeller = "label_both") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = rel(0.9)))

bullets <- bullets %>% mutate(grooves = ccdata %>% purrr::map(.f = cc_locate_grooves, 
                                                              method = "middle", adjust = 30, return_plot = TRUE))

# do.call(gridExtra::grid.arrange, lapply(bullets$grooves, `[[`, 2))



# obtain the signiture
bullets <- bullets %>% mutate(sigs = purrr::map2(.x = ccdata, .y = grooves, 
                                                 .f = function(x, y) {
                                                   cc_get_signature(ccdata = x, grooves = y, span1 = 0.75, span2 = 0.03)
                                                 }))

# signatures <- bullets %>% select(source, sigs) %>% tidyr::unnest()
# bullet_info <- bullets %>% select(source, bullet, land)

# signatures %>% filter(!is.na(sig), !is.na(raw_sig)) %>%
#   left_join(bullet_info, by = "source") %>% 
#   ggplot(aes(x = x)) + 
#   geom_line(aes(y = raw_sig), colour = "grey70") + 
#   geom_line(aes(y = sig), colour = "grey30") + 
#   facet_grid(bullet ~ land, labeller = "label_both") + 
#   ylab("value") + ylim(c(-5, 5)) + theme_bw()

bullets$bulletland <- paste0(bullets$bullet, "-", bullets$land)

land1 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
land2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]
land1$bullet <- "first-land"
land2$bullet <- "second-land"
aligned <- sig_align(land1$sig, land2$sig)


##############################################################
# demo starts here
# please load all functions in func_collection.R first

segments <- get_segs(aligned$lands$sig1, 25)
y <- aligned$lands$sig2

nseg <- 15

ccr <- get_ccf3(y, segments$segs[[nseg]], min.overlap = length(segments$segs[[nseg]]))
fp <- findpeaks(ccr$ccf[!is.na(ccr$ccf)], npeaks = 5,
                threshold = 0, sortstr = T, nups = 0, ndowns = 1)
pos1 <- ccr$lag - segments$index[[nseg]][1] + 1

# level 2
tt <- get_seg_level(segments, nseg, level = 2)
ccr2 <- get_ccf3(y, tt$aug_seg, min.overlap = length(tt$aug_seg))
fp2 <- findpeaks(ccr2$ccf[!is.na(ccr2$ccf)], npeaks = 5, 
                 threshold = 0, sortstr = T, nups = 0, ndowns = 1)
pos2 <- ccr2$lag - tt$aug_idx[1] + 1

# level 3
tt3 <- get_seg_level(segments, nseg, level = 3)
ccr3 <- get_ccf3(y, tt3$aug_seg, min.overlap = length(tt3$aug_seg))
fp3 <- findpeaks(ccr3$ccf[!is.na(ccr3$ccf)], npeaks = 5, 
                 threshold = 0, sortstr = T, nups = 0, ndowns = 1)
pos3 <- ccr3$lag - tt3$aug_idx[1] + 1

sort(pos1[!is.na(ccr$ccf)][fp[,2]])
sort(pos2[!is.na(ccr2$ccf)][fp2[,2]])
sort(pos3[!is.na(ccr3$ccf)][fp3[,2]])

plot(pos1, ccr$ccf, type = 'l') 
plot(pos2, ccr2$ccf, type = 'l')
plot(pos3, ccr3$ccf, type = 'l')

rr <- sig_get_peaks(ccr2$ccf, smoothfactor = 1, window = F, striae = F)
rr$extrema[order(rr$peaks)[1:5]]

rr$dframe$x[order(rr$dframe$smoothed[rr$peaks])[1:5]]

########
# this one works
rr$peaks[order(rr$peaks.heights,decreasing = T)][1:5] - tt$aug_idx[1]

rr
sort(ccr$lag[!is.na(ccr$ccf)][fp[,2]])
rr$peaks

corrr <- get_ccf3(y, segments$segs[[nseg]], min.overlap = length(segments$segs[[nseg]]))
corrr$lag
plot(corrr$lag - segments$index[[nseg]][1], corrr$ccf, type = 'l')


head(lag(y, 10), 10)
