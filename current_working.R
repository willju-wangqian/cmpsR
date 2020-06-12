# loading packages
library(tidyverse)
library(x3ptools)
library(randomForest)
library(bulletxtrctr)
library(assertthat)

source("func_collection.R")

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

# codes for plotting

# crosscuts <- bullets %>% tidyr::unnest(ccdata)
# ggplot(data = crosscuts, aes(x = x, y = value)) + 
#   geom_line() + 
#   facet_grid(bullet ~ land, labeller = "label_both") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = rel(0.9)))

bullets <- bullets %>% mutate(grooves = ccdata %>% purrr::map(.f = cc_locate_grooves, 
                                                              method = "middle", adjust = 30, return_plot = TRUE))

# codes for plotting
# do.call(gridExtra::grid.arrange, lapply(bullets$grooves, `[[`, 2))



# obtain the signiture
bullets <- bullets %>% mutate(sigs = purrr::map2(.x = ccdata, .y = grooves, 
                                                 .f = function(x, y) {
                                                   cc_get_signature(ccdata = x, grooves = y, span1 = 0.75, span2 = 0.03)
                                                 }))

# codes for plotting
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
# CMPS Algorithm
nseg <- 25

segments <- get_segs(aligned$lands$sig1, nseg)
y <- aligned$lands$sig2

seg_scale_max <- 3
npeaks.set <- c(5, 3, 1)

system.time({
  ccp.list <- lapply(1:nseg, function(nseg) {
    ccr.list <- lapply(1:seg_scale_max, function(seg_scale) {
      get_ccr_peaks(y, segments, seg_scale = seg_scale, nseg = nseg, npeaks = npeaks.set[seg_scale])
    })
    
    get_ccp(ccr.list)
  })
})


ccp.list.one <- lapply(1:nseg, function(nseg) {
  ccr <- get_ccr_peaks(y, segments, seg_scale = 1, nseg = nseg, npeaks = 5)
  ccr$peaks.pos
})

cmps <- get_CMPS(ccp.list, Tx = 25)
cmps$pos.df %>% head()
cmps$rec.position

extract_feature_cmps <- function(x, y, nseg = 25, seg_scale_max = 3, Tx = 25, npeaks.set = c(5, 3, 1),
                                 full_result = FALSE) {
  if (length(npeaks.set) != seg_scale_max) { 
    print("Need to specify the number of peaks for each segment scale.")
    return(NULL)
  }
  
  segments <- get_segs(x, nseg)
  
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
  if(full_result) { return(cmps) } 
  else { return(cmps$CMPS.score) }
}

land1.name <- unique(bullets$bulletland)[1:6]
land2.name <- unique(bullets$bulletland)[7:12]

extract_feature_cmps(aligned$lands$sig1, aligned$lands$sig2)

comparisons.cmps <- data.frame(expand.grid(land1 = land1.name, land2 = land2.name), stringsAsFactors = FALSE)

comparisons.cmps <- comparisons.cmps %>% mutate(aligned = purrr::map2(.x = land1, .y = land2, 
                                                            .f = function(xx, yy) {
                                                              land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
                                                              land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
                                                              land1$bullet <- "first-land"
                                                              land2$bullet <- "second-land"
                                                              
                                                              sig_align(land1$sig, land2$sig)
                                                            }))


cmps.collect <- rep(-1, 36)

for(i in 1:36) {
  s <- paste("comparing ", comparisons.cmps$land1[i], " and ", comparisons.cmps$land2[i], ", loop ", i, sep = '')
  print(s)
  cmps.collect[i] <- extract_feature_cmps(comparisons.cmps$aligned[[i]]$lands$sig1, comparisons.cmps$aligned[[i]]$lands$sig2)
}



comparisons.cmps %>% select(land1, land2) %>% slice(6)
comparisons.cmps %>% select(land1, land2) %>% slice(14)

comparisons.cmps$cmps <- cmps.collect
comparisons.cmps %>% select(land1, land2, cmps)

# system.time({
#   
#   comparisons.cmps <- comparisons.cmps %>% 
#   mutate(cmps.score = aligned %>% purrr::map_dbl(.f = function(a) {
#     extract_feature_cmps(a$lands$sig1, a$lands$sig2)
#   }))
# 
# })

# [1] "comparing 1-1 and 2-1, loop 1"
# [1] "comparing 1-2 and 2-1, loop 2"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-3 and 2-1, loop 3"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-4 and 2-1, loop 4"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-5 and 2-1, loop 5"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-6 and 2-1, loop 6"
# [1] "comparing 1-1 and 2-2, loop 7"
# [1] "comparing 1-2 and 2-2, loop 8"
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-3 and 2-2, loop 9"
# [1] "comparing 1-4 and 2-2, loop 10"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-5 and 2-2, loop 11"
# [1] "comparing 1-6 and 2-2, loop 12"
# [1] "comparing 1-1 and 2-3, loop 13"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-2 and 2-3, loop 14"
# [1] "comparing 1-3 and 2-3, loop 15"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-4 and 2-3, loop 16"
# [1] "comparing 1-5 and 2-3, loop 17"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-6 and 2-3, loop 18"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-1 and 2-4, loop 19"
# [1] "the length of the highest level should be 1."
# [1] "comparing 1-2 and 2-4, loop 20"
# [1] "the length of the highest level should be 1."
# [1] "the length of the highest level should be 1."





