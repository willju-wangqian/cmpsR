# loading packages
library(tidyverse)
library(x3ptools)
library(randomForest)
library(bulletxtrctr)
library(assertthat)
library(zoo)

source("func_collection.R")

# Following codes are obtained from the case study section of Chapter 3

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

bullets <- bullets %>% mutate(grooves = ccdata %>% purrr::map(.f = cc_locate_grooves, 
                                                              method = "middle", adjust = 30, return_plot = TRUE))

# codes for plotting
# do.call(gridExtra::grid.arrange, lapply(bullets$grooves, `[[`, 2))



# obtain the signiture
bullets <- bullets %>% mutate(sigs = purrr::map2(.x = ccdata, .y = grooves, 
                                                 .f = function(x, y) {
                                                   cc_get_signature(ccdata = x, grooves = y, span1 = 0.75, span2 = 0.03)
                                                 }))



bullets$bulletland <- paste0(bullets$bullet, "-", bullets$land)

lands <- unique(bullets$bulletland)




##############################################################
# CMPS Algorithm



comparisons.cmps <- data.frame(expand.grid(land1 = lands, land2 = lands), stringsAsFactors = FALSE)

comparisons.cmps <- comparisons.cmps %>% mutate(aligned = purrr::map2(.x = land1, .y = land2, 
                                                            .f = function(xx, yy) {
                                                              land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
                                                              land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
                                                              land1$bullet <- "first-land"
                                                              land2$bullet <- "second-land"
                                                              
                                                              sig_align(land1$sig, land2$sig)
                                                            }))


#########################################
system.time({
  comparisons.cmps <- comparisons.cmps %>% 
    mutate(cmps = aligned %>% purrr::map_dbl(.f = function(a) {
      extract_feature_cmps(a$lands$sig1, a$lands$sig2)
    }))
})

###
# user  system elapsed 
# 204.14    0.12  205.38 


comparisons.cmps %>% select(land1, land2, cmps)


comparisons <- data.frame(expand.grid(land1 = lands[1:6], land2 = lands[7:12]), stringsAsFactors = FALSE)

comparisons <- comparisons %>% mutate(aligned = purrr::map2(.x = land1, .y = land2, 
                                                            .f = function(xx, yy) {
                                                              land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
                                                              land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
                                                              land1$bullet <- "first-land"
                                                              land2$bullet <- "second-land"
                                                              
                                                              sig_align(land1$sig, land2$sig)
                                                            }))
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
