library(tidyverse)
library(x3ptools)
library(randomForest)
library(bulletxtrctr)
library(assertthat)
library(pracma)



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

for (i in 1:6) {
  x3p_write(b1$x3p[[i]], file = paste0("example_x3p_br1_blt1_l", i, ".x3p"))
}

for (i in 1:6) {
  x3p_write(b2$x3p[[i]], file = paste0("example_x3p_br1_blt2_l", i, ".x3p"))
}


b1$bullet <- 1
b2$bullet <- 2
b1$land <- 1:6
b2$land <- 1:6
bullets <- rbind(b1, b2)

bullets <- bullets %>% mutate(x3p = x3p %>% purrr::map(.f = x3p_m_to_mum)) %>% 
  mutate(x3p = x3p %>% purrr::map(.f = function(x) x %>% rotate_x3p(angle = -90) %>% 
                                    y_flip_x3p()))
# contrastive learning

bullets$land

x3ptools::x3p_image(bullets$x3p[[2]])
x3ptools::x3p_image(bullets$x3p[[9]])

for (i in 1:12) {
  x3p_write(bullets$x3p[[i]], 
            file = paste0("proc_example_x3p_br1_blt", bullets$bullet[i], "_l", bullets$land[i], ".x3p"))
}


bullets <- bullets %>% mutate(crosscut = x3p %>% purrr::map_dbl(.f = x3p_crosscut_optimize))

bullets <- bullets %>% mutate(ccdata = purrr::map2(.x = x3p, .y = crosscut, 
                                                   .f = x3p_crosscut))

crosscuts <- bullets %>% tidyr::unnest(ccdata)

ggplot(data = crosscuts, aes(x = x, y = value)) + 
  geom_line() + 
  facet_grid(bullet ~ land, labeller = "label_both") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = rel(0.9)))

bullets <- bullets %>% mutate(grooves = ccdata %>% purrr::map(.f = cc_locate_grooves, 
                                                              method = "middle", adjust = 30, return_plot = TRUE))

# do.call(gridExtra::grid.arrange, lapply(bullets$grooves, `[[`, 2))



# obtain the signiture
bullets <- bullets %>% mutate(sigs = purrr::map2(.x = ccdata, .y = grooves, 
                                                 .f = function(x, y) {
                                                   cc_get_signature(ccdata = x, grooves = y, span1 = 0.75, span2 = 0.03)
                                                 }))

signatures <- bullets %>% select(source, sigs) %>% tidyr::unnest()
bullet_info <- bullets %>% select(source, bullet, land)

signatures %>% filter(!is.na(sig), !is.na(raw_sig)) %>% 
  left_join(bullet_info, by = "source") %>% 
  ggplot(aes(x = x)) + 
  geom_line(aes(y = raw_sig), colour = "grey70") + 
  geom_line(aes(y = sig), colour = "grey30") + 
  facet_grid(bullet ~ land, labeller = "label_both") + 
  ylab("value") + ylim(c(-5, 5)) + theme_bw()

bullets$bulletland <- paste0(bullets$bullet, "-", bullets$land)
lands <- unique(bullets$bulletland)
comparisons <- data.frame(expand.grid(land1 = lands, land2 = lands), stringsAsFactors = FALSE)

comparisons <- comparisons %>% mutate(aligned = purrr::map2(.x = land1, .y = land2, 
                                                            .f = function(xx, yy) {
                                                              land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
                                                              land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
                                                              land1$bullet <- "first-land"
                                                              land2$bullet <- "second-land"
                                                              
                                                              sig_align(land1$sig, land2$sig)
                                                            }))

###########################
# sig_align?
land1 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
land2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]
land1$bullet <- "first-land"
land2$bullet <- "second-land"
aligned <- sig_align(land1$sig, land2$sig)

sig1 <- aligned$lands$sig1
sig2 <- aligned$lands$sig2


seg1 <- get_segs(sig1, 25)

seg1 <- split(sig1, ceiling(seq_along(sig1)/(length(sig1)/25)))

cors <- get_ccf(sig2, seg1[[1]], round(0.75 * min(length(seg1[[1]]), length(sig2))))
tt <- sig2[6:50]

cors <- get_ccf(seg1[[1]], tt)
get_ccf(seg1[[1]], sig2[6:50], round(0.75 * min(length(seg1[[1]]), length(sig2))))
cors


get_ccf(seg1, sig2, min.overlap = round(0.75 * 45))



cors %>% as.data.frame() %>% ggplot() +
  geom_line(aes(x=lag, y=ccf))

cors$lag[which.max(cors$ccf)]
cors$ccf[which.max(cors$ccf)]



aligned$lands$sig3 <- aligned$lands$sig2
extract_feature_ccf(aligned$lands)

get_ccf(land1$sig, land2$sig)

sig1 <- land1 %>% select(x, sig)
sig2 <- land2 %>% select(x, sig)

sig1 %>% ggplot() +
  geom_line(aes(x = x, y = sig))

sig2 %>% ggplot() +
  geom_line(aes(x = x, y = sig))

# rr <- get_ccf(sig2$sig, sig1$sig, round(0.75 * min(length(sig1), length(sig2))))
rr <- get_ccf(sig2$sig, sig1$sig)
# rr <- get_ccf(aligned$lands$sig1, aligned$lands$sig2, round(0.75 * min(length(sig1), length(sig2))))

which.max(rr$ccf)
rr$ccf[1352]
rr$lag[1352]


x <- runif(20)
get_ccf(x, lead(x, 5))
get_ccf(x, lag(x, 5), min.overlap = 3)
x <- runif(100)
get_ccf(x[45:50], x, min.overlap = 6)

get_ccf(x, x[45:50], min.overlap = 6)

cor(sig1$sig[1:500], sig2$sig[1:500])

###########################
subset(comparisons, land1 == "2-3" & land2 == "1-2")$aligned[[1]]$lands %>% 
  mutate(`b2-l3` = sig1, `b1-l2` = sig2) %>% select(-sig1, -sig2) %>% tidyr::gather(sigs, 
                                                                                    value, `b2-l3`, `b1-l2`) %>% ggplot(aes(x = x, y = value, colour = sigs)) + 
  geom_line() + theme_bw() + scale_color_brewer(palette = "Dark2")


############################
comparisons <- comparisons %>% 
  mutate(ccf0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_ccf(x$lands)), 
         lag0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_lag(x$lands)), 
         D0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_D(x$lands)), 
         length0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_length(x$lands)), 
         overlap0 = aligned %>% purrr::map_dbl(.f = function(x) extract_feature_overlap(x$lands)), 
         striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75), 
         cms_per_mm = purrr::map2(striae, aligned, .f = function(s, a) {
           extract_feature_cms_per_mm(s$lines, a$lands, resolution = 1.5625)
         }), 
         matches0 = striae %>% purrr::map_dbl(.f = function(s) {
           bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", 
                                                          match = TRUE)
         }), 
         mismatches0 = striae %>% purrr::map_dbl(.f = function(s) {
           bulletxtrctr:::extract_helper_feature_n_striae(s$lines, type = "peak", 
                                                          match = FALSE)
         }), 
         bulletA = gsub("([1-2])-([1-6])", "\\1", land1), 
         bulletB = gsub("([1-2])-([1-6])","\\1", land2), 
         landA = gsub("([1-2])-([1-6])", "\\2", land1), 
         landB = gsub("([1-2])-([1-6])", "\\2", land2))

comparisons <- comparisons %>% 
  mutate(features = purrr::map2(.x = aligned, .y = striae, 
                                .f = extract_features_all, resolution = 1.5625), 
         legacy_features = purrr::map(striae, extract_features_all_legacy, resolution = 1.5625)) %>% 
  tidyr::unnest(legacy_features)

comparisons %>% ggplot(aes(x = landA, y = landB, fill = ccf)) + geom_tile() + 
  scale_fill_gradient2(low = "grey80", high = "darkorange", midpoint = 0.5) + 
  facet_grid(bulletB ~ bulletA, labeller = "label_both") + xlab("Land A") + 
  ylab("Land B") + theme(aspect.ratio = 1)

comparisons$rfscore <- predict(bulletxtrctr::rtrees, newdata = comparisons, 
                               type = "prob")[, 2]

comparisons %>% ggplot(aes(x = landA, y = landB, fill = rfscore)) + geom_tile() + 
  scale_fill_gradient2(low = "grey80", high = "darkorange", midpoint = 0.5) + 
  facet_grid(bulletB ~ bulletA, labeller = "label_both") + xlab("Land A") + 
  ylab("Land B") + theme(aspect.ratio = 1)
