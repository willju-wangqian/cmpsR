library(bulletxtrctr)
library(tidyverse)
library(x3ptools)
library(CMPS)

data("bullets")
# land2_3 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
# land1_2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]

# compute cmps

# # algorithm with multi-peak insepction at three different segment scales
# cmps_with_multi_scale <- extract_feature_cmps(land2_3$sig, land1_2$sig, full_result = TRUE)
# 
# # algorithm with multi-peak inspection at the basis scale only
# cmps_without_multi_scale <- extract_feature_cmps(land2_3$sig, land1_2$sig, seg_scale_max = 1, 
#                                                  npeaks.set = 5, full_result = TRUE)
## Not run: 
library(tidyverse)
library(bulletxtrctr)

lands <- unique(bullets$bulletland)

comparisons <- data.frame(expand.grid(land1 = lands[1:6], land2 = lands[7:12]), 
                          stringsAsFactors = FALSE)

comparisons <- comparisons %>% mutate(
  aligned = purrr::map2(.x = land1, .y = land2, 
                        .f = function(xx, yy) {
                          land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
                          land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
                          land1$bullet <- "first-land"
                          land2$bullet <- "second-land"
                          
                          sig_align(land1$sig, land2$sig)
                        }))

comparisons <- comparisons %>% 
  mutate(cmps = aligned %>% purrr::map(.f = function(a) {
    extract_feature_cmps(a$lands$sig1, a$lands$sig2, full_result = TRUE)
  }))

# comparisons.cmps <- comparisons.cmps %>% 
#   mutate(cmps = aligned %>% purrr::map_dbl(.f = function(a) {
#     extract_feature_cmps(a$lands$sig1, a$lands$sig2, full_result = FALSE)
#   }))
# comparisons.cmps %>% select(land1, land2, cmps) 

comparisons <- comparisons %>% 
  mutate(
    cmps_score = sapply(comparisons$cmps, function(x) x$CMPS.score),
    cmps_nseg = sapply(comparisons$cmps, function(x) x$nseg)
  )

comparisons <- comparisons %>%
  mutate(
    bulletA = gsub("(\\d)-\\d", "\\1", land1),
    landA = gsub("\\d-(\\d)", "\\1", land1),
    bulletB = gsub("(\\d)-\\d", "\\1", land2),
    landB = gsub("\\d-(\\d)", "\\1", land2)
  )

comparisons %>%   
  group_by(bulletA, bulletB) %>% tidyr::nest() %>%
  mutate(
    cmps_max_bar = data %>% purrr::map_dbl(
      .f = function(d) max(compute_average_scores(land1 = d$landA,
                                                  land2 = d$landB,
                                                  d$cmps_score)))
  )

comparisons %>% select(-aligned, -cmps)

comparisons %>% ggplot(aes(x = landB, y = landA, fill = cmps_score)) + geom_tile() + 
  scale_fill_gradient2(low = "gray80", high = "darkorange", midpoint = 6) + 
  facet_grid(bulletA ~ bulletB, labeller = "label_both") + xlab("Land B") + 
  ylab("Land A") + theme(aspect.ratio = 1)

# cp1 <- comparisons %>% select(land1, land2, cmps_score, cmps_nseg)
# cp1
# 
# compute_average_scores(cp1$land1, cp1$land2, cp1$cmps_score)
# 
# # Set the data up to be read in, cleaned, etc.
# library(bulletxtrctr)
# library(x3ptools)
# 
# bullets <- bullet_pipeline(
#   location = list(
#     Bullet1 = c(hamby252demo$bullet1),
#     Bullet2 = c(hamby252demo$bullet2)
#   ),
#   x3p_clean = function(x) x %>%
#     x3p_scale_unit(scale_by=10^6) %>%
#     rotate_x3p(angle = -90) %>%
#     y_flip_x3p()
# ) %>%
#   mutate(land = paste0(rep(1:2, each = 6), "-", rep(1:6, times = 2)))
# 
# comparisons <- data.frame(
#   expand.grid(land1 = bullets$land, land2 = bullets$land),
#   stringsAsFactors = FALSE)
# comparisons <- comparisons %>%
#   mutate(
#     aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
#       land1 <- bullets$sigs[bullets$land == xx][[1]]
#       land2 <- bullets$sigs[bullets$land == yy][[1]]
#       land1$bullet <- "first-land"
#       land2$bullet <- "second-land"
#       
#       sig_align(land1$sig, land2$sig)
#     }),
#     striae = purrr::map(aligned, sig_cms_max),
#     features = purrr::map2(.x = aligned, .y = striae, extract_features_all),
#     rfscore = purrr::map_dbl(features, rowMeans) # This is a hack until the new RF is fit...
#   )
# 
# # Clean up a bit
# comparisons <- comparisons %>%
#   mutate(
#     bulletA = gsub("(\\d)-\\d", "\\1", land1),
#     landA = gsub("\\d-(\\d)", "\\1", land1),
#     bulletB = gsub("(\\d)-\\d", "\\1", land2),
#     landB = gsub("\\d-(\\d)", "\\1", land2)
#   ) %>%
#   
# comparisons %>%   
#   group_by(bulletA, bulletB) %>% tidyr::nest() %>%
#   mutate(
#     bullet_score = data %>% purrr::map_dbl(
#       .f = function(d) max(compute_average_scores(land1 = d$landA,
#                                                   land2 = d$landB,
#                                                   d$rfscore)))
#   )
# 
# cp1 <- cp1 %>% mutate(
#   bulletA = gsub("(\\d)-\\d", "\\1", land1),
#   landA = gsub("\\d-(\\d)", "\\1", land1),
#   bulletB = gsub("(\\d)-\\d", "\\1", land2),
#   landB = gsub("\\d-(\\d)", "\\1", land2)
# )
# 
# max(compute_average_scores(cp1$landA, cp1$landB, cp1$cmps_score))














