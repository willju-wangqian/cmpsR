library(tidyverse)
library(x3ptools)
library(bulletxtrctr)
library(CMPS)

manual <- read_csv("D:\\Research/Bullet/grooves-manual-id/grooves_manual_hamby44.csv")

manual.m <- manual %>% 
  group_by(scan_id) %>% 
  summarise(
    crosscut.m = mean(crosscut),
    groove.lm = mean(groove_left_manual), 
    groove.rm = mean(groove_right_manual)
  ) %>%
  nest(grooves = c(groove.lm, groove.rm))
# %>% 
#   select(scan_id, crosscut.m, groove.lm, groove.rm) %>% group_by(scan_id) %>% 
#   mutate(
#     grooves = list(groove = c(groove.lm, groove.rm), b = NULL)
#   ) %>% ungroup()

manual.m$grooves <- lapply(manual.m$grooves, function(g) {
  list(groove = c(g$groove.lm, g$groove.rm))
})


