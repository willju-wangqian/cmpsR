
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CMPS

<!-- badges: start -->

<!-- badges: end -->

The CMPS package is an implementation of the Congruent Matching Profile
Segments (CMPS) method \[cite?\]. In general, it can be used for
objective comparison of striated tool marks, but in our examples, we
mainly use it for bullet profiles/signatures comparison. The CMPS number
is expected to be large if two signatures are similar. So it can also be
considered as a feature that measures the similarity of two bullets.

## Installation

You can install the released version of CMPS from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("CMPS")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("willju-wangqian/CMPS")
```

## Example

This is a basic example which shows you how to compute the CMPS number
of a comparison between two bullet signatures.

In this example, we are taking the signature of the third land engraved
area (LEA) of the second bullet as the reference signature and the
second LEA of the first bullet as the comparison signature. This is a KM
comparison.

``` r
library(CMPS)

data("bullets")
land23 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
land12 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]

x <- land23$sig
y <- land12$sig

cmps <- extract_feature_cmps(x, y, full_result = TRUE)
cmps
#> $CMPS.score
#> [1] 14
#> 
#> $rec.position
#> [1] -5
#> 
#> $pos.df
#>   position cmps
#> 1       -9   14
#> 2       -8   14
#> 3       -7   14
#> 4       -6   14
#> 5       -5   14
#> 6       -4   14
#> 7       -3   14
#> 8       -2   14
#> 9       -1   14
#> 
#> $nseg
#> [1] 22
```

Plot of x and y

<img src="image/step0.png" width="100%" />

## Main Idea

The main idea of the CMPS method is that:

1.  we take the first signature as the reference signature (`x`) and cut
    it into consecutive and non-overlapping basis segments of the same
    length. In this case, we have 22 basis segments in total.

<img src="image/step1_1.png" width="80%" />

2.  for each basis segment, we compute the cross-correlation function
    (ccf) between the basis segment and the comparison signature (`y`)

<img src="image/step2_1.png" width="80%" /><img src="image/step2_2.png" width="80%" />

  - for the `ccf` curve, the `position` represents the shift of the
    segment. A negative value means a shift to the left, and a positive
    value means a shift to the right;
  - we are interested in the peaks in the ccf curve and the positions of
    those peaks (as indicated by the red vertical line in the plot
    above). In other words, if we shift the segment, which position
    would give us the “best fit”?

<!-- end list -->

3.  ideally, if two signatures are identical, we are expecting the
    position of the highest peak in the ccf curve remains the same
    across all ccf curves;

<!-- end list -->

  - ideal case: compare `x` to itself. The highest peak has value 1 and
    is marked by the blue dot.

<img src="image/step3_1.png" width="80%" />

  - real case: compare `x` to `y`. We mark the 5 highest peaks for each
    ccf curve because the position of the “highest peak” might not be
    the best one.

<img src="image/step3_2.png" width="80%" />

4.  each ccf curve votes for 5 candidate positions, then we ask two
    questions in order to obtain the CMPS number/score:

<!-- end list -->

  - which position receives the most votes? -\> the best position
    (indicated by the red vertical line)

  - how many segments have voted for the best position? -\> CMPS number
    
    If we focus on these 7 segments only, the CMPS number is 6.

<!-- end list -->

5.  how can the segment vote more wisely?

<!-- end list -->

  - by increasing the segment length (scale), one can reduce the number
    of “false positive” peaks

<img src="image/step5_1.png" width="80%" /><img src="image/step5_2.png" width="80%" />

  - we choose 5 peaks at scale 1; 3 peaks at scale 2; 1 peak at scale 3
    
    the peak shared by all three scales is a consistent correlation peak
    (ccp). And the position of the ccp is our best choice. This is
    called a “multi segment lengths” strategy.

  - and then, each ccf curve votes for only 1 best condidate position if
    a ccp can be found; again, we ask two quesitons:
    
      - which position receives the most votes?
      - how many segments have voted for this position? -\> CMPS number

<!-- end list -->

6.  if we consider all segments

<!-- end list -->

  - when comparing x and y, a KM comparison

<!-- end list -->

``` r
extract_feature_cmps(x, y, seg_length = 50, seg_scale_max = 3, Tx = 25, 
                     npeaks.set = c(5, 3, 1), full_result = TRUE)
#> $CMPS.score
#> [1] 14
#> 
#> $rec.position
#> [1] -5
#> 
#> $pos.df
#>   position cmps
#> 1       -9   14
#> 2       -8   14
#> 3       -7   14
#> 4       -6   14
#> 5       -5   14
#> 6       -4   14
#> 7       -3   14
#> 8       -2   14
#> 9       -1   14
#> 
#> $nseg
#> [1] 22
```

  - the result of a known non match (KNM) comparison:

<!-- end list -->

``` r
land23 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
land13 <- bullets$sigs[bullets$bulletland == "1-3"][[1]]

result <- extract_feature_cmps(land23$sig, land13$sig, seg_length = 50, seg_scale_max = 3, Tx = 25, 
                     npeaks.set = c(5, 3, 1), full_result = TRUE)
result$CMPS.score
#> [1] 1
```

## Full Comparison Between Two Bullets

`extract_feature_cmps()` can also be used in a pipeline fashion

``` r
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

cp1 <- comparisons %>% select(land1, land2, cmps_score, cmps_nseg)
cp1
```

## Future Work

  - improve and manage to submit the CMPS package to CRAN

  - cross-validate some hyper-parameters

  - apply CMPS to
    
      - the dataset used in the paper to see if we can reproduce their
        results
      - other datasets
