
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

cmps <- extract_feature_cmps(x, y, full_result = T)
cmps
#> $CMPS.score
#> [1] 0
#> 
#> $rec.position
#> NULL
#> 
#> $pos.df
#> NULL
#> 
#> $nseg
#> [1] 22
```

    #> - Session info ---------------------------------------------------------------
    #>  setting  value                       
    #>  version  R version 3.6.3 (2020-02-29)
    #>  os       Windows 10 x64              
    #>  system   x86_64, mingw32             
    #>  ui       RTerm                       
    #>  language (EN)                        
    #>  collate  English_United States.1252  
    #>  ctype    English_United States.1252  
    #>  tz       America/Chicago             
    #>  date     2020-07-09                  
    #> 
    #> - Packages -------------------------------------------------------------------
    #>  package     * version    date       lib source        
    #>  assertthat    0.2.1      2019-03-21 [1] CRAN (R 3.6.2)
    #>  backports     1.1.5      2019-10-02 [1] CRAN (R 3.6.1)
    #>  callr         3.4.3      2020-03-28 [1] CRAN (R 3.6.3)
    #>  cli           2.0.2      2020-02-28 [1] CRAN (R 3.6.3)
    #>  CMPS        * 0.0.0.9000 2020-07-09 [1] local         
    #>  crayon        1.3.4      2017-09-16 [1] CRAN (R 3.6.2)
    #>  desc          1.2.0      2018-05-01 [1] CRAN (R 3.6.3)
    #>  devtools      2.3.0      2020-04-10 [1] CRAN (R 3.6.3)
    #>  digest        0.6.25     2020-02-23 [1] CRAN (R 3.6.3)
    #>  ellipsis      0.3.0      2019-09-20 [1] CRAN (R 3.6.2)
    #>  evaluate      0.14       2019-05-28 [1] CRAN (R 3.6.2)
    #>  fansi         0.4.1      2020-01-08 [1] CRAN (R 3.6.2)
    #>  fs            1.4.0      2020-03-31 [1] CRAN (R 3.6.3)
    #>  glue          1.4.1      2020-05-13 [1] CRAN (R 3.6.3)
    #>  htmltools     0.4.0      2019-10-04 [1] CRAN (R 3.6.2)
    #>  knitr         1.28       2020-02-06 [1] CRAN (R 3.6.3)
    #>  lattice       0.20-38    2018-11-04 [2] CRAN (R 3.6.3)
    #>  magrittr      1.5        2014-11-22 [1] CRAN (R 3.6.2)
    #>  memoise       1.1.0      2017-04-21 [1] CRAN (R 3.6.3)
    #>  pkgbuild      1.0.6      2019-10-09 [1] CRAN (R 3.6.3)
    #>  pkgload       1.0.2      2018-10-29 [1] CRAN (R 3.6.3)
    #>  prettyunits   1.1.1      2020-01-24 [1] CRAN (R 3.6.2)
    #>  processx      3.4.2      2020-02-09 [1] CRAN (R 3.6.3)
    #>  ps            1.3.3      2020-05-08 [1] CRAN (R 3.6.3)
    #>  R6            2.4.1      2019-11-12 [1] CRAN (R 3.6.2)
    #>  Rcpp          1.0.4.6    2020-04-09 [1] CRAN (R 3.6.3)
    #>  remotes       2.1.1      2020-02-15 [1] CRAN (R 3.6.3)
    #>  rlang         0.4.6      2020-05-02 [1] CRAN (R 3.6.3)
    #>  rmarkdown     2.1        2020-01-20 [1] CRAN (R 3.6.2)
    #>  rprojroot     1.3-2      2018-01-03 [1] CRAN (R 3.6.3)
    #>  sessioninfo   1.1.1      2018-11-05 [1] CRAN (R 3.6.3)
    #>  stringi       1.4.6      2020-02-17 [1] CRAN (R 3.6.2)
    #>  stringr       1.4.0      2019-02-10 [1] CRAN (R 3.6.2)
    #>  testthat      2.3.2      2020-03-02 [1] CRAN (R 3.6.3)
    #>  usethis       1.6.1      2020-04-29 [1] CRAN (R 3.6.3)
    #>  withr         2.1.2      2018-03-15 [1] CRAN (R 3.6.2)
    #>  xfun          0.12       2020-01-13 [1] CRAN (R 3.6.2)
    #>  yaml          2.2.1      2020-02-01 [1] CRAN (R 3.6.2)
    #>  zoo           1.8-8      2020-05-02 [1] CRAN (R 3.6.3)
    #> 
    #> [1] C:/Users/juwan/Documents/R/win-library/3.6
    #> [2] C:/Program Files/R/R-3.6.3/library

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
    ccf curve.

<img src="image/step3_2.png" width="80%" />

4.  each ccf curve votes for 5 candidate positions, then we ask two
    questions:

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

  - if we consider all segments when comparing x and y, a KM comparison

<!-- end list -->

``` r
extract_feature_cmps(x, y, seg_length = 50, seg_scale_max = 3, Tx = 25, 
                     npeaks.set = c(5, 3, 1), full_result = T)
#> $CMPS.score
#> [1] 0
#> 
#> $rec.position
#> NULL
#> 
#> $pos.df
#> NULL
#> 
#> $nseg
#> [1] 22
```
