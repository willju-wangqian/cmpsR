#' Function to calculate the cross-correlation between two sequences
#'
#' This function is used for CMPS algorithm. 
#' @param x numeric sequence of values
#' @param y numeric sequence of values
#' @param min.overlap integer, minimal number of values in the overlap between sequences x and y to calculate a correlation value. Set to 10 percent of the maximum length of either sequence (HH: this might be problematic for CMPS)
#'
#' @return list consisting of the lag where the maximum correlation is achieved, and the maximum correlation value.
#' @export
#' @importFrom assertthat assert_that
#' @importFrom stats cor complete.cases
#' @importFrom dplyr lag
#' @examples
#' data("bullets")
#' land2_3 <- bullets$sigs[bullets$bulletland == "2-3"][[1]]
#' land1_2 <- bullets$sigs[bullets$bulletland == "1-2"][[1]]
#' x <- land2_3$sig
#' y <- land1_2$sig
#' 
#' segments <- get_segs(x, len = 50)
#' 
#' ccr <- get_ccf4(y, segments$segs[[7]], 
#'                 min.overlap = length(segments$segs[[7]]))
get_ccf4 <- function (x, y, min.overlap = round(0.1 * max(length(x), length(y)))) 
{
  # requires x to be the longer signature
  x <- as.vector(unlist(x))
  y <- as.vector(unlist(y))
  nx <- length(x)
  ny <- length(y)
  assert_that(is.numeric(x), is.numeric(y))
  assert_that(nx > 0, ny > 0, nx >= ny) # this is the only change
  xx <- c(rep(NA, ny - min.overlap), x, rep(NA, ny - min.overlap))
  yy <- c(y, rep(NA, length(xx) - ny))
  lag.max <- length(yy) - length(y)
  lags <- 0:lag.max
  
  cors <- sapply(lags, function(lag) {
    cor(xx, dplyr::lag(yy, lag), use = "pairwise.complete.obs")
  })
  
  # running time improved with this chunk
  ns <- sapply(lags, function(lag) {
    sum(complete.cases(xx, dplyr::lag(yy, lag)))
  })
  
  cors[ns < min.overlap] <- NA
  lag <- lags - (ny - min.overlap)
  return(list(lag = lag, ccf = cors))
}