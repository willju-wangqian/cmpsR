#' Function to calculate the cross-correlation between two sequences
#' 
#' some more description, x is the longer (reference) sequence
#' @param x numeric sequence of values
#' @param y numeric sequence of values
#' @param min.overlap integer, minimal number of values in the overlap between sequences x and y to calculate a correlation value. Set to 10 percent of the maximum length of either sequence (HH: this might be problematic for CMPS)
#' @export
#' @return list consisting of the lag where the maximum correlation is achieved, and the maximum correlation value.
get_ccf3 <- function (x, y, min.overlap = round(0.1 * max(length(x), length(y)))) 
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
    cor(xx, lag(yy, lag), use = "pairwise.complete")
  })
  ns <- sapply(lags, function(lag) {
    dim(na.omit(cbind(xx, lag(yy, lag))))[1]
  })
  
  cors[ns < min.overlap] <- NA
  lag <- lags - (ny - min.overlap)
  return(list(lag = lag, ccf = cors))
}
