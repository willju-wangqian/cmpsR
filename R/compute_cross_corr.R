#' @useDynLib CMPS COMPUTE_CROSS_CORR_
compute_cross_corr <- function(x, y, min.overlap) .Call(COMPUTE_CROSS_CORR_, x, y, min.overlap)