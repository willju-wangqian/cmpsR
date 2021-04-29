#' @useDynLib CMPS NA_TRIM_
na_trim_cmps <- function(x) .Call(NA_TRIM_, x)