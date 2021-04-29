#' @useDynLib CMPS add_
add <- function(x, y) .Call(add_, x, y)