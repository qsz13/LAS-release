# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

stl_sort <- function(x) {
    .Call('LAS_stl_sort', PACKAGE = 'LAS', x)
}

normalizeInput <- function(x) {
    invisible(.Call('LAS_normalizeInput', PACKAGE = 'LAS', x))
}

#' @export
normalizeInputMatrix <- function(x) {
    .Call('LAS_normalizeInputMatrix', PACKAGE = 'LAS', x)
}

