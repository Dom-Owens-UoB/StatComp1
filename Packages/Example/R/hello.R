# An example package for Statistical Computing 1
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Title
#' Example package
#' @param x
#'
#' @return
#' SUms entries of vector
#' @export
#'sum
#' @examples
sum <- function(x) {
  s <-0
  for (i in 1:length(x)) {
    s <- s + x[i]
  }
  s
}
