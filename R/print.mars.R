source("R/mars.r")

#' Print a mars object
#'
#' @param object A mars object
#' @param ... Additional arguments (not used)
#' @export
print.mars <- function(object, ...) {
  cat("Call:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  coef_values <- object$coefficients
  print(coef_values)
}

