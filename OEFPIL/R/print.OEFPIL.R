#' @name print.OEFPIL
#' @title Print function for an object of class 'OEFPIL'
#' @description Function prints the information about an object of class \code{"OEFPIL"}.
#'
#' @param x an object of class \code{"OEFPIL"} (a result of a call to \code{\link{OEFPIL}}).
#' @param ... other arguments.
#'
#' @details Function prints the formula and estimated parameters of the model from an \code{OEFPIL} object.
#'
#' @seealso \code{\link{OEFPIL}}
#'
#' @examples
#' \dontshow{
#' utils::example("coef.OEFPIL",echo=FALSE)}
#' ##-- Continuing the coef.OEFPIL(.) example:
#'
#' ##Use of print function
#' print(st1)
#'
#' @export


print.OEFPIL <- function(x,...) {
  ## A print method for "OEFPIL".

  l <- (length(x) - 8) / 3
  ## number of parameters

  cat("OEFPIL\n\n")
  cat("Formula:\n")
  cat(x$contents$input.form.string)
  cat("\n\n")
  cat("Estimated parameters:\n")
  cat(unlist(x[1:l]))

  invisible(x)
}
