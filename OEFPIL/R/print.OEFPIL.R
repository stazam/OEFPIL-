#' @name print.OEFPIL
#' @title Print function for an object of class 'OEFPIL'
#' @description Function will print out the information about object of class 'OEFPIL'.
#' @usage ## S3 method for class 'OEFPIL'
#'    print(object)
#'
#' @param object object of class 'OEFPIL'.
#'
#' @details Function will print the formula and estimated parameters of model from an object \code{OEFPIL}.
#'
#' @seealso \code{\link{OEFPIL}}
#'
#' @examples
#' library(MASS)
#'
#' ##Creating a data file
#' steamdata <- steam
#' colnames(steamdata) <- c("x","y")
#' startsteam <- list(b1 = 5, b2 = 8, b3 = 200)
#' k <- nrow(steamdata)
#' CM <- diag(rep(0.1,2*k))
#'
#' ##Creating OEFPIL object
#' st1 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), startsteam, CM, useNLS = F)
#'
#' ##Use of print function
#' print(st1)
#'
#' @export


print.OEFPIL <- function(output.form) {
  ## A print method for "OEFPIL".

  l <- (length(output.form) - 8) / 3
  ## number of parameters

  cat("OEFPIL\n\n")
  cat("Formula:\n")
  cat(output.form$contents$input.form.string)
  cat("\n\n")
  cat("Estimated parameters:\n")
  cat(unlist(output.form[1:l]))

  invisible(output.form)
}
