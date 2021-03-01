#' @name coef.OEFPIL
#' @title Extract model coefficients from an 'OEFPIL' object
#' @description Function which extracts the estimated model coefficients from an object of class 'OEFPIL'.
#' @usage ## S3 method for class 'OEFPIL'
#'    coef(object)
#'
#' @param object an object of class 'OEFPIL' (a result of a call to \code{OEFPIL})
#'
#' @return A named vector of estimated model coefficients extracted from an 'OEFPIL' object.
#'
#' @seealso \code{\link{OEFPIL}}
#'
#' @examples
#'
#' ##Creating a data file (using steam data from MASS library)
#' library(MASS)
#' steamdata <- steam
#' colnames(steamdata) <- c("x","y")
#' startsteam <- list(b1 = 5, b2 = 8, b3 = 200)
#' k <- nrow(steamdata)
#' CM <- diag(rep(0.1,2*k))
#'
#' ##Creating an OEFPIL object
#' st1 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), startsteam, CM, useNLS = F)
#'
#' ##Use of coef function
#' coef(st1)
#'
#' @export

coef.OEFPIL <- function(output.form) {
  ## Function for extracting the estimated coefficients from OEFPIL object.

  l <- (length(output.form) - 8) / 3
  ## number of parameters

  output <- unlist(output.form[1:l])
  names(output) <- output.form$contents$names.of.parameters
  return(output)
}

