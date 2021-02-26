#' @name coef.OEFPIL
#' @title Extract model coefficients from an object of class 'OEFPIL'
#' @description Function which extracts the estimated model coefficients from 'OEFPIL' object.
#' @usage ## S3 method for class 'OEFPIL'
#'    coef(object)
#'
#' @param object object of class 'OEFPIL'
#'
#' @return Vector of estimated model coefficients from on object \code{OEFPIL}
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

