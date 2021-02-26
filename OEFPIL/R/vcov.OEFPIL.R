#' @name vcov.OEFPIL
#' @title Extract covariance matrix from an object of class 'OEFPIL'
#' @description Function for extracting the estimated covariance matrix from 'OEFPIL' object.
#' @usage ## S3 method for class 'OEFPIL'
#'    vcov(object)
#'
#' @param object object of class 'OEFPIL'
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
#' ##Use of vcov function
#' vcov(st1)
#'
#' @export




vcov.OEFPIL <- function(output.form) {
  ## Function for extracting the estimated covariance matrix from 'OEFPIL' object.

  return(output.form$cov.m_Est)
}
