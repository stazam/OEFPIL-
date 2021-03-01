#' @name vcov.OEFPIL
#' @title Covariance matrix from an OEFPIL object
#' @description Function for extracting the estimated covariance matrix from an object of class \code{"OEFPIL"}.
#' @usage ## S3 method for class 'OEFPIL'
#'    vcov(object)
#'
#' @param object an object of class \code{"OEFPIL"} (a result of a call to \code{\link{OEFPIL}}).
#'
#' @return A matrix of the estimated covariances between the parameter estimates from an \code{"OEFPIL"} object.
#'
#' @seealso \code{\link{OEFPIL}}
#'
#' @examples
#' \dontshow{
#' utils::example("coef.OEFPIL",echo=FALSE)}
#' ##-- Continuing the coef.OEFPIL(.) example:
#'
#' ##Use of vcov function
#' vcov(st1)
#'
#' @export




vcov.OEFPIL <- function(output.form) {
  ## Function for extracting the estimated covariance matrix from 'OEFPIL' object.

  return(output.form$cov.m_Est)
}
