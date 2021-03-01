#' @name confint.OEFPIL
#' @title Confidence intervals for OEFPIL parameters
#' @description Function computes confidence intervals for the parameters counted by \code{OEFPIL} function.
#' @usage ## S3 method for class 'OEFPIL'
#'    confint(object, signif.level = output.form$contents$signif.level)
#'
#' @param object an object of class \code{"OEFPIL"} (a result of a call to \code{\link{OEFPIL}}).
#'
#' @param signif.level a numerical value or a vector of significance levels for confidence intervals. If missing, a value from the input \code{"OEFPIL"} object is used.
#'
#' @details The confidence intervals are computing under normality assumption.
#'
#' @return A matrix of estimated confidence intervals for model coefficients from an \code{"OEFPIL"} object. The matrix contains lower and upper confidence limits (columns) for each parameter (rows).
#'
#' @seealso \code{\link{OEFPIL}}
#'
#' @examples
#' \dontshow{
#' utils::example("coef.OEFPIL",echo=FALSE)}
#' ##-- Continuing the coef.OEFPIL(.) example:
#'
#' ##Use of confint function
#' #one numerical value
#' confint(st1)
#'
#' #vector of numerical values
#' confint(st1, signif.level = c(0.01,0.05,0.1))
#'
#' @export


confint.OEFPIL <- function(output.form, signif.level = output.form$contents$signif.level) {
  ## Function calculate confidence intervals for parameters counted by OEFPIL function.

  if (!( is.vector(signif.level) && is.numeric(signif.level))){
    stop("Input of significance levels is not a numerical vector" )
  } #check if the input is numerical vector

  if ( !( all( signif.level <  1) & all( signif.level > 0) ) ){
    stop("Values for significance level should be between 0 and 1!")
  } #check if the values are between zero and one

  cov_m <- output.form$cov.m_Est ## Estimate of covariance matrix
  l <- dim(cov_m)[1] ## number of parameters

  lst.parameters <- output.form[1:l]

  if (IsListOK(lst.parameters) && IsListOK(cov_m)) {

    d <- length(signif.level)

    sl <- sort(c(signif.level/2, 1 - signif.level/2), decreasing = F)

    vec.parameters <- unlist(lst.parameters)

    CI.matrix <- cbind(vec.parameters + matrix(rep(qnorm(sl[1:d]), l), l, d, byrow = T) * sqrt(diag(cov_m)),
                       vec.parameters + matrix(rep(qnorm(sl[(d+1):(2*d)]), l), l, d, byrow = T) * sqrt(diag(cov_m)))

    row.names(CI.matrix) <- names(vec.parameters)
    colnames(CI.matrix) <- paste(round(sl * 100, 2), "%")

    return(CI.matrix)

  } else {

    return(NA)

  }

}
