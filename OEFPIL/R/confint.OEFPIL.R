#' @name confint.OEFPIL
#' @title Calculate confidence intervals for an object of class 'OEFPIL'
#' @description Function calculate confidence intervals for parameters counted by OEFPIL function.
#' @usage ## S3 method for class 'OEFPIL'
#'    confint(object, signif.level = output.form$contents$signif.level)
#'
#' @param object object of class 'OEFPIL'.
#'
#' @param signif.level numerical value or vector of significance levels for confidence intervals.
#'
#' @details We can add one numerical value or vector of numerical values of significance levels for confidence intervals. The default case for parameter \code{signif.level} is value from \code{OEFPIL} object.
#'
#' @return Matrix of estimated confidence intervals for model coefficients from on object \code{OEFPIL}
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
#' ##Use of confint function
#' #one numerical value
#' confint(st1)
#'
#' #vector of numerical values
#' confint(st1,signif.level = c(0.01,0.05,0.1))
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
