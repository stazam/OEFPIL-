#' @name confBands.OEFPIL
#' @title Calculate confidence bands for an object of class 'OEFPIL'
#' @description Function calculates pointwise confidence bands of estimated function from 'OEFPIL'.
#' @usage ## S3 method for class 'OEFPIL'
#'    confBands(object, xx, signif.level = 0.05)
#'
#' @param object object of class 'OEFPIL'.
#' @param xx numerical vector of points, where confidence bands will be calculated.
#' @param signif.level numerical value or vector of significance levels for confidence intervals (default value is 0.05)
#'
#' @details We can add one numerical value or vector of numerical values of significance levels for confidence intervals.
#'
#' @return Matrix with named columns of estimated pointwise confidence bands of estimated function from on object \code{OEFPIL}. And also points where the bands are calculated.
#'
#' @seealso \code{\link{OEFPIL}}
#'
#' @examples
#' \dontshow{
#' utils::example("coef.OEFPIL",echo=FALSE)}
#' ##-- Continuing the coef.OEFPIL(.) example:
#'
#' ##Use of confBands function with default parameters
#' (a <- confBands(st1))
#'
#' #vector of numerical values
#' (b <- confBands(st1,signif.level = c(0.01,0.05)))
#'
#' @export
confBands <- function(x, xx, signif.level) {
  UseMethod("confBands")
}

#' @export
confBands.OEFPIL <- function(output.form, xx, signif.level = 0.05, new.obs.variance) {
  ## This is for calculating confidence bands of estimated function from OEFPIL.
  ## output.form      . . . output from OEFPIL()
  ## xx               . . . in these points we calculate CI (confidence intervals) or
  ##                        CB (conf. bands)
  ## new.obs.variance . . . a variance of the new observation;
  ##                        it is needed for prediction intervals.

  LOF <- output.form$contents$LOF ## list of functions
  x <- output.form$contents[[3]] ## x data
  y <- output.form$contents[[4]] ## y data
  CM <- output.form$contents$CM ## covariance matrix of the data

  cov_m <- output.form$cov.m_Est ## estimate of covariance matrix of parameters
  l <- length(output.form$contents$names.of.parameters) ## number of parameters

  lst.parameters <- output.form[1:l] ## parameter estimation
  names(lst.parameters) <- output.form$contents$names.of.parameters

  lst.parameters_previous.step <- output.form[(2*l+6):(3*l+5)]
  names(lst.parameters_previous.step) <- output.form$contents$names.of.parameters
  ## estimate from the previous step

  if (IsListOK(lst.parameters) && IsListOK(lst.parameters_previous.step) && IsListOK(cov_m)) {

    if (missing(xx)) {
      xx <- seq(from = min(x), to = max(x), by = 0.1)
    }
    yy <- sapply(xx, function(val, LP){do.call(LOF[[1]], args=c(val, LP))}, lst.parameters)

    Omega <- sapply(1:l, function(i) {
      sapply(xx, function(val, LP){do.call(LOF[[2+i]], args=c(val, LP))}, lst.parameters_previous.step)
    })
    ## i-th row of matrix is value of vector omega in the point xx[i]

    variance <- apply(((Omega %*% cov_m) * Omega), 1, sum)

    if (missing(new.obs.variance)) {

      n <- length(diag(CM)) / 2
      new.obs.variance <- mean(diag(CM)[1:n]) + mean(diag(CM)[(n+1):(2*n)])
      ## "estimation" of variance of the new observation (needed for prediction interval)

    }


    sl <- sort(c(signif.level/2, 1 - signif.level/2), decreasing = F)

    d <- length(signif.level)
    k <- length(xx)

    PCB_lwr <- matrix(rep(yy, d), k, d) + matrix(rep(qnorm(sl[1:d]), k), k, d, byrow = T) * sqrt(variance)
    PCB_upr <- matrix(rep(yy, d), k, d) + matrix(rep(qnorm(sl[(d+1):(2*d)]), k), k, d, byrow = T) * sqrt(variance)
    ## pointwise confidence band

    PredictCB_lwr <- matrix(rep(yy, d), k, d) + matrix(rep(qnorm(sl[1:d]), k), k, d, byrow = T) * sqrt(variance + new.obs.variance)
    PredictCB_upr <- matrix(rep(yy, d), k, d) + matrix(rep(qnorm(sl[(d+1):(2*d)]), k), k, d, byrow = T) * sqrt(variance + new.obs.variance)
    ## pointwise confidence band

    PointwiseCB <- cbind(PCB_lwr, PCB_upr)
    PredictCB <- cbind(PredictCB_lwr, PredictCB_upr)

    colnames(PointwiseCB) <- paste(round(sl * 100, 2), "%")
    colnames(PredictCB) <- paste(round(sl * 100, 2), "%")

    return(invisible(list(xx = xx, yy = yy, PointwiseCB = PointwiseCB, PredictCB = PredictCB)))
  }
}



