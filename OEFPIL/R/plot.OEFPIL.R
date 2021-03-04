#' @name plot.OEFPIL
#' @title Plot the estimate from an OEPFIL object
#' @description Plot of the iterated linearization estimate of a function from an \code{"OEFPIL"} object with pointwise confidence bands.
#' @usage ## S3 method for class 'OEFPIL'
#'     plot(object, xx, signif.level,...)
#'
#' @param object an object of class \code{"OEFPIL"} (a result of a call to \code{\link{OEFPIL}}).
#' @param xx  a sequence of x-coordinates of points for computing and plotting confidence bands. If missing, the default sequence \code{seq(from = min(x), to = max(x), length.out = 301)} is used.
#' @param signif.level a numerical value or a vector of significance levels for confidence bands.
#' @param ... additional arguments (same as in \link{plot} function) affecting the plot.
#' @details If the \code{signif.level} argument is missing, the value is set to 0.05, but the confidence bands are not plotted.
#'          The confidence bands are computing under normality assumption.
#'
#'@return Returns an object of type list containing at least the following components
#'  \itemize{ \item \code{xx} a numerical vector of points where bands are calculated.
#'            \item \code{yy} a numerical vector with values of estimated function in \code{xx}.
#'            \item \code{PointwiseCB} a matrix of pointwise confidence bands at points \code{xx}.
#'            }
#'
#' @seealso \code{\link{OEFPIL}}
#'
#' @examples
#' \dontshow{
#' utils::example("coef.OEFPIL",echo=FALSE)}
#' ##-- Continuing the coef.OEFPIL(.) example:
#'
#' ##Use of plot function with default parameters, signif.level is not set up...No confidence bands are plotted
#' plot(st1)
#'
#' ##Use of plot function with different parameters
#' plot(st1, seq(0,113,0.1), signif.level = c(0.01,0.05), main = "Graph of estimated function")
#'
#' ##Return values of plot function
#' (a <- plot(st1, signif.level = 0.05))
#'
#'
#' @export



plot.OEFPIL <- function(output.form, xx, signif.level, ...) {
  ## Function plots confidence bands of list from OEFPIL() function.
  ## output.form . . . output from OEFPIL()
  ## xx          . . . in these points we calculate and plot CI (confidence intervals) or
  ##                    CB (conf. bands)

  draw.CB <- T

  if (missing(signif.level)) {
    signif.level <- 0.05
    draw.CB <- F
  }
  d <- length(signif.level)

  x <- output.form$contents[[3]] ## x data
  y <- output.form$contents[[4]] ## y data

  dep.var.name <- output.form$contents$dep.var.name ## name of dependant variabe
  idp.var.name <- output.form$contents$idp.var.name ## name of independant variable

  if (missing(xx)) {
    xx <- seq(from = min(x), to = max(x), length.out = 301)
  }

  CB <- confBands.OEFPIL(output.form, xx = xx, signif.level = signif.level)

  plot(x, y, xlab = idp.var.name, ylab = dep.var.name, ... = ...)

  lines(CB$xx, CB$yy, lwd = 2, col = "black")

  if (draw.CB == T) {
    for (i in 0:(d-1)) {
      lines(CB$xx, CB$PointwiseCB[, i+1], lwd = 2, lty = 2, col = i + 2)
      lines(CB$xx, CB$PointwiseCB[, 2*d - i], lwd = 2, lty = 2, col = i + 2)

      lines(CB$xx, CB$PredictCB[, i+1], lwd = 2, lty = 3, col = i + 2)
      lines(CB$xx, CB$PredictCB[, 2*d - i], lwd = 2, lty = 3, col = i + 2)
    }
  }

  return(invisible(list(xx = CB$xx, yy = CB$yy, PointwiseCB = CB$PointwiseCB,
                        PredictCB = CB$PredictCB)))
}
