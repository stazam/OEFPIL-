#' @name plot.OEFPIL
#' @title Plot of the iterated linearization estimate from an OEPFIL object
#' @description Plot of the iterated linearization estimate from an OEPFIL object.
#' @usage ## S3 method for class 'OEFPIL'
#'     plot(object, xx, signif.level,...)
#'
#' @param object object of class 'OEFPIL'
#' @param xx numerical vector of points, where confidence bands will be calculated and plotted.
#' @param signif.level numerical value or vector of significance levels for confidence intervals.
#' @param ... additional arguments (same as \link{plot} function) affecting the plot.
#'
#' @details If the signif.level parameter is missing, we do not plot any confidence bands, although the values will be set to 0.05.
#'
#'@return Returns an object of type list containing at least the following components
#'  \itemize{ \item  \code{xx} points where bands are calculated.
#'            \item \code{yy} values of estimated function.
#'            \item \code{PointwiseCB} matrix of pointwise confidence bands at points \code{xx}.
#'            }
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
#' CM <- diag(rep(12,2*k))
#'
#' ##Creating OEFPIL object
#' st1 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), startsteam, CM, useNLS = F)
#'
#' ##Use of plot function with default parameters, signif.level is not set up...No confidence bands
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
  ## signig.level . . . significance level
  ## ...         . . . additional arguments

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
    }
  }

  return(invisible(list(xx = CB$xx, yy = CB$yy, PointwiseCB =CB$PointwiseCB)))
}
