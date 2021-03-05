#' @name curvplot.OEFPIL
#' @title Ggplot function of estimated curve of object of class 'OEFPIL'
#' @description Function for plotting estimated curve with pointwise confidence bands for 'OEFPIL' object.
#'
#' @usage curvplot.OEFPIL(object, signif.level, xx)
#'
#' @param object object of class 'OEFPIL'
#' @param signif.level numeric value or vector of significance levels for confidence interval
#' @param print print out result summaries in the console (default \code{TRUE})
#'
#' @details If the \code{signif.level} parameter is missing, the estimated curve will be plotted without confidence bands.
#'          The rest of parameters can be added similarly as in ggplot.
#'
#' @seealso \code{\link{OEFPIL}}, \code{\link{paramplot.OEFPIL}} and code{\link{plot.OEFPIL}}.
#'
#' @import ggplot2
#' @examples
#' library(MASS)
#' library(ggplot2)
#'
#' ##Creating a data file
#' steamdata <- steam
#' colnames(steamdata) <- c("x","y")
#' n <- nrow(steamdata)
#' CM1 <- diag(rep(10,2*n))
#' CM2 <- diag(c(rep(12,n), rep(14,n)))
#'
#' ##Creating OEFPIL objects
#' st1 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 5, b2 = 8, b3 = 200),
#'              CM1, useNLS = F)
#' st2 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 5, b2 = 8, b3 = 200),
#'              CM2, useNLS = F)
#'
#' ##Use of curvplot.OEFPIL function on object of class 'OEFPIL'
#' curvplot.OEFPIL(st1, signif.level = 0.05)
#'
#' ##Use of curvplot.OEFPIL function on object of class 'OEFPIL' with different arguments
#' curvplot.OEFPIL(st2, signif.level = c(0.01,0.05), xx = seq(0,110,1))
#'
#' ##Use of curvplot.OEFPIL function with additional arguments as for ggplot2
#' curvplot.OEFPIL(st1, signif.level = 0.05) + labs(x = "New x label") + labs(title = "New graph title")
#'
#' @export



curvplot.OEFPIL <- function(output.form, signif.level, xx){
  ## Function for plotting estimated curve with pointwise confidence bands
  ## Input: output.form ... 'OEFPIL' object
  ##        signif.level ... vector of significance levels
  ##        (if missing, only estimated curve is plotted)
  ##        xx ... sequence of x-coordinates of points for computing and plotting confidence bands
  ##        (if missing, the default sequence seq(from = min(x), to = max(x), length.out = 301) is used)

  if(!is.list(output.form)){
    stop("Input has to be an 'OEFPIL' object.")
  } #input control

  out.form <- output.form[names(output.form) != "logs" & names(output.form) != "cov.m_nlsLM"]

  if(!IsListOK(out.form)){
    stop("There are NA or NaN in estimated parameter values or in the covariance matrix.")
  } #control of NaN values

  x <- output.form$contents[[3]]
  y <- output.form$contents[[4]]

  dep.var.name <- output.form$contents$dep.var.name ## name of dependent variable
  idp.var.name <- output.form$contents$idp.var.name ## name of independent variable

  if (missing(xx)) {
    xx <- seq(from = min(x), to = max(x), length.out = 301)
  } ## default sequence for drawing

  n <- length(xx)
  data <- data.frame(x = x, y = y) ## data frame with original data

  ## graph of data points and estimated curve without conf. bands
  if (missing(signif.level)) {
    LOF <- output.form$contents$LOF ## list of functions
    l <- length(output.form$contents$names.of.parameters)
    lst.parameters <- output.form[1:l] ## odhad parametru
    names(lst.parameters) <- output.form$contents$names.of.parameters
    yy <- sapply(xx, function(val, LP){do.call(LOF[[1]], args=c(val, LP))}, lst.parameters)
    datax <- data.frame(x = xx, y = yy)

    ggplot(data) +
      geom_point(aes(x,y),pch = 1, cex = 2) +
      labs(y = dep.var.name, x = idp.var.name) +
      geom_line(data = datax, aes(x,y), colour = "blue", size = 1)
  } else{

    ## creating data frame for the cofidence bands plotting
    CB <- confBands.OEFPIL(output.form, xx = xx, signif.level = signif.level) ## computing CB
    d <- length(signif.level) ## number of different confidence bands
    labels.vec <- paste((1-signif.level) * 100, "% CB", sep = "") ## name of CB for legend

    int <- c()
    for (i in 1:d){
      int[(2*n*(i-1) +1):(2*i*n)]<- c(CB$PointwiseCB[,i],rev(CB$PointwiseCB[,2*d - i+1]))
    }
    siglevel <- c()
    for(i in 1:d){
      siglevel[((i-1)*n*2 +1):(2*n*i)] <- rep(signif.level[i],2*n)
    }
    dataxx <- data.frame(x = rep(c(xx,rev(xx)), d), y = rep(c(CB$yy, rev(CB$yy)),d),
                         int = int, signif.level = as.factor(siglevel))

    ## graph of estimated curve with pointwise confidence bands
    ggplot(data) +
      geom_point(aes(x,y),pch = 1, cex = 2) +
      labs(y = dep.var.name, x = idp.var.name) +
      geom_polygon(data = dataxx, aes(x = x, y = int, fill = signif.level),
                   alpha = 0.3) +
      geom_line(data = dataxx, aes(x,y), colour = "blue", size = 1) +
      scale_fill_discrete(name = "Pointwise CB", labels = labels.vec)
  }
}
