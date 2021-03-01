#' @name paramplot.OEFPIL
#' @title Plot parameters of an OEFPIL object
#' @description Function for plotting the estimated values of the parameters with error bars (plus minus standard deviation) using \code{ggplot} for an object (or list of objects) of class \code{"OEFPIL"}.
#'
#' @usage paramplot.OEFPIL(object)
#'
#' @param object an object or a \code{list} of objects of class \code{"OEFPIL"} (a result of a call to \code{\link{OEFPIL}}).
#'
#' @details The input list has to be without \code{NaN}, \code{NA}, \code{Inf} or \code{-Inf} values in the estimated parameters or covariance matrix in the source \code{"OEFPIL"} object. In that case the function returns a warning message and no graph is plotted (see Example 3).
#'
#' @note Due to possible large differences in units of estimated parameters, the \code{scale} argument for facetting in the \code{ggplot} graph is set to \code{"free"}. It should be taken into account when interpreting the results.
#'
#' @seealso \code{\link{OEFPIL}}, \code{\link{curvplot.OEFPIL}} and \code{\link{plot.OEFPIL}}.
#'
#' @examples
#' \dontshow{
#' utils::example("coef.OEFPIL",echo=FALSE)}
#' ##-- Continuing the coef.OEFPIL(.) example:
#'
#' CM2 <- diag(c(rep(0.2^2,n), rep(0.1^2,n)))
#' st2 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 5, b2 = 8, b3 = 200),
#'              CM2, useNLS = F)
#'
#' ##Example 1 - Use of paramplot.OEFPIL function on an object of class 'OEFPIL'
#' paramplot.OEFPIL(st2)
#'
#' ##Example 2 - Use of paramplot.OEFPIL function on a list of objects of class 'OEFPIL'
#' paramplot.OEFPIL(list(st1,st2))
#'
#' ##Example 3 - Use of paramplot.OEFPIL function on an object with NaN values (i. e. OEPFIL function does not converge)
#' startsteam <- list(b1 = 0.1, b2 = 5, b3 = 200)
#' st3 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), startsteam,
#'               CM1, useNLS = F)
#' paramplot.OEFPIL(st3)
#'
#' @import ggplot2
#' @export



paramplot.OEFPIL <- function(output.list){
  ## Function for graph of estimated parameters with error bars (plus minus standard deviation)
  ## with ggplot
  ## output.list . . . list of 'OEFPIL' objects

  if(!is.list(output.list)){
    stop("Input has to be a list.")
  } # control of an input

  out.list <- lapply(output.list, function(i){i[names(i) != "logs" & names(i) != "cov.m_nlsLM"]})

  if(!IsListOK(out.list)){
    stop("There are NA or NaN in estimated parameter values or in the covariance matrix.")
  } # control of NaN values

  # creating data structure for ggplot
  N <- length(output.list)

  if (is.vector(output.list[[1]])){
    # data structure for only one 'OEFPIL' object input
    o.form <- output.list[names(output.list) != "logs" & names(output.list) != "cov.m_nlsLM"]
    if(!IsListOK(o.form)){
      stop("There are NA or NaN in estimated parameter values or in the covariance matrix.")
    }
    coefnames <- output.list$contents$names.of.parameters
    k <- length(coefnames)
    data <- matrix(NA, nrow = k, ncol = 4)
    data <- as.data.frame(data)
    data[,1] <- rep(1, each = k)
    data[,2] <- coefnames
    data[,3] <- unlist(output.list[1:k])
    data[,4] <- sqrt(diag(output.list$cov.m_Est))
  } else {
    # data structure for list of 2 or more 'OEFPIL' objects
    coefnames <- output.list[[1]]$contents$names.of.parameters
    k <- length(coefnames)
    data <- matrix(NA, nrow = k*N, ncol = 4)
    data <- as.data.frame(data)
    data[,1] <- rep(1:N, each = k)
    data[,2] <- rep(coefnames,N)
    for (i in 1:N){
      l <- k*i -k +1
      data[l:(k*i),4] <- sqrt(diag(output.list[[i]]$cov.m_Est))
      data[l:(k*i),3] <- unlist(output.list[[i]][1:k])
    }
  }

  names(data[,3]) <- NULL
  colnames(data) <- c("model","cf", "est", "sdest")
  data$model <- as.factor(data$model)

  # plotting graph with ggplot
  ggplot(data, aes(x = model, y = est, col = model)) +
    geom_pointrange(aes(ymin = est - sdest, ymax =  est + sdest))+
    labs(x = "", y = "") +
    facet_wrap(~ cf, scale = "free") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
}
