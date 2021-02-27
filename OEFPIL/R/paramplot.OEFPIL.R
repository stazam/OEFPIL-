#' @name paramplot.OEFPIL
#' @title Parameter estimation graph with ggplot of object of class 'OEFPIL'
#' @description Function for parameter estimation graph with error bars (plus minus standard deviation) with ggplot.
#'
#' @usage paramplot.OEFPIL(object)
#'
#' @param object object or \code{list} of class 'OEFPIL'
#'
#' @details In the case, if we add estimated parameters with NaN, NA, Inf or -Inf, we get empty plot.
#'
#' @note Due to possible large differences in units of estimated parameters, the \code{scale} argument for facetting is set to "free". It should be taken into account when interpreting the results.
#'
#' @seealso \code{\link{OEFPIL}}, \code{\link{curvplot.OEFPIL}} and code{\link{plot.OEFPIL}}.
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
#' ##Use of paramplot.OEFPIL function on object of class 'OEFPIL'
#' paramplot.OEFPIL(st2)
#'
#' ##Use of paramplot.OEFPIL function on list of object of class 'OEFPIL'
#' paramplot.OEFPIL(list(st1,st2))
#'
#' ##Now we look at the case when we plot parameters with Nan values (OEPFIL function does not converge).
#' startsteam <- list(b1 = 0.1, b2 = 5, b3 = 200)
#' st3 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 0.1, b2 = 5, b3 = 200),
#'               CM1, useNLS = F)
#' paramplot.OEFPIL(st3)
#'
#' @import ggplot2
#' @export



paramplot.OEFPIL <- function(output.list){
  ## Function for parameter estimation graph with erorr bars (plus minus standard deviation)
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
