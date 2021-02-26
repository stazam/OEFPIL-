#' @name summary.OEFPIL
#' @title Summary function from an object of class 'OEFPIL'
#' @description Function for fast and clean output of all basic information of 'OEFPIL' object.
#' @usage ## S3 method for class 'OEFPIL'
#'    summary(object, signif.level = output.form$contents$signif.level, print = T)
#'
#' @param object object of class 'OEFPIL'
#' @param signif.level significance level for confidence interval
#' @param print print out result summaries in the console (default \code{TRUE})
#'
#' @details The default case for parameter \code{signif.level} is value from \code{OEFPIL} object.
#'
#' @return Returns an object of type list containing at least the following components
#'  \itemize{\item \code{param_Est} is the (numerical) vector of estimated model parameters
#'            \item  \code{sd} standard deviations for estimated model parameters.
#'            \item \code{cov.m_Est} covariance matrix of estimated model parameters.
#'            \item \code{it_num} number of iterations
#'            \item \code{CI_parameters} matrix of lower and upper bounds for confidence intervals.
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
#' CM <- diag(rep(0.1,2*k))
#'
#' ##Creating OEFPIL object
#' st1 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), startsteam, CM, useNLS = F)
#'
#' ##Use of summary function with default parameters
#' summary(st1)
#'
#' ##Use of summary function with different parameters
#' summary(st1, signif.level = 0.05, print = F)
#'
#' @export


summary.OEFPIL <- function(output.form, signif.level = output.form$contents$signif.level,
                           print = T) {
  ## Function for fast and clean output of all basic information
  ## output.form  . . . output from OEFPIL
  ## signif.level . . . significance level
  ## print        . . . states, that if we want to output table as well


  if (IsListOK(output.form$cov.m_Est)) {

    l <- dim(output.form$cov.m_Est)[1] ## number of parameters

    summary.form <- output.form[c(1:l, 2*l+1, 2*l+5, 3*l+6)]
    ## summary.form is portion of list, which contains just components required
    ## for summary function
    summary.form$input.form.string <- output.form$contents$input.form.string
    ## choice of inevitable parameters from list

    if (IsListOK(summary.form)) {

      l <- dim(summary.form$cov.m_Est)[1] ## number of parameters

      if (length(signif.level) != length(output.form$contents$signif.level)) {
        summary.form$CI_parameters <- confint.OEFPIL(output.form, signif.level)
        ## We are calculating new CI, in case of the confidence level was assigned

      } else if (all(signif.level != output.form$contents$signif.level)) {
        summary.form$CI_parameters <- confint.OEFPIL(output.form, signif.level)
        ## We are calculating new CI, in case of the confidence level was assigned

      }

      d <- dim(summary.form$CI_parameters)[2]

      Param.mat <- matrix(0, nrow = l, ncol = 2 + d)
      Param.mat[,1] <- unlist(summary.form[1:l]) ## parameters
      Param.mat[,2] <- sqrt(diag(summary.form$cov.m_Est)) ## sd
      Param.mat[,3:(2 + d)] <- summary.form$CI_parameters ## CI

      rownames(Param.mat) <- rownames(summary.form$cov.m_Est)
      colnames(Param.mat) <- c("Param Est", "        Std Dev",
                               paste("  CI Bound", colnames(summary.form$CI_parameters)))

      if (print == T) {
        cat(paste(c("Summary of the result: ", "\n", "\n"), sep = ""))
        cat(paste(summary.form$input.form.string, "\n\n", sep = ""))
        print(Param.mat)
        cat(paste(c("\n", "Estimated covariance matrix:", "\n"), sep = ""))
        print(summary.form$cov.m_Est)
        cat(paste(c("\n", "Number of iterations:", summary.form$it_num, "\n"), sep = ""))
      }

      output <- list(param_Est = Param.mat[,1], sd =  Param.mat[,2],
                     cov.m_Est = summary.form$cov.m_Est, it_num = summary.form$it_num,
                     CI_parameters = Param.mat[,3:(2 + d)])

      return(invisible(output))

    } else {
      logg <- paste("The summary cannot be calculated because of some NaN, NA, Inf or -Inf values",
                    "\n", "in the OEFPIL object. Logs from OEFPIL:", sep = "")
      message(logg)
      message(output.form$logs)
    }
  }
}
