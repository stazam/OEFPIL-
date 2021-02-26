#' @name ortresiduals.OEFPIL
#' @title Calculation of orthogonal residuals from an object of class 'OEFPIL'
#' @description Function for calculating orthogonal residuals of an 'OEFPIL' object.
#' @usage ortresiduals.OEFPIL(object, a)
#'
#' @param object object of class 'OEFPIL'
#' @param a numeric value, for defining minimization interval for the \code{\link{optimize}} function (if not defined, default value 0.05 * range(x) is used). Must be positive
#'
#' @return Returns an object of type list containing at least the following components
#'  \itemize{\item \code{x.ores} numerical vector of x coordinates of points, where the minimal distance is realized
#'            \item  \code{o.resid} numerical vector of orthogonal residuals (minimal Euclidean distances between data points and estimated function)
#'            \item \code{SSort} orthogonal sum of squares.
#'            }
#'
#' @seealso \code{\link{OEFPIL}}
#'
#' @examples
#'library(MASS)
#'
#' ##Creating a data file
#' steamdata <- steam
#' colnames(steamdata) <- c("x","y")
#' n <- nrow(steamdata)
#' CM1 <- diag(rep(0.1,2*n))
#'
#' ##Creating OEFPIL object
#'st1 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 5, b2 = 8, b3 = 200),
#'              CM1, useNLS = F)
#'
#' ##Use ortresiduals.OEFPIL function on the OEFPIL object, with specified value 'a'
#' ortresiduals.OEFPIL(st1,5)
#'
#' ##Use ortresiduals.OEFPIL function without value 'a' (defaut.value = 0.05 * range(x))
#' ortresiduals.OEFPIL(st1)
#'
#' ##Choice of too narrow interval. Result can be misleading!
#'ortresiduals.OEFPIL(st1,0.5)
#' @export

ortresiduals.OEFPIL <- function(output.form, a){
  ## Function for computing orthogonal residuals of an 'OEFPIL' object
  ## input: output.form ... 'OEFPIL' object
  ##        a ... constant for defining minimization interval for the optimize() function
  ##              (if not defined, default value 0.05 * range(x) is used), a must be positive
  ## output: x.ores ... x coordinates of points, where the minimal distance is realized
  ##         o.resid ... orthogonal residuals (minimal Euclidean distances between data points
  ##                     and estimated function)
  ##         SSort ... orthogonal sum of squares

  if(!is.list(output.form)){
    stop("Input has to be a list.")
  }

  l <- dim(output.form$cov.m_Est)[1] ## number of parameters

  if(!IsListOK(output.form[1:l])){
    stop("There are NA or NaN in estimated parameter values.")
  }

  x <- output.form$contents[[3]] ## x-ova data
  y <- output.form$contents[[4]] ## y-ova data

  ## setting value of argument a (if it is not defined in input)
  if(missing(a)){
    rngx <- max(x) - min(x)
    a <- 0.05 * rngx
    MES <- paste("Argument 'a' was not defined, value ", a, " was used for the calculation.",
                 sep = "")
    message(MES)
  }

  if(a <= 0){
    stop("a must be a positive number")
  }


  for (i in 1:l){
    assign(output.form$contents$names.of.parameters[i],output.form[[i]])
  } ## assigning estimated parameter values

  ## rewriting main function in optimize requested format - ftomin() function
  formstring <- strsplit(output.form$contents$input.form.string, "~")[[1]][2]
  eval(parse(text = paste("ftomin <- function(x){", formstring, "}", sep = "")))

  fceoptim <- function(f, x0, y0, a){
    ## input: f - function to minimize
    ##        x0, y0 - coordinates of the data point
    ##        a - constant for defining minimization interval
    ## output: clasicall output for optimize()

    eudist <- function(x){sqrt((x - x0)^2 + (f(x) - y0)^2)}
    optimize(eudist, c(x0 - a, x0 + a))

  } ## function for finding minimal Euclidean distance from point [x0,y0] to the curve

  ortsumvec <- c()
  xoptvec <- c()
  for(i in 1:length(x)){
    OPT <- fceoptim(ftomin, x[i], y[i], a)
    ortsumvec[i] <- OPT$objective ## vector of minimal Euclidean distances (orthogonal residuals)
    xoptvec[i] <- OPT$minimum ## x coordinates of points, where the minimal distance is realized
  }
  return(list(x.ores = xoptvec, o.resid = ortsumvec, SSort = sum(ortsumvec^2)))
}




