#' @name OEFPIL
#' @title Estimation of Parameters by Iterated Linearization
#' @description Function for computing optimal estimate of function parameters by iterated linearization.
#' @usage OEFPIL(data, form, start.val, CM = diag(dim(data)[1] * 2),  max.iter = 100,
#'     see.iter.val = F, save.file.name, th = 0.001, signif.level = 0.05, useNLS = T)
#'
#' @param data Data file can be any object of type \code{data.frame} with 2 columns or \code{list} with 2 elements.
#' @param form an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param start.val starting values of estimating parameters.
#' @param CM covariance matrix (does not change in iteration process).
#' @param max.iter maximum number of iterations
#' @param see.iter.val logical value (default \code{FALSE}) indicating if we want to display and save partial results of algorithm.
#' @param save.file.name name of the file for saving results.
#' @param th numerical value, indicating threshold necessary for estimation stoppage.
#' @param signif.level significance level for confidence interval.
#' @param useNLS logical. If \code{TRUE} (the default value), function will set up starting parameters calculated by \code{\link{nlsLM}} function (nonlinear least square estimation).
#'
#' @details Models for OEPFIL function are specified symbolically. A typical model has the form \code{y ~ f(x, a_1,...,a_n)}, where
#'
#'  \itemize{\item \code{y} is the (numerical) response vector
#'            \item  \code{x} is a vector which specifies non-linear predictor  for \code{y}
#'            \item terms \code{a_1,...,a_n} are parameters of specified model. Function \code{f} must have continuous second partial derivatives with respect to \code{x} and parameters \code{a_1,...a_n}.}
#'
#'     In the \code{data} entry of type \code{data.frame}, both columns must be named as variables in formula. The same holds for elements of \code{list}.
#'
#' @return Returns an object of class \code{'OEPFIL'}, which is a list containing at least the following components
#' \itemize{
#'    \item \code{name_est} estimations of model parameters
#'    \item \code{name_upgraded.start.val} estimations of model parameters
#'    \item \code{it_num} number of iterations
#'    \item \code{CI_parameters} list of confidence intervals for estimated parameters (significance level is based on parameter signif.level)
#'    \item \code{logs} warnings or messages of events, which happen during the run of our function
#'    \item \code{name_previous} values from the previous iterative step
#' }
#' In addition we get \code{contents}, which is list of outputs as original values and entered parameters, which are usable in other estimation process.
#' If we set up parameter \code{useNLS} to TRUE, the start.values which enters in the estimation process will be calculated by \code{\link{nlsLM}} function. Otherwise the \code{start.values} and \code{name_upgraded.start.val} will be the same.
#' In the case, if we add estimated parameters with NaN, NA, Inf or -Inf, we end up with error message.
#'
#' @note The symbol \code{pi} is reserved for Ludolf's constant. So naming one of the model´s parameters by this symbol results in constant entry of the model.
#'
#' @references Kubáček, L. and Kubáčková, L. (2000) \emph{Statistika a metrologie}, Univerzita Palackého v Olomouci.
#'
#'    Köning1, R., Wimmer, G. and Witkovský, V. (2014) \emph{Ellipse fitting by nonlinear constraints to demodulate quadrature homodyne interferometer signals and to determine the statistical uncertainty of the interferometric phase}, Measurement Science and Technology.
#'
#' @seealso \code{\link{NanoIndent.OEFPIL}} and \code{\link{nls}} for nonlinear least square package. And especially function \code{\link{nlsLM}}.
#'
#' @examples
#' ## We use data file: "silica2098.RData" which is part of OEFPIL package
#' colnames(silica2098) <- c('x','y')
#'
#' ## Preparing parameter for OEFPIL function
#' max.iter = 100
#' see.iter.val = F
#' th = 0.001
#' signif.level = 0.05
#' useNLS = T
#'
#'
#' m.start <- 1.5
#' hp.start <- 0.9
#' alpha.start <- 5
#' start.val <- list(alpha=alpha.start, m=m.start, hp=hp.start)
#' names(start.val) <- c("alpha", "m", "hp")
#'
#' ## Imputed formula
#' form <- y ~ alpha * (x - hp) ^ m
#'
#' ##Now we can use OEFPIL function with defined starting values
#' output.form <- OEFPIL(silica2098, form, start.val,  max.iter = max.iter, see.iter.val = see.iter.val, th = th, signif.level = signif.level, useNLS = useNLS)
#'
#' ## Finally we can use summary function
#' summary(output.form)
#'
#' @examples
#' library(MASS)
#' steamdata <- steam
#' colnames(steamdata) <- c("x","y")
#' n <- nrow(steamdata)
#' CM1 <- diag(rep(0.1,2*n))
#' CM2 <- diag(c(rep(0.2^2,n), rep(0.1^2,n)))
#'
#' ##Use of OEFPIL function with defined parameters
#' st1 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 5, b2 = 8, b3 = 200),
#'              CM1, useNLS = F)
#'
#'
#' ## Finally we can use summary function
#' summary(st1)
#'
#' ## Plot of estimated functions
#' plot(output.form, signif.level = signif.level)
#'
#' ##Now we look at the case when algorithm of OEPFIL function does not converge
#' startsteam <- list(b1 = 0.1, b2 = 5, b3 = 200)
#' st3 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 0.1, b2 = 5, b3 = 200),
#'               CM1, useNLS = F)
#'
#' ##If set up useNLS = T, we get better upgraded starting values and algorithm converges
#' startsteam <- list(b1 = 0.1, b2 = 5, b3 = 200)
#' st3 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 0.1, b2 = 5, b3 = 200),
#'               CM1, useNLS = T)

#'
#' @import MASS
#' @import minpack.lm
#' @import Deriv
#' @import matrixcalc
#' @import plyr
#' @import ggplot2
#'
#' @export



OEFPIL <- function(data, form, start.val, CM = diag(dim(data)[1] * 2),  max.iter = 100,
                   see.iter.val = F, save.file.name, th = .Machine$double.eps ^ (2 / 3),
                   signif.level = 0.05, useNLS = T) {
  ##   Optimum Estimate of Function Parameters by Iterated Linearization.
  ## Users will enter initial formula and list of starting values of all parameters from
  ## start.val. It's important to satisfy following conditions instructions:
  ## (1) Initial formula has left and right hand side, where on the left side we specify
  ##     dependant variable, i.e. : y ~ x + a + 3. Notation: ~ x + a + 1 is not possible.
  ## (2) In the entry data table "data" are all columns named. This names corresponds
  ##     to names in entered formula.
  ## (3) On the right hand side of the formula, there must be entered valid functional notation
  ##     - it´s not typical R-style formula, since we does model response with function f().
  ##     Correct notation should be y ~ f(x, a_1, ..., a_n), where y stands for dependant and
  ##     x for independant variable, a_1, ..., a_n are parameters and f() is differentiable
  ##     function (with respect to x and a_1,...,a_n). So our correct notation
  ##     is y ~ "funkcny predpis for y".
  ##     The symbol \code{pi} is reserved for Ludolf's constant. So naming one of the
  ##     model´s parameters by this symbol results in constant entry of the model.


  ######################################################################################
  ## Further parameters:
  ## data           . . . entry data file, where the columns must be named
  ## CM             . . . covariance matrix - does not change in iteration process
  ## max.iter       . . . maximum number of iterations (default = 100)
  ## see.iter.val   . . . logical value (default \code{TRUE}) indicating if we want to display and save partial results of algorithm.
  ## save.file.name . . . name of file, where user want to save iterations of algorithm in, if none
  ##                      is given, nothing is going to save; recommended format is .Rdata
  ## th             . . . threshold necessary for iteration's stoppage
  ## signif.level   . . . significance level for confidence interval
  ## useNLS         . . . T, for calculation of initial estimates with the function nlsLM()
  ##                      F, for calculation of initial estimates from the starting values parameters


  logs <- NA

  if (ncol(data) != 2) {
    stop("The data has to contain two named columns - one for a dependent variable and one for an independent variable!")
  }
  ## stopifnot( ncol(data) == 2 ) # Another variant of program stoppage

  if (IsListOK(data) == F) {

    n <- dim(data)[1]
    odr <- sort(unique(c(which(is.infinite(data[, 1])),
                         which(is.na(data[, 1])),
                         which(is.infinite(data[, 2])),
                         which(is.na(data[, 2])))))

    data <- data[ - odr,]
    CM <- CM[ - c(odr, odr + n), - c(odr, odr + n)]

    logg <- paste("The data rows ", paste(odr, collapse = ", "),
                  " contained NaN, NA, Inf or -Inf values.", "\n",
                  "These rows were removed and given covariance matrix was upgraded accordingly.", "\n", sep = "")
    message(logg) ## display error massage in console
    logs <- paste(na.omit(logs), logg, sep = "//")
  }

  if (length(colnames(data)) != 2 | sum(is.na(colnames(data))) != 0) {
    stop("Both columns of the data have to be named!")
  }
  ## Control, if both columns of data are named.

  if (is.formula(form) == FALSE) {
    stop("There has to be a formula as an input value.")
  }
  ## Control of formula correction

  if ( !(is.matrix(CM)) | !(all(dim(CM) == c(dim(data)[1] * 2, dim(data)[1] * 2))) | !is.positive.semi.definite(CM)) {

    logg <- paste("'CM' has to be a covariance matrix.", "\n",
                  "dim(CM) = c(dim(data)[1] * 2, dim(data)[1] * 2)", sep="")

    stop(logg)
  }
  ## Control of input of covariance matrix

  ## is.positive.semi.definite(CM) # We haven't checked that.

  if (th < .Machine$double.eps ^ (2 / 3)) {
    th <- .Machine$double.eps ^ (2 / 3)
    logg <- paste("Threshold has to be greater than ", signif(.Machine$double.eps^(2 / 3), 4),
                  ". It was set to this value instead. ", "\n", sep = "")
    message(logg) ## display error massage in console
    logs <- paste(na.omit(logs), logg, sep = "//")
  }
  ## Control of threshold

  if (!(is.logical(see.iter.val))) {
    stop("Parameter 'see.iter.val' has to be logical.")
  }
  if (!(is.logical(useNLS))) {
    stop("Parameter 'useNLS' has to be logical.")
  }
  ## Type control of parameters being logical


  input.form.string <- SafeDeparse(form)
  ## Transformation of formula to string

  main.func.call <- str2lang(strsplit(input.form.string, "~")[[1]][2])
  ## Extraction of functional prescription and creation of "call" object

  if (strsplit(input.form.string, "~")[[1]][1] == "") {
    stop("There has to be the name of a dependent variable on the left side of the formula!")
  }
  ## Formula must have prescribed left side and the dependant variable must be named
  ## I.e. y ~ ... is OK, but formula cannot start with ~ ...

  if ("pi" %in% all.vars(form)) {
    orig.vars <- all.vars(form)[- which(all.vars(form) == "pi")]
  } else {
    orig.vars <- all.vars(form)
  }
  ## If the function contain "pi" sign, we make sure that the function will take it an account
  ## as a constatnt and not one of fucntional parameters.


  l <- length(orig.vars) - 2
  ## l is a number of parameters; the number of derivativs which we want to calculate is l + 1

  dep.var.name <- orig.vars[1]

  if (!(dep.var.name %in% colnames(data))) {
    logg <- paste("The dependent variable  '", dep.var.name,
                  "' specified on the left side of a formula is not found in the data (columns of the data have to be named).",
                  sep="")
    stop(logg)
  }
  ## Control of presence of name of dependant variable

  col.dep.var <- which(colnames(data) == dep.var.name)
  col.idp.var <- which(colnames(data) != dep.var.name)

  idp.var.name <- colnames(data)[col.idp.var]

  y <- data[,col.dep.var]

  if (!(idp.var.name %in% orig.vars)) {
    logg <- paste("The variable '", idp.var.name,
                  "' is not found on the right side of a given formula.",
                  sep = "")
    stop(logg)
  }
  ## Control the presence of name of independant variable on the right hand side of formula.

  x <- data[,col.idp.var]

  if (which(orig.vars == idp.var.name) != 2) {
    idp.var.order <- which(orig.vars == idp.var.name)
    vars <- c(orig.vars[1], idp.var.name, orig.vars[-c(1, idp.var.order)])
  } else {
    vars <- orig.vars
  }
  ## We are changing the order of independant variable and parameters
  ## The first spot will take dependant variable, then independent and finally parameters.

  if (!(is.list(start.val))) {
    start.val <- as.list(start.val)
  }
  ## Transformation of starting values to type list, if needed.

  if (!(all(vars[3:length(vars)] %in% names(start.val)))) {
    jm.index <- which(!(vars[3:length(vars)] %in% names(start.val)))
    logg <- paste("There are some parameters on the right side of the formula, but there are no starting values for them: ",
                  paste(vars[2 + jm.index], collapse = ", "), ".", sep = "")
    stop(logg)
  }
  ## Test the presence of all starting values of parameters stated in formula

  if (!(all(names(start.val) %in% vars[3:length(vars)]))) {
    start.val[which(!(names(start.val) %in% vars[3:length(vars)]))] <- NULL
  }
  ## Removing starting values of parameters which are not present in formula.
  ## I.e. removing the starting values, which are extra (we have start value but
  ## no parameter in formula)

  pom.factor <- factor(vars[3:length(vars)], levels = names(start.val),
                       ordered = T)
  pom.factor <- sort(pom.factor)
  vars[3:length(vars)] <- as.vector(pom.factor)
  ## Ordering of parameter's name, in the order of users's list of starting values

  args <-  vector(mode = "list", length = (l + 1))
  names(args) <- vars[2:length(vars)]
  ## Creation of blank list; arguments are entry values of function

  # args[vars[2]] <- start.val[vars[2]]
  ## Eventual assigning of starting values

  MainFunction <- function(args, body, env = parent.frame()) {
    args <- as.pairlist(args)
    eval(call("function", args, body), env)
  }

  LOF <- list()
  LOF[[1]] <- MainFunction(args, main.func.call)

  LOF[2:(l+2)] <- sapply(1:(l+1), function(i) {
    Deriv(MainFunction(args, main.func.call), vars[i+1])
  })
  ## Derivative of entered function with respect to independent variable and parameters

  # Choice of modified starting values with usage of NLS
  if (useNLS == T) {

    lst.start.val <- tryCatch( expr = {

      nlc <- nls.control(maxiter = 1000)

      formula.input <- as.formula(input.form.string)

      nonlin_model <- nlsLM(formula.input, data = data, control = nlc,
                            start = start.val)
      ## If we set up "form" instead of "formula.input", it will fail

      coef.vec <- coefficients(nonlin_model)
      logs <- NA # warning protocol is set to be empty

      L <- list()
      L[names(coef.vec)] <- as.list(coef.vec)
      L[["CM0"]] <- vcov(nonlin_model)
      L[["logs"]] <- logs
      ## the output from this part (i.e. from nlsLM) is list containing all parameters

      L
      ## thanks to the final row is the output from this part list containing all the parameters
    }, error = function(e) {

      err_R <- NULL
      err_R <- conditionMessage(e) # save the error messages from R
      logg <- paste("Cannot use nls(), different start values were used",
                    "\n", "Original nls error: ", err_R, "\n", sep="")
      message(logg) ## display error massage in console
      logs <- paste(na.omit(logs), logg, sep = "//")

      L <- start.val
      CM0 <- NA
      L[["CM0"]] <- CM0
      L[["logs"]] <- logs

      return(L)
    })

    ## ErrM <- paste("Problem with nls(). Original error: ", err_R, sep="")
  } else {

    lst.start.val <- start.val
    CM0 <- NA
    logs <- NA
    lst.start.val[["CM0"]] <- CM0
    lst.start.val[["logs"]] <- logs

  }

  CM0 <- lst.start.val$CM0
  logs <- lst.start.val$logs

  # Function calculate iterations
  lst.iteration <- OEFPILIter(y0 = y, x0 = x, L = lst.start.val, CM = CM, max.iter = max.iter,
                              see.iter.val = see.iter.val, save.file.name = save.file.name,
                              th = th, LOF = LOF, logs = logs)

  contents <- list(input.form.string = input.form.string, LOF = LOF, x = x, y = y,
                   dep.var.name = dep.var.name, idp.var.name = idp.var.name,
                   names.of.parameters = names(lst.iteration$lst.parameters),
                   signif.level = signif.level)
  ## list contains rest of parameters usable in following functions

  names(contents)[3] <- idp.var.name
  names(contents)[4] <- dep.var.name

  lst.output <- list()
  length(lst.output) <- 3*l + 8
  lst.output[1:l] <- lst.iteration$lst.parameters
  lst.output[(l+1):(2*l)] <- lst.start.val[1:l]
  lst.output[[2*l+1]] <- lst.iteration$cov_m
  lst.output[[(2*l+2)]] <- CM0
  lst.output[[(2*l+3)]] <- lst.iteration$yEst
  lst.output[[(2*l+4)]] <- lst.iteration$xEst
  lst.output[[(2*l+5)]] <- lst.iteration$it_num
  lst.output[(2*l+6):(3*l+5)] <- lst.iteration$lst.parameters_previous.step
  lst.output[[(3*l+6)]] <- NA
  lst.output[[(3*l+7)]] <- contents
  lst.output[[(3*l+8)]] <- lst.iteration$logs

  names(lst.output) <- c(paste(names(lst.iteration$lst.parameters),"_Est",sep=""),
                         paste(names(lst.start.val)[1:l],"_upgraded.start.val",sep=""), "cov.m_Est", "cov.m_nlsLM",
                         paste(vars[1], "_Est", sep=""), paste(vars[2], "_Est", sep=""), "it_num",
                         paste(names(lst.iteration$lst.parameters_previous.step), "_previous.step", sep=""),
                         "CI_parameters", "contents", "logs")

  class(lst.output) <- c("OEFPIL", class(lst.output))

  lst.output[[(3*l+6)]] <- confint.OEFPIL(lst.output, signif.level = signif.level)

  return(lst.output)
}

################################################################################

OEFPILIter <- function(y0, x0, L, CM, max.iter = 100, see.iter.val = F,
                       save.file.name, th, LOF,  logs = NA) {
  ##   Iteration for the calculation of parametr estimation in generalised version of algorythm
  ## y0             . . . entry vector y (mesasure of strength)
  ## x0             . . . entry vector x (measure of depth)
  ## L              . . . list of original parameters estimations
  ## CM             . . . covariance matrix -  this part does not change in the process of
  ##                      algorythm
  ## max.iter       . . . maximum number of iterations of the algorithm (default = 100)
  ## see.iter.val   . . . T, if we want to continuously save and output results from the process of algorithm
  ## save.file.name . . . name of file, where user want to save iterations of algorithm in, if none
  ##                      is given, nothing is going to save; recommended format is .Rdata
  ## th             . . . treshold necessary for iteration's stoppage
  ## LOF            . . . entry list of user's defined functions
  ## logs           . . . string containing the messages and warnings about the whole
  ##                      process of function

  l <- length(L) - 2
  ## number of parameters

  L0 <- list()
  length(L0) <- l
  L0 <- L[1:l]
  ## original values

  L1 <- list()
  length(L1) <- l
  ## new values

  # parameter's value will not change in the iterative process; "L1" sign is current
  # "L0" sign is the value from the previous step
  # We set them to be same at the beginning
  L1 <- L0

  Q22 <- NA

  # x1 and y1 are values, which will be updated in the process
  x1 <- x0
  y1 <- y0

  N <- length(y0)

  # set of number of parameter for number of iterations; final number of iterations
  # so it will be considered with zero step
  it_num <- 0

  while((ConditionForIteration(L0, L1, th) || it_num == 0) && it_num < max.iter) {

    fval <-  sapply(x1, function(val, LP){do.call(LOF[[1]], args=c(val, LP))}, L1)
    ## It contains results of function of x1 defined by user, for parameters contained in L1


    b_vec <- y1 - fval

    B11 <- - diag(sapply(x1, function(val, LP){do.call(LOF[[2]], args=c(val, LP))}, L1))
    ## Quantification of derivatives of user's defined function with respect to x1 for
    ## value of parameters contained in L1 and set as a diagonal entries in matrix

    B1 <- cbind(B11, diag(N))

    if ((sum(is.infinite(B1)) != 0) || (sum(is.na(B1)) != 0)) {
      logg <- paste("There is NaN, NA, Inf or -Inf value in matrix B1. Iteration aborted!!!", "\n",
                    "The following results (if any) are only partial.", "\n", sep="")
      message(logg)
      logs <- paste(na.omit(logs), logg, sep = "//")
      break
    }
    ## while cyclus will be stopped in case of NA, NaN, Inf or -Inf value in B1

    B2 <- CreateB2(LOF, L1, x1, l)
    ## creation of matric B2

    if ((sum(is.infinite(B2)) != 0) || (sum(is.na(B2)) != 0)) {
      logg <- paste("There is NaN, NA, Inf or -Inf value in matrix B2. Iteration aborted!!!", "\n",
                    "The following results (if any) are only partial.", "\n", sep="")
      message(logg)
      logs <- paste(na.omit(logs), logg, sep = "//")
      break
    }
    ## while cyclus will be stopped in case of NA, NaN, Inf or -Inf value in B2

    M11 <- B1 %*% CM %*% t(B1) # CM is covariance matrix

    lst.chol <- tryCatch(expr = {
      Lmat <- t(chol(M11)) # Cholesky decomposition of matrix M11

      logg <- NA

      L <- list(Lmat = Lmat, logg = logg)
      ## output from this part is the final row (i.e. L)
    }, error = function(e) {

      err_R <- NULL
      err_R <- conditionMessage(e) # ulozeni chybove hlasky z R
      logg <- paste("Problems with computing Cholesky decomposition (chol() function). Iteration aborted!!!",
                    "\n", "Original chol error: ", err_R, "\n", sep="")
      message(logg) ## display error message in the console
      Lmat <- NA

      L <- list(Lmat = Lmat, logg = logg)
      return(L)
    })
    ## safe startup of function chol. We provide, so that the function will not break,
    ## just display warning.

    if (grepl("Iteration aborted!!!", lst.chol$logg)) {
      logs <- paste(na.omit(logs), lst.chol$logg, sep = "//")
      break
    }
    ## Interruption of cycle, in case of error message of the chol() function

    Lmat <- lst.chol$Lmat
    E <- forwardsolve(Lmat, B2) # Solves a triangular system of linear equations.
    Esing <- svd(E)
    Ue <- Esing$u
    Ve <- Esing$v
    Se <- diag(Esing$d)
    Seinv <- diag(1 / (Esing$d))
    Fmat <- Ve %*% Seinv
    G <- forwardsolve(t(Lmat), Ue)
    Q21 <- Fmat %*% t(G)
    Q11 <- chol2inv(Lmat) - G %*% t(G)
    Q22 <- - Fmat %*% t(Fmat)

    ############################################################################
    # Old algorythm with help of function ginv (with warning message) ###
    ############################################################################
    #
    # M <- rbind(cbind(M11, B2), cbind(t(B2), matrix(0,l,l)))
    # ## M is the matrix, which inverse is matrix Q
    #
    # if ((sum(is.infinite(M)) > 0) | (sum(is.na(M)) > 0)) {
    #   logg <- paste("Infinite, NA or NaN input into ginv(). Iteration aborted!!!", "\n",
    #                 "The following results (if any) are only partial.", "\n", sep="")
    #   message(logg)
    #   logs <- paste(na.omit(logs), logg, sep = "//")
    #   break
    # }
    # ## stoppage of while cycle in case of Inf, NaN or NA input in ginv
    #
    # Q <- ginv(M)
    # ## calculation of inverse matrix
    #
    # Q11 <- Q[1:N, 1:N]
    # Q22 <- Q[(N+1):(N+l),(N+1):(N+l)]
    # Q21 <- Q[(N+1):(N+l), 1:N]

    colnames(Q22) <- names(L1)
    rownames(Q22) <- names(L1)

    output01 <- (diag(2*N) - CM %*% t(B1) %*% Q11 %*% B1) %*% c(x0 - x1, y0 - y1) - CM %*% t(B1) %*% Q11 %*% b_vec
    xdiff_hat <- output01[1:N]
    ydiff_hat <- output01[(N+1):(2*N)]

    output02 <- -Q21 %*% B1 %*% c(x0 - x1, y0 - y1) - Q21 %*% b_vec
    ## vector containing increase of all l parameters in current step

    L0 <- L1

    x1 <- x1 + xdiff_hat
    y1 <- y1 + ydiff_hat
    ## update (improve) values of vectors x1,y1

    L1 <- as.list(unlist(L0) + output02)
    names(L1) <- names(L0)
    ## update (improvement) of values of the rest of the parametra; we have to save either
    ## original and updated values due to condition from the beginning of cycle

    if (IsListOK(L1) == F) {
      logg <- paste("During the iteration, the estimated parameter equals NaN, NA, Inf or -Inf.", "\n",
                    "Iteration aborted!!! The following results (if any) are only partial.", "\n", sep="")
      message(logg)
      logs <- paste(na.omit(logs), logg, sep = "//")
      break
    }
    ## control of correctness of updated (improved) values (L1)

    # The final number of iterations including zero step.
    it_num <- it_num + 1

    ## Condition which save and output all values of estimates in process of iteration
    ## if it's required by user.
    if (see.iter.val == T) {

      print(paste("it_num ", it_num, ": ###############################",
                  sep = ""))
      print(data.frame(L1, row.names = ""))
      print(- Q22)
      print("##########################################")
    }

    if (!missing(save.file.name)) {
      list.to.save <- list(lst.parameters = L1, cov_m = - Q22, yEst = y1, xEst = x1,
                           it_num = it_num)
      list.name <- paste("semi.results_", it_num, sep = "")
      assign(list.name, list.to.save)
      AddObjectToRdata(list.to.save, list.name, rda_file = save.file.name, overwrite = T)

    }
  }

  return(list(lst.parameters = L1, cov_m = -Q22, yEst = y1, xEst = x1,
              it_num = it_num, lst.parameters_previous.step = L0, logs = logs))
}



################################################################################

SafeDeparse <- function(expr){
  ## Function for "safe" split of expression (in our case formula) in the way, that in the final
  ## string will be no jump over the next line


  ret <- paste(deparse(expr), collapse = "")
  ## Add more rows into one with whitespaces

  # Removing of "whitespace"
  gsub("[[:space:]][[:space:]]+", " ", ret)
}

################################################################################

ConditionForIteration <- function(L0, L1, th) {
  ## The output is logical value TRUE, in case the condition of the input into "do while"
  ## is satisfied.
  ## L0 . . . list containing values of parameters from previous step
  ## L0 . . . list containing current values of parameter
  ## th . . . threshold

  if (sum(abs(unlist(L1) - unlist(L0)) / abs(unlist(L0)) > th) > 0) {
    output <- T
  } else {
    output <- F
  }

  return(output)
}

###############################################################################

CreateB2 <- function(LOF, LP, xvec, m){
  ## Function creates B2 matrix in algorithm
  ## LOF  ... List of all functions passed by user.
  ##          Derivatives with respect to parameter are on the positions 3 and later in the correct order
  ## LP   ... list containing values of parameters in the correct order
  ## xvec ... vector of x values
  ## m    ... number of parameters

  B2 <- sapply(1:m, function(i) {
    sapply(xvec, function(val, LPP){do.call(LOF[[2+i]], args=c(val, LPP))}, LP)
  })

  return(-B2)
}

################################################################################

IsListOK <- function(List) {
  ## Function which checks in list data structure for presence of NaN, NA, Inf or -Inf.
  ## It can be used either for list's components as well as for particular elements of
  ## components of list.

  pv <- sum(sapply(List, `%in%`, x = NA), sapply(List, `%in%`, x = NaN),
            sapply(List, `%in%`, x = Inf), sapply(List, `%in%`, x = - Inf))

  if (pv == 0) {

    return(T)

  } else {

    return(F)
  }
}

################################################################################

AddObjectToRdata <- function(obj, name_obj, rda_file, overwrite = FALSE) {
  ## Function adds on another objects into the RData file. Copied (and modified) from:
  ## https://stackoverflow.com/questions/38364745/add-a-data-frame-to-an-existing-rdata-file
  .dummy <- NULL
  if (!file.exists(rda_file)) save(.dummy, file = rda_file)

  old_e <- new.env()
  new_e <- new.env()

  load(file = rda_file, envir = old_e)

  # name_obj <- deparse(substitute(obj))   # get the name of the object
  # from the previous (internet) version

  # new_e[[name_obj]] <- get(name_obj)     # use this only outside a function
  new_e[[name_obj]] <- obj

  # merge object from old environment with the new environment
  # ls(old_e) is a character vector of the object names
  if (overwrite) {

    # the old variables take precedence over the new ones
    invisible(sapply(ls(new_e), function(x)
      assign(x, get(x, envir = new_e), envir = old_e)))

    # And finally we save the variables in the environment
    save(list = ls(old_e), file = rda_file, envir = old_e)
  }
  else {
    invisible(sapply(ls(old_e), function(x)
      assign(x, get(x, envir = old_e), envir = new_e)))

    # And finally we save the variables in the environment
    save(list = ls(new_e), file = rda_file, envir = new_e)
  }
}

