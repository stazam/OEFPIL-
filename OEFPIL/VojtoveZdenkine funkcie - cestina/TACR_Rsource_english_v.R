################################################################################
# Funkce ###
################################################################################

library(MASS)
library(minpack.lm)
library(plyr) ## nutne pro funkci is.formula()
library(matrixcalc) ## nutne pro is.positive.semi.definite()
library(Deriv)

################################################################################

OEFPIL <- function(data, form, start.val, CM = diag(dim(data)[1] * 2),  max.iter = 100,
                   see.iter.val = F, save.file.name, th = 0.001, signif.level = 0.05,
                   useNLS = T) {
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


  #is.positive.definite(A)

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
    puvodni.jmena <- all.vars(form)[- which(all.vars(form) == "pi")]
  } else {
    puvodni.jmena <- all.vars(form)
  }
  ## If the function contain "pi" sign, we make sure that the function will take it an account
  ## as a constatnt and not one of fucntional parameters.


  l <- length(puvodni.jmena) - 2
  ## l is a number of parameters; the number of derivativs which we want to calculate is l + 1

  dep.var.name <- puvodni.jmena[1]

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

  if (!(idp.var.name %in% puvodni.jmena)) {
    logg <- paste("The variable '", idp.var.name,
                  "' is not found on the right side of a given formula.",
                  sep="")
    stop(logg)
  }
  ## Control the presence of name of independant variable on the right hand side of formula.

  x <- data[,col.idp.var]

  if (which(puvodni.jmena == idp.var.name) != 2) {
    idp.var.order <- which(puvodni.jmena == idp.var.name)
    jmena <- c(puvodni.jmena[1], idp.var.name, puvodni.jmena[-c(1, idp.var.order)])
  } else {
    jmena <- puvodni.jmena
  }
  ## We are changing the order of independant variable and parameters
  ## The first spot will take dependant variable, then independent and finally parameters.

  if (!(is.list(start.val))) {
    start.val <- as.list(start.val)
  }
  ## Transformation of starting values to type list, if needed.

  if (!(all(jmena[3:length(jmena)] %in% names(start.val)))) {
    jm.index <- which(!(jmena[3:length(jmena)] %in% names(start.val)))
    logg <- paste("There are some parameters on the right side of the formula, but there are no starting values for them: ",
                  paste(jmena[2 + jm.index], collapse = ", "), ".", sep = "")
    stop(logg)
  }
  ## Test the presence of all starting values of parameters stated in formula

  if (!(all(names(start.val) %in% jmena[3:length(jmena)]))) {
    start.val[which(!(names(start.val) %in% jmena[3:length(jmena)]))] <- NULL
  }
  ## Removing starting values of parameters which are not present in formula.
  ## I.e. removing the starting values, which are extra (we have start value but
  ## no parameter in formula)


  # jmena.sh <- names(start.val) ## the older version, in case of non functionality of down row
  pom.factor <- factor(jmena[3:length(jmena)], levels = names(start.val),
                       ordered = T)
  pom.factor <- sort(pom.factor)
  jmena[3:length(jmena)] <- as.vector(pom.factor)
  ## Ordering of parameter's name, in the order of users's list of starting values

  args <-  vector(mode = "list", length = (l + 1))
  names(args) <- jmena[2:length(jmena)]
  ## Creation of blank list; arguments are entry values of function

  # args[jmena[2]] <- start.val[jmena[2]]
  ## Eventual assigning of starting values

  MainFunction <- function(args, body, env = parent.frame()) {
    args <- as.pairlist(args)
    eval(call("function", args, body), env)
  }

  LOF <- list()
  LOF[[1]] <- MainFunction(args, main.func.call)

  LOF[2:(l+2)] <- sapply(1:(l+1), function(i) {
    Deriv(MainFunction(args, main.func.call), jmena[i+1])
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

  # CI.matrix <- CIOEFPIL(lst.iteration$lst.parameters, lst.iteration$cov_m,
  #                          signif.level = signif.level)
  ## Function calculates confidence interavals for parameters

  contents <- list(input.form.string = input.form.string, LOF = LOF, x = x, y = y,
                   dep.var.name = dep.var.name, idp.var.name = idp.var.name,
                   names.of.parameters = names(lst.iteration$lst.parameters),
                   signif.level = signif.level)
  ## list contains rest of parameters usable in following functions

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
                         paste(jmena[1], "_Est", sep=""), paste(jmena[2], "_Est", sep=""), "it_num",
                         paste(names(lst.iteration$lst.parameters_previous.step), "_previous.step", sep=""),
                         "CI_parameters", "contents", "logs")

  class(lst.output) <- "OEFPIL"

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

confint.OEFPIL <- function(output.form, signif.level = output.form$contents$signif.level) {
  ## Function calculate confidence intervals for parametrs counted by OEFPIL funciton.


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

################################################################################

plot.OEFPIL <- function(output.form, xx, signif.level = 0.05,
                        add.plot = T) {
  ## Function plots confidence bands of list from OEFPIL() function.
  ## output.form . . . output from OEFPIL()
  ## xx          . . . in these points we calculate and plot CI (confidence intervals) or
  ##                    CB (conf. bands)



  LOF <- output.form$contents$LOF ## list of functions
  x <- output.form$contents$x ## x-ova data
  y <- output.form$contents$y ## y-ova data

  cov_m <- output.form$cov.m_Est ## estimate of covariance matrix
  l <- dim(cov_m)[1] ## number of parameters

  lst.parameters <- output.form[1:l] ## parameter estimation
  names(lst.parameters) <- output.form$contents$names.of.parameters

  lst.parameters_previous.step <- output.form[(2*l+6):(3*l+5)]
  names(lst.parameters_previous.step) <- output.form$contents$names.of.parameters
  ## estimate from the previous step

  dep.var.name <- output.form$contents$dep.var.name ## name of dependant variabe
  idp.var.name <- output.form$contents$idp.var.name ## name of independant variable

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

    sl <- sort(c(signif.level/2, 1 - signif.level/2), decreasing = F)

    d <- length(signif.level)
    k <- length(xx)

    PCB_lwr <- matrix(rep(yy, d), k, d) + matrix(rep(qnorm(sl[1:d]), k), k, d, byrow = T) * sqrt(variance)
    PCB_upr <- matrix(rep(yy, d), k, d) + matrix(rep(qnorm(sl[(d+1):(2*d)]), k), k, d, byrow = T) * sqrt(variance)
    ## pointwise confidence band

    if (add.plot == T) {
      plot(x, y, xlab = idp.var.name, ylab = dep.var.name)

      lines(xx, yy, lwd = 2, col = "black")
      for (i in 1:d) {
        lines(xx, PCB_lwr[, d - i + 1], lwd = 2, lty = 2, col = i + 1)
        lines(xx, PCB_upr[, i], lwd = 2, lty = 2, col = i + 1)
      }

    }

    PointwiseCB <- cbind(PCB_lwr, PCB_upr)
    colnames(PointwiseCB) <- paste(round(sl * 100, 2), "%")

    return(invisible(list(PointwiseCB = PointwiseCB)))
  }
}

################################################################################

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
    }
  }

  output <- list(param_Est = Param.mat[,1], sd =  Param.mat[,2],
                 cov.m_Est = summary.form$cov.m_Est, it_num = summary.form$it_num,
                 CI_parameters = Param.mat[,3:(2 + d)])

  return(invisible(output))
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

NanoIndent <- function(data, alpha.start, m.start, hp.start, unload.data = F,
                       ucut = 0.98, lcut = 0.4, CM, uh = 0.5, uF = 0.001,
                       max.iter = 100, see.iter.val = F, save.file.name,
                       th = 0.001, signif.level = 0.05, useNLS = T) {
  ##  Function for calculation of nanoindentation. It uses OEPFIL, which she prepared
  ##data (cut them off) and formula for. Funtional dependence of parameters is fixed in the form:
  ##f = alpha * (h - hp) ^ m, where f is strength and h depth (on the unloaded part of curve).
  ## User has an option to enter his own starting values of parameters. If he does not make it,
  ## function will warn him, that the starting values will be set as stated in the algorithm.
  ## Function NanoIndent has some extra features compared to OEFPIL. It can output some extra
  ## warnings. For example if "hp >= min(h)" (does not make sence from phicical prospective).
  ## unload.data . . . F, in case we have both complete load and unload curve, in this case the function will make cut off by its own
  ##                   T, if we know, that we only have unloaded curve;
  ## ucut      .  .  . it's the upper bound of cut off F_max
  ## lcut      .  .  . it's the lower (?? nizsi ??) bound of cut off F_max, i.e. if ucut = 0.98,
  ##                   lcut = 0.4 we consider 40 - 98 % from F_max (standard/norm recommendation)

  logs <- NA

  ICL <- FindIndentCurve(data, lcut = lcut, ucut = ucut, unload.data = unload.data)
  ## If we do not use "clean" data, we have to use cut off function
  ## If we have unload function, some cut off must be used

  data.cut <- data.frame(ICL$f_orig, ICL$h_orig)
  colnames(data.cut) <- c("f", "h")
  ## Creation of data file, which has to have named columns - important for future
  ##  work with OEFPIL variable.

  Fmax <- ICL$Fmax
  hmax <- ICL$hmax
  hmin <- ICL$hmin

  # If any of starting values are missing, we will use our
  if (missing(alpha.start) || missing(m.start) || missing(hp.start)) {

    index.miss <- which(c(missing(alpha.start), missing(m.start), missing(hp.start)))
    par.miss <- paste(c("alpha", "m", "hp")[index.miss], collapse = ", ")


    # starting values for alpha, hp and m by proposition of algorithm
    # (choice by proposition of algorithm and correct from the physical prosepctive)
    m.start <- 1.5
    hp.start <- 0.9 * hmin
    alpha.start <- Fmax / ((hmax - hp.start) ^ m.start)
    start.val <- list(alpha = alpha.start, m = m.start, hp = hp.start)

    logg <- paste("Starting values for these parameters were not given: ",
                  par.miss, ".", "\n", "These starting values were used instead:\n",
                  "alpha.start = ", alpha.start, "\n","m.start     = ", m.start,
                  "\n", "hp.start    = ", hp.start, "\n", sep = "")
    message(logg)
    logs <- paste(na.omit(logs), logg, sep = "// ")

  } else {
    start.val <- list(alpha = alpha.start, m = m.start, hp = hp.start)
  }

  form <- f ~ alpha * (h - hp) ^ m

  if (missing(CM)) {
    ## If the user does not define covariance matrix, we will use this one

    vec.uni <- rep(1, length(data.cut$h))
    CM <- CovMatrix(uh, uF, vec.uni, vec.uni)
    ## covariance matrix calculated by the proposition of algorythm
    ## Does not change in the provess of estimation.
  }

  output.form <- OEFPIL(data.cut, form, start.val, CM = CM, max.iter = max.iter,
                        see.iter.val = see.iter.val, save.file.name = save.file.name,
                        th = th, signif.level = signif.level, useNLS = useNLS)

  if (IsListOK(output.form$hp_Est) && IsListOK(output.form$h_Est)) {

    if (output.form$hp_Est >= min(output.form$h_Est)) {
      logg <- "Warning: After finishing the iteration, the value of hp is greater or equal than min of upgraded h values. \n"
      message(logg)
      output.form$logs <- paste(na.omit(c(logs, output.form$logs)), logg, sep = "//")
    }

  }

  # Warningy pro parametry m a hp z ParEst
  # if (m > 20) {
  #   logg <- "Warning: After finishing the iteration, the value of m is greater than 20. \n"
  #   message(logg)
  #   logs <- paste(logs, logg, sep = "")
  # }
  #
  # if (m < 0.001) {
  #   logg <- "Warning: After finishing the iteration, the value of m is lower than 0.001 \n"
  #   message(logg)
  #   logs <- paste(logs, logg, sep = "")
  # }
  #
  # if (hp >= min(h)) {
  #   logg <- "Warning: After finishing the iteration, the value of hp is greater or equal than min of upgraded h values. \n"
  #   message(logg)
  #   logs <- paste(logs, logg, sep = "")
  # }

  # output.form$input_cut_data <- data.cut

  return(output.form)
}

################################################################################

FindIndentCurve <- function(data, lcut = 0.4, ucut = 0.98, unload.data = F) {
  ## Function calculate the origin and the end of unloaded curve and cut off from this curve.
  ## First option is, if we have both loaded and unloaded data:
  ## (1) We consider, that the origin is the point, where the maximum strength is measured.
  ## If there are is more than one point, we consider the last of these points.
  ## Measurement containing negative values of force and depth are cutted off. The last step
  ## is to cut off the curve by the percentage ranfe from F_max, entered by lcut, ucut
  ## (2) We have got data just from unloaded curve (unload.data = T). We cut off
  ## the negative values and afterwards modification by ucut and lcut.
  ## data        . . . entry data file; 1. column is h, 2. col. is F
  ## unload.data . . . F, in case we have both complete load and unload curve, in this case the function will make cut off by its own
  ##                   T, if we know, that we only have unloaded curve;
  ## ucut      .  .  . it's the upper bound of cut off F_max
  ## lcut      .  .  . it's the lower (?? nizsi ??) bound of cut off F_max, i.e. if ucut = 0.98,
  ##                   lcut = 0.4 we consider 40 - 98 % from F_max (standard/norm recommendation)

  if (unload.data == F) {

    # cutting off the unload curve
    pos_Fmax <- dim(data)[1] - which.max(rev(data[,2])) + 1
    Fmax <- data[pos_Fmax,2]
    hmax <- data[pos_Fmax,1]
    data <- data[pos_Fmax:dim(data)[1],]

    # cutting off the unload curve - the old algorythm
    #pos_hmax <- dim(data)[1] - which.max(rev(data[,1])) + 1
    #hmax <- data[pos_hmax,1]
    #Fmax <- data[pos_hmax,2]
    #data <- data[pos_hmax:dim(data)[1],]

    # cutting off the negative values
    if (any(data < 0)) {
      mz <- min(which(data[,1] < 0), which(data[,2] < 0))
      data <- data[1:(mz - 1),]
    }

  } else {

    # cutting off the negative values
    if (any(data < 0)) {
      mz <- min(which(data[,1] < 0), which(data[,2] < 0))
      data <- data[1:(mz - 1),]
    }

    pos_Fmax <- 1
    Fmax <- data[pos_Fmax,2]
    hmax <- data[pos_Fmax,1]

  }

  hmin <- min(data[,1])

  # cutting off the right part of the unload  curve
  UB <- Fmax * ucut
  LB <- Fmax * lcut

  index <- which( data[,2] <= UB & data[,2] >= LB )

  f_orig <- data[index,2]
  h_orig <- data[index,1]

  return(list(f_orig = f_orig, h_orig = h_orig, Fmax = Fmax, hmax = hmax, hmin = hmin))
}

################################################################################

CovMatrix <- function(uh, uF, hvec, Fvec) {
  # uh, uF - uncertainties of depth and load
  # Fvec - vector of measured values of loading
  ## output: covariance matrix

  N <- length(Fvec)
  M0 <- matrix(0, nrow = N, ncol = N)
  M1 <- uh^2 * diag(hvec)
  M2 <- uF^2 * diag(Fvec)^2
  covM <- cbind(rbind(M1, M0), rbind(M0, M2))
  return(covM)
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
  ## Function creates B2 matrix in algorythm
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

################################################################################
