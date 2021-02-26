#' @name NanoIndent.OEFPIL
#' @title Calculation of pestimation of parameters of unloaded curve during Nanoindentation
#' @description Function calculates parameters of unloaded curve during nanoindentation by iterated linearization.
#' @usage NanoIndent.OEFPIL(data, alpha.start, m.start, hp.start, unload.data = F,ucut = 0.98,
#'                          lcut = 0.4, CM, uh = 0.5, uF = 0.001, max.iter = 100, see.iter.val = F,
#'                          save.file.name, th = .Machine$double.eps ^ (2 / 3), signif.level = 0.05,
#'                          useNLS = T)
#'
#' @param data Data file can be any object of type \code{data.frame} with 2 columns or \code{list} with 2 elements.
#' @param form an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param name.start starting values of estimating parameters.
#' @param unload.data logical value (default FALSE ), which we set up to F in the case if we have both loaded and unloaded curve. In case of TRUE, we have only unloaded curve.
#' @param ucut numerical value, indicating the upper bound of cut off.
#' @param lcut numerical value, indicating the lower bound of cut off.
#' @param CM covariance matrix (does not change in iteration process).
#' @param uh standard deviation of depth
#' @param uF standard deviation of force
#' @param max.iter maximum number of iterations
#' @param see.iter.val logical value (default \code{TRUE}) indicating if we want to display and save partial results of algorithm.
#' @param save.file.name name of the file for saving results.
#' @param th numerical value, indicating threshold necessary for estimation stoppage.
#' @param signif.level significance level for confidence interval
#' @param useNLS logical. If \code{TRUE} (the default value), function will set up starting parameters calculated by \code{\link{nlsLM}} function (nonlinear least square estimation).
#'
#' @details Function prepared data (clean and cut them off) for \code{\link{OEPFIL}}, for parameter estimation.
#'  Functional dependence of parameters is fixed in
#' the form: \code{f = alpha * (h - hp) ^ m}, where \code{f} is strength and \code{h} depth
#' (on the unloaded part of curve). User has an option to enter his own starting values of
#' parameters, or the values will be set up by algorithm.
#'
#' In case, we do not add our own covariance matrix, it will be calculated by the algorithm with use of \code{u_h} and \code{u_f}.
#'
#' @return Returns an object of class \code{'OEPFIL'} is a list containing at least the following components
#' \itemize{
#'    \item name_est estimations of model parameters
#'    \item name_upgraded.start.val estimations of model parameters
#'    \item it_num number of iterations
#'    \item CI_parameters list of confidence intervals for estimated parameters (significance level is based on parameter signif.level)
#'    \item logs warnings or messages of events, which happen during the run of our function
#'    \item \code{name_previous} values from the previous iterative step
#' }
#' In addition we get \code{contents}, which is list of outputs as original values and entered parameters, which are usable in other estimation process.
#' If we set up parameter \code{useNLS} to TRUE, the start.values which enters in the estimation process will be claculated by \code{\link{nlsLM}} function. Otherwise the \code{start.values} and \code{name_upgraded.start.val} will be the same.
#' In the case, if we add estimated parameters with NaN, NA, Inf or -Inf, we end up with error message.
#'
#' @references Kubáček, L. and Kubáčková, L. (2000) \emph{Statistika a metrologie}, Univerzita Palackého v Olomouci.
#'
#'    Köning1, R., Wimmer, G. and Witkovský, V. (2014) \emph{Ellipse fitting by nonlinear constraints to demodulate quadrature homodyne interferometer signals and to determine the statistical uncertainty of the interferometric phase}, Measurement Science and Technology.
#'
#' @seealso \code{\link{OEFPIL}}
#'
#' @examples
#' ## We use "uncut" data file: "silicaberk.txt" which is part of OEFPIL package
#' ## Preparing parameter for OEFPIL function
#' unload.data = T
#' ucut = 0.98
#' lcut = 0.2
#' uh = 0.5
#' uF = 0.001
#' max.iter = 100
#' see.iter.val = F
#' th = 0.001
#' signif.level = 0.05
#' useNLS = T
#'
#'
#' ##use of NanoIndent function with default parameters
#' output.form.NI <- NanoIndent.OEFPIL(data)
#'
#' ##use of summary function for clear output
#' summary(output.form.NI)
#'
#'
#' ##use of NanoIndent function with our predefined parameters
#' output.form.NI <- NanoIndent.OEFPIL(silicaBerk, unload.data = unload.data, ucut = ucut,
#'                             lcut = lcut, uh = uh, uF = uF, max.iter = max.iter,
#'                             see.iter.val = see.iter.val, th = th,
#'                             signif.level = signif.level, useNLS = useNLS)
#'
#' ##use of summary function for clear output
#' summary(output.form.NI)
#'
#' ##plot of estimated functions
#' plot(output.form.NI, signif.level = signif.level)
#'
#' @export


################################################################################

NanoIndent.OEFPIL <- function(data, alpha.start, m.start, hp.start, unload.data = F,
                       ucut = 0.98, lcut = 0.4, CM, uh = 0.5, uF = 0.001,
                       max.iter = 100, see.iter.val = F, save.file.name,
                       th = .Machine$double.eps ^ (2 / 3), signif.level = 0.05,
                       useNLS = T) {
  ##  Function for calculation of nanoindentation. It uses OEPFIL, which she prepared
  ## data (cut them off) and formula for. Funtional dependence of parameters is fixed in the form:
  ## f = alpha * (h - hp) ^ m, where f is strength and h depth (on the unloaded part of curve).
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
  ##   The function ignores NA values in the data.
  ##
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
      mz <- min(which(data[,1] < 0), which(data[,2] < 0), na.rm = T)
      data <- data[1:(mz - 1),]
    }

  } else {

    # cutting off the negative values
    if (any(data < 0)) {
      mz <- min(which(data[,1] < 0), which(data[,2] < 0), na.rm = T)
      data <- data[1:(mz - 1),]
    }

    pos_Fmax <- 1
    Fmax <- data[pos_Fmax,2]
    hmax <- data[pos_Fmax,1]

  }

  hmin <- min(data[,1], na.rm = T)

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



