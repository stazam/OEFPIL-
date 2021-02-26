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
  ## Uzivatel zada vstupni formuli a startovaci list hodnot vsech parametru
  ## start.val. Je dulezite, aby bylo splneno:
  ## (1) Zadana formule ma levou a pravou stranu, kdy na leve strane je specifikovana
  ##     zavisla promenna, napr.: y ~ x + a + 3. Zapis: ~ x + a + 1 neni mozny.
  ## (2) Ve vstupni tabulce "data" jsou pojmenovany sloupecky. Tato jmena odpovidaji
  ##     znaceni v zadane formuli.
  ## (3) Na prave strane formule musi byt zadan validni funkcni predpis - nejde o formuli
  ##     v klasickem R-kovskem smyslu, jelikoz modelujeme odezvu pomoci funkce.
  ##     Korektni zapis by tedy byl y ~ f(x, a_1, ..., a_n), kde y je zavisla a x nezavisla
  ##     promenna, a_1, ..., a_n jsou parametry a f() je diferencovatelna funkce (podle x
  ##     i podle a_1, ..., a_n). Nas zapis vsak je y ~ "funkci predpis pro y".
  ##     Jmeno parametru "pi" je rezervovano pro danou konstantu. Pokud uzivatel
  ##     pojmenuje nejaky parametr "pi", funkce ho nebude brat jako parametr, ale jako
  ##     Ludolfovo cislo.

  ######################################################################################
  ## Dalsi parametry:
  ## data           . . . vstupni soubor dat, musi byt pojmenovany sloupce
  ## CM             . . . kov. matice - v prubehu alg. se nemeni
  ## max.iter       . . . maximalni pocet iteraci algoritmu (default = 100)
  ## see.iter.val   . . . T, pokud chceme vypisovat a ukladat prubezne vysledky v prubehu iterace
  ## save.file.name . . . jmeno souboru do ktereho chce uzivtel ukladat prubeh iterace; pokud jej
  ##                      nezadame, neexportuje se nic; doporuceny format: .Rdata
  ## th             . . . treshold nutny pro zastaveni iterace
  ## signif.level   . . . hladina vyznamnosti pro IS
  ## useNLS         . . . T, pro vypocet pocatecnich odhadu pomoci nlsLM()
  ##                      F, nastavi poc. odhady z dat podle start.val

  logs <- NA

  if (ncol(data) != 2) {
    stop("The data has to contain two named columns - one for a dependent variable and one for an independent variable!")
  }
  ## stopifnot( ncol(data) == 2 ) # Jina varianta zastaveni programu

  if (length(colnames(data)) != 2 | sum(is.na(colnames(data))) != 0) {
    stop("Both columns of the data have to be named!")
  }
  ## Overujeme, ze oba sloupecky dat jsou pojmenovany.

  if (is.formula(form) == FALSE) {
    stop("There has to be a formula as an input value.")
  }
  ## Kontrola, ze jde o formuli

  if ( !(is.matrix(CM)) | !(all(dim(CM) == c(dim(data)[1] * 2, dim(data)[1] * 2))) | !is.positive.semi.definite(CM)) {

    logg <- paste("'CM' has to be a covariance matrix.", "\n",
                  "dim(CM) = c(dim(data)[1] * 2, dim(data)[1] * 2)", sep="")

    stop(logg)
  }
  ## Kontrola zadani kovariancni matice


  #is.positive.definite(A)

  if (!(is.logical(see.iter.val))) {
    stop("Parameter 'see.iter.val' has to be logical.")
  }
  if (!(is.logical(useNLS))) {
    stop("Parameter 'useNLS' has to be logical.")
  }
  ## Kontrola zda jsou dane parametry logickeho typu


  input.form.string <- SafeDeparse(form)
  ## Prevod zadane formule na "string"

  main.func.call <- str2lang(strsplit(input.form.string, "~")[[1]][2])
  ## vyjmuti funkcniho predpisu a vytvoreni objektu "call"

  if (strsplit(input.form.string, "~")[[1]][1] == "") {
    stop("There has to be the name of a dependent variable on the left side of the formula!")
  }
  ## Formule musi mit uvedenou levou stranu a musi byt identifikovana zavisla promenna.
  ## Tj. y ~ ... je v poradku, ale formule nemuze zacinat ~ ...

  if ("pi" %in% all.vars(form)) {
    puvodni.jmena <- all.vars(form)[- which(all.vars(form) == "pi")]
  } else {
    puvodni.jmena <- all.vars(form)
  }
  ## Pokud se ve funkci vyskytuje "pi", timto osetrime, aby toto "pi" funkce nezahrnula
  ## mezi parametry, ale nahlizela na nej jako na konstantu.

  l <- length(puvodni.jmena) - 2
  ## l je pocet parametru; pocet derivaci, ktere chceme spocitat, je vsak l + 1

  dep.var.name <- puvodni.jmena[1]

  if (!(dep.var.name %in% colnames(data))) {
    logg <- paste("The dependent variable  '", dep.var.name,
                  "' specified on the left side of a formula is not found in the data (columns of the data have to be named).",
                  sep="")
    stop(logg)
  }
  ## Podminka na overeni pritomnosti jmena zavisle promenne

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
  ## Podminka na overeni pritomnosti jmena nez. prom. na prave strane formule.

  x <- data[,col.idp.var]

  if (which(puvodni.jmena == idp.var.name) != 2) {
    idp.var.order <- which(puvodni.jmena == idp.var.name)
    jmena <- c(puvodni.jmena[1], idp.var.name, puvodni.jmena[-c(1, idp.var.order)])
  } else {
    jmena <- puvodni.jmena
  }
  ## Prehazujeme v poradi nezavislou promennou a parametry.
  ## Na 1. miste chceme zavislou prom., na 2. nezavislou prom., pak parametry.

  if (!(is.list(start.val))) {
    start.val <- as.list(start.val)
  }
  ## Prevod startovacich hodnot do formatu list, pokud uz v nem nejsou

  if (!(all(jmena[3:length(jmena)] %in% names(start.val)))) {
    jm.index <- which(!(jmena[3:length(jmena)] %in% names(start.val)))
    logg <- paste("There are some parameters on the right side of the formula, but there are no starting values for them: ",
                  paste(jmena[2 + jm.index], collapse = ", "), ".", sep = "")
    stop(logg)
  }
  ## Overujeme podminku, ze vsechny parametry uvedene ve formuli na prave strane,
  ## maji zadanou startovaci hodnotu uzivatelem.

  if (!(all(names(start.val) %in% jmena[3:length(jmena)]))) {
    start.val[which(!(names(start.val) %in% jmena[3:length(jmena)]))] <- NULL
  }
  ## Odstraneni startovacich hodnot tech parametru, jenz se nevyskytuji ve formuli.
  ## Tj. odstarneni startovacich hodnot, ktere jsou navic (mame st. hodnotu, ale
  ## nemame parametr).

  # jmena.sh <- names(start.val) ## starsi varianta, pokud by radek nize nefungoval
  pom.factor <- factor(jmena[3:length(jmena)], levels = names(start.val),
                       ordered = T)
  pom.factor <- sort(pom.factor)
  jmena[3:length(jmena)] <- as.vector(pom.factor)
  ## Serazeni jmen parametru podle poradi, jak je uzivatel zadal do listu st. hodnot

  args <-  vector(mode = "list", length = (l + 1))
  names(args) <- jmena[2:length(jmena)]
  ## vytvoreni prazdneho listu; argumenty jsou vstupni hodnoty funkce

  # args[jmena[2]] <- start.val[jmena[2]]
  ## pripadne prirazeni startovacich hodnot

  MainFunction <- function(args, body, env = parent.frame()) {
    args <- as.pairlist(args)
    eval(call("function", args, body), env)
  }

  LOF <- list()
  LOF[[1]] <- MainFunction(args, main.func.call)

  LOF[2:(l+2)] <- sapply(1:(l+1), function(i) {
    Deriv(MainFunction(args, main.func.call), jmena[i+1])
  })
  ## derivace zadane funkce podle nezavisle promenne a podle parametru

  # volba upravenych startovacich hodnot s vyuzitim NLS
  if (useNLS == T) {

    lst.start.val <- tryCatch( expr = {

      nlc <- nls.control(maxiter = 1000)

      formula.input <- as.formula(input.form.string)

      nonlin_model <- nlsLM(formula.input, data = data, control = nlc,
                            start = start.val)
      ## Pokud bychom zadali "form" misto "formula.input", nebude to fungovat.
      ## Netusim proc.

      coef.vec <- coefficients(nonlin_model)
      logs <- NA # protokol hlaseni je nastaven na prazdnou hodnotu

      L <- list()
      L[names(coef.vec)] <- as.list(coef.vec)
      L[["CM0"]] <- vcov(nonlin_model)
      L[["logs"]] <- logs
      ## vystup z teto casti (tj. z nlsLM) je list obsahujici vsechny parametry

      L
      ## diky poslednimu radku je vystup z teto casti list obsahujici vsechny parametry
    }, error = function(e) {

      err_R <- NULL
      err_R <- conditionMessage(e) # ulozeni chybove hlasky z R
      logg <- paste("Cannot use nls(), different start values were used",
                    "\n", "Original nls error: ", err_R, "\n", sep="")
      message(logg) ## vypise chybovou hlasku na obrazovku
      logs <- paste(na.omit(logs), logg, sep = "//")

      L <- start.val
      CM0 <- NA
      L[["CM0"]] <- CM0
      L[["logs"]] <- logs

      return(L)
    })

    ## ErrM <- paste("Problem with nls(). Original error: ", err_R, sep="")
    ## v pripade chyby sice vypisu hlasku na obrazovku, ale jeste nevim, jak ji ulozit jako vystup funkce
    ## pripadne reseni je naznaceno v:
    ## https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function

  } else {

    lst.start.val <- start.val
    CM0 <- NA
    logs <- NA
    lst.start.val[["CM0"]] <- CM0
    lst.start.val[["logs"]] <- logs

  }

  CM0 <- lst.start.val$CM0
  logs <- lst.start.val$logs

  # Nasleduje funkce zajistujici vlastni iteraci
  lst.iteration <- OEFPILIter(y0 = y, x0 = x, L = lst.start.val, CM = CM, max.iter = max.iter,
                                 see.iter.val = see.iter.val, save.file.name = save.file.name,
                                 th = th, LOF = LOF, logs = logs)

  # CI.matrix <- CIOEFPIL(lst.iteration$lst.parameters, lst.iteration$cov_m,
  #                          signif.level = signif.level)
  ## funkce zajistujici vypocet intervalu spolehlivosti pro parametry

  contents <- list(input.form.string = input.form.string, LOF = LOF, x = x, y = y,
                   dep.var.name = dep.var.name, idp.var.name = idp.var.name,
                   names.of.parameters = names(lst.iteration$lst.parameters),
                   signif.level = signif.level)
  ## list obsahujici potrebne zbyle a uzitecne veci, ktere vyuzijeme zejmena
  ## v navazujicich funkcich

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
                         paste(names(lst.start.val)[1:l],"_start.val",sep=""), "cov.m_Est", "cov.m_nlsLM",
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
  ##   Iterace na vypocet odhadu parametru v zobecne verzi algoritmu.
  ## y0             . . . vstupni vektor y (mereni sily)
  ## x0             . . . vstupni vektor x (mereni hloubky)
  ## L              . . . list s pocatecnimi odhady paramtru
  ## CM             . . . kov. matice - tato cast se v prubehu alg. nemeni
  ## max.iter       . . . maximalni pocet iteraci algoritmu (default = 100)
  ## see.iter.val   . . . T, pokud chceme vypisovat a ukladat prubezne vysledky v prubehu iterace
  ## save.file.name . . . jmeno souboru do ktereho chce uzivtel ukladat prubeh iterace; pokud jej
  ##                      nezadame, neexportuje se nic; doporuceny format: .Rdata
  ## th             . . . treshold nutny pro zastaveni iterace
  ## LOF            . . . vstupni list funkci zadanych uzivatelem
  ## logs           . . . textovy retezec obsahujici zpravy a hlaseni
  ##                      o dosavadnim prubehu funkce

  l <- length(L) - 2
  ## pocet parametru

  L0 <- list()
  length(L0) <- l
  L0 <- L[1:l]
  ## puvodni hodnoty

  L1 <- list()
  length(L1) <- l
  ## nove hodnoty

  # hodnoty parametru se budou v prubehu menit; oznaceni "L1" je aktualni
  # hodnota; oznaceni "L0" bude hodnota par. z minuleho kroku;
  # na zacatku nastavime obe stejne
  L1 <- L0

  Q22 <- NA

  # x1 a y1 jsou hodnoty, ktere budeme v prubehu vylepsovat
  x1 <- x0
  y1 <- y0

  N <- length(y0)

  # nastaveni iteracniho pocitadla; vysledny pocet iteraci
  # tak bude vcetne nuloveho kroku
  it_num <- 0

  while((ConditionForIteration(L0, L1, th) || it_num == 0) && it_num < max.iter) {

    fval <-  sapply(x1, function(val, LP){do.call(LOF[[1]], args=c(val, LP))}, L1)
    ## Obsahuje vycisleni uzivatelem zadane funkce v hodnotach x1 pro hodnoty parametru
    ## obsazene v listu L1

    b_vec <- y1 - fval

    B11 <- - diag(sapply(x1, function(val, LP){do.call(LOF[[2]], args=c(val, LP))}, L1))
    ## Vycisleni derivace uzivatelem zadane funkce v hodnotach x1 pro hodnoty parametru
    ## obsazene v listu L1 a naskladano do diagonalni matice

    B1 <- cbind(B11, diag(N))

    if ((sum(is.infinite(B1)) != 0) || (sum(is.na(B1)) != 0)) {
      logg <- paste("There is NaN, NA, Inf or -Inf value in matrix B1. Iteration aborted!!!", "\n",
                    "The following results (if any) are only partial.", "\n", sep="")
      message(logg)
      logs <- paste(na.omit(logs), logg, sep = "//")
      break
    }
    ## preruseni while cyklu v pripade NA, NaN, Inf nebo -Inf hodnot v B1

    B2 <- CreateB2(LOF, L1, x1, l)
    ## vytvoreni matice B2

    if ((sum(is.infinite(B2)) != 0) || (sum(is.na(B2)) != 0)) {
      logg <- paste("There is NaN, NA, Inf or -Inf value in matrix B2. Iteration aborted!!!", "\n",
                    "The following results (if any) are only partial.", "\n", sep="")
      message(logg)
      logs <- paste(na.omit(logs), logg, sep = "//")
      break
    }
    ## preruseni while cyklu v pripade NA, NaN, Inf nebo -Inf hodnot v B2

    M11 <- B1 %*% CM %*% t(B1) # CM je kovariancni matice

    lst.chol <- tryCatch(expr = {
      Lmat <- t(chol(M11)) # Choleskeho rozklad matice M11

      logg <- NA

      L <- list(Lmat = Lmat, logg = logg)
      ## vystup z teto casti je posledni radek (tj. L)
    }, error = function(e) {

      err_R <- NULL
      err_R <- conditionMessage(e) # ulozeni chybove hlasky z R
      logg <- paste("Problems with computing Cholesky decomposition (chol() function). Iteration aborted!!!",
                    "\n", "Original chol error: ", err_R, "\n", sep="")
      message(logg) ## vypise chybovou hlasku na obrazovku
      Lmat <- NA

      L <- list(Lmat = Lmat, logg = logg)
      return(L)
    })
    ## Bezpecne spusteni funkce chol. Zajistime, aby funkce nespadla, ale jen vypsala
    ## varovani.

    if (grepl("Iteration aborted!!!", lst.chol$logg)) {
      logs <- paste(na.omit(logs), lst.chol$logg, sep = "//")
      break
    }
    ## Preruseni cyklu, pokud chol() hlasi chybu

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
    # Stary zpusob vypoctu pomoci ginv (i s varovnou hlaskou) ###
    ############################################################################
    #
    # M <- rbind(cbind(M11, B2), cbind(t(B2), matrix(0,l,l)))
    # ## M je matice, jejiz inverzi vznikne matice Q
    #
    # if ((sum(is.infinite(M)) > 0) | (sum(is.na(M)) > 0)) {
    #   logg <- paste("Infinite, NA or NaN input into ginv(). Iteration aborted!!!", "\n",
    #                 "The following results (if any) are only partial.", "\n", sep="")
    #   message(logg)
    #   logs <- paste(na.omit(logs), logg, sep = "//")
    #   break
    # }
    # ## preruseni while cyklu v pripade Inf, NaN nebo NA vstupu do ginv
    #
    # Q <- ginv(M)
    # ## vypocet inverzni matice
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
    ## vektor obsahujici prirustky vsech l parametru v danem kroku

    L0 <- L1

    x1 <- x1 + xdiff_hat
    y1 <- y1 + ydiff_hat
    ## vylepseni hodnot vektoru x1, y1

    L1 <- as.list(unlist(L0) + output02)
    names(L1) <- names(L0)
    ## vylepseni hodnot zbylych parametru; musime uchovat jak puvodni tak vylepsenou
    ## hodnotu kvuli podmince na zacatku cyklu

    if (IsListOK(L1) == F) {
      logg <- paste("During the iteration, the estimated parameter equals NaN, NA, Inf or -Inf.", "\n",
                    "Iteration aborted!!! The following results (if any) are only partial.", "\n", sep="")
      message(logg)
      logs <- paste(na.omit(logs), logg, sep = "//")
      break
    }
    ## kontrola korektnosti vylepsenych hodnot (L1)

    # vysledny pocet iteraci je vcetne nulteho kroku
    it_num <- it_num + 1

    # Podminka ktera vypisuje a uklada vsechny hodnoty odhadu v prubehu iterace,
    # pokud to uzivatel vyzaduje.
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
  ##   Funkce pro vypocet konfidencnich intervalu pro
  ## parametry vypoctene pomoci OEFPIL.

  cov_m <- output.form$cov.m_Est ## odhad kovariancni matice
  l <- dim(cov_m)[1] ## pocet parametru

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
  ## Funkce vypocte a vykreli bodove pasy spolehlivosti. Aplikujeme na
  ## vystupni list z funkce OEFPIL().
  ## output.form . . . vystup z OEFPIL()
  ## xx          . . . v techto bodech vykreslujeme a pocitame IS, pripadne PS

  LOF <- output.form$contents$LOF ## list of functions
  x <- output.form$contents$x ## x-ova data
  y <- output.form$contents$y ## y-ova data

  cov_m <- output.form$cov.m_Est ## odhad kovariancni matice
  l <- dim(cov_m)[1] ## pocet parametru

  lst.parameters <- output.form[1:l] ## odhad parametru
  names(lst.parameters) <- output.form$contents$names.of.parameters

  lst.parameters_previous.step <- output.form[(2*l+6):(3*l+5)]
  names(lst.parameters_previous.step) <- output.form$contents$names.of.parameters
  ## odhad paramtru z predchoziho kroku

  dep.var.name <- output.form$contents$dep.var.name ## jmeno zavisle promenne
  idp.var.name <- output.form$contents$idp.var.name ## jmeno nezavisle promenne

  if (IsListOK(lst.parameters) && IsListOK(lst.parameters_previous.step) && IsListOK(cov_m)) {

    if (missing(xx)) {
      xx <- seq(from = min(x), to = max(x), by = 0.1)
    }
    yy <- sapply(xx, function(val, LP){do.call(LOF[[1]], args=c(val, LP))}, lst.parameters)

    Omega <- sapply(1:l, function(i) {
      sapply(xx, function(val, LP){do.call(LOF[[2+i]], args=c(val, LP))}, lst.parameters_previous.step)
    })
    ## i-ty radek matice je vycisleny vektor omega v bode xx[i]

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
    ## Predpokladme, ze pribudou i dalsi vystupy, proto je format list
  }
}

################################################################################

summary.OEFPIL <- function(output.form, signif.level = output.form$contents$signif.level,
                             print = T) {
  ## Funkce slouzici k prehledenmu a rychlemu vypisu zakladnich informaci
  ## output.form  . . . vystup z OEFPIL
  ## signif.level . . . hladina vyznamnosti
  ## print        . . . udava, zdali chceme vypisovat tabulku

  if (IsListOK(output.form$cov.m_Est)) {

    l <- dim(output.form$cov.m_Est)[1] ## pocet parametru

    summary.form <- output.form[c(1:l, 2*l+1, 2*l+5, 3*l+6)]
    ## summary.form je oklesteny list, ktery obsahuje jen potrebne veci
    ## pro funkci summary
    summary.form$input.form.string <- output.form$contents$input.form.string
    ## vyber potrebnych veci ze seznamu

    if (IsListOK(summary.form)) {

      l <- dim(summary.form$cov.m_Est)[1] ## pocet parametru

      if (length(signif.level) != length(output.form$contents$signif.level)) {
        summary.form$CI_parameters <- confint.OEFPIL(output.form, signif.level)
        ## pocitame nove CI, pokud je zadana nova hl. vyznamosti

      } else if (all(signif.level != output.form$contents$signif.level)) {
        summary.form$CI_parameters <- confint.OEFPIL(output.form, signif.level)
        ## pocitame nove CI, pokud je zadana nova hl. vyznamosti

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
  ##   Funkce na "bezpecne" rozdeleni vyrazu expr (v nasem pripade formule) tak,
  ## aby ve vyslednem retezci nebyly prechody na dalsi radek

  ret <- paste(deparse(expr), collapse = "")
  ## Spojeni casti z vice radku do jedne pomoci mezer (whitespace)

  # Odstraneni "whitespace"
  gsub("[[:space:]][[:space:]]+", " ", ret)
}

################################################################################

NanoIndent <- function(data, alpha.start, m.start, hp.start, unload.data = F,
                       ucut = 0.98, lcut = 0.4, CM, uh = 0.5, uF = 0.001,
                       max.iter = 100, see.iter.val = F, save.file.name,
                       th = 0.001, signif.level = 0.05, useNLS = T) {
  ##   Funkce pro nanoindentaci. Pouziva OEFPIL, pro kterou pripravi data (osekne je),
  ## a formuli. Funkci zavislot parametru je pevne zvolena ve tvaru:
  ## f = alpha * (h - hp) ^ m, kde f je sila a h hloubka (na odtezovaci casti krivky).
  ##   Uzivatel ma moznost zadat si svoje vlastni startovaci hodnoty parametru.
  ## Pokud tak neucini, funkce ho upozorni, ze zvolila jine pocatecni hodnoty (jsou
  ## zvoleny podle navrhu v algoritmu).
  ##  Funkce NanoIndent navic oproti OEFPIL muze vypisovat zpravy navic.
  ## Napriklad kdyz "hp >= min(h)" (tato situace nedava fyzikalne smysl).
  ## unload.data . . . F, pokud mame k dispozici kompletni odtezovaci i zatezovaci krivku,
  ##                   v takovem pripade si funkce sama odtezovaci krivku vysekne;
  ##                   T, uz vime, ze mame k dispozici jen odtezovaci krivku;
  ## ucut        . . . udava horni hranici orezu F_max
  ## lcut        . . . udava nizsi hranici orezu F_max, tj. pokud ucut = 0.98,
  ##                   lcut = 0.4 uvazujeme 40 - 98 % z F_max (doporuceni normy)

  logs <- NA

  ICL <- FindIndentCurve(data, lcut = lcut, ucut = ucut, unload.data = unload.data)
  ## nemame-li k dispozici cista data, musime je vyseknout
  ## mame-li uz odtezovaci krivku, jiste oseknuti musime provest take

  data.cut <- data.frame(ICL$f_orig, ICL$h_orig)
  colnames(data.cut) <- c("f", "h")
  ## vytvoreni datoveho souboru, ktery ma pojmenovane sloupce - nutno provest
  ## pro dalsi praci s OEFPIL

  Fmax <- ICL$Fmax
  hmax <- ICL$hmax
  hmin <- ICL$hmin

  # Pokud kterakoli startovaci hodnota parametru chybi, pouzijeme vlastni
  if (missing(alpha.start) || missing(m.start) || missing(hp.start)) {

    index.miss <- which(c(missing(alpha.start), missing(m.start), missing(hp.start)))
    par.miss <- paste(c("alpha", "m", "hp")[index.miss], collapse = ", ")

    # startovaci volba pro alpha, hp, m a podle navrhu v algoritmu
    # (volba podle navrhu a korektni z fyzikalnho hlediska)
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
    ## pokud uzivatel nezada kovariancni matici, pouzije se tato

    vec.uni <- rep(1, length(data.cut$h))
    CM <- CovMatrix(uh, uF, vec.uni, vec.uni)
    ## kovariancni matice vytvorena podle navrhu v algoritmu;
    ## v prubehu se nemeni
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
  ## Nepotrebujeme, protoze data najdeme pod output.form.NI$contents$x a output.form.NI$contents$y

  return(output.form)
}

################################################################################

FindIndentCurve <- function(data, lcut = 0.4, ucut = 0.98, unload.data = F) {
  ##   Funkce slouzi k najiti pocatku a konce odtezovaci krivky a
  ## vyseknuti prislusne jeji casti. Prvni varianta je, ze na vstupu mame
  ## data ze zatezovani i odtezovani:
  ## (1) Uvazujeme za jeji pocatek ten bod, kde je namerena maximalni
  ## sila. Pokud je takovych bodu vice, uvazujeme takovy posledni
  ## namereny bod. Mereni obsahujici zaporne hodnoty sily a hloubky
  ## jsou odseknuta. Poslednim krokem je orezani krivky na procentualni
  ## rozsah z Fmax zadany pomoci lcut a ucut.
  ## (2) Mame data jen z odtezovaci krivky (unload.data = T). Dojde
  ## k odseknuti zapornych hodnot a nasledne uprava pomoci lcut a ucut.
  ## data        . . . vstupni soubor dat; 1. sloupecek je h, 2. sl. je F
  ## unload.data . . . F, pokud mame k dispozici kompletni odtezovaci i zatezovaci krivku,
  ##                   v takovem pripade si funkce sama odtezovaci krivku vysekne;
  ##                   T, uz vime, ze mame k dispozici jen odtezovaci krivku;
  ## ucut        . . . udava horni hranici orezu F_max
  ## lcut        . . . udava nizsi hranici orezu F_max, tj. pokud ucut = 0.98,
  ##                   lcut = 0.4 uvazujeme 40 - 98 % z F_max (doporuceni normy)

  if (unload.data == F) {

    # odseknuti odtezovaci krivky
    pos_Fmax <- dim(data)[1] - which.max(rev(data[,2])) + 1
    Fmax <- data[pos_Fmax,2]
    hmax <- data[pos_Fmax,1]
    data <- data[pos_Fmax:dim(data)[1],]

    # odseknuti odtezovaci krivky - stary zpusob
    #pos_hmax <- dim(data)[1] - which.max(rev(data[,1])) + 1
    #hmax <- data[pos_hmax,1]
    #Fmax <- data[pos_hmax,2]
    #data <- data[pos_hmax:dim(data)[1],]

    # odseknuti zapornych hodnot
    if (any(data < 0)) {
      mz <- min(which(data[,1] < 0), which(data[,2] < 0))
      data <- data[1:(mz - 1),]
    }

  } else {

    # odseknuti zapornych hodnot
    if (any(data < 0)) {
      mz <- min(which(data[,1] < 0), which(data[,2] < 0))
      data <- data[1:(mz - 1),]
    }

    pos_Fmax <- 1
    Fmax <- data[pos_Fmax,2]
    hmax <- data[pos_Fmax,1]

  }

  hmin <- min(data[,1])

  # vyseknuti prislusne casti odtezovaci krivky
  UB <- Fmax * ucut
  LB <- Fmax * lcut

  index <- which( data[,2] <= UB & data[,2] >= LB )

  f_orig <- data[index,2]
  h_orig <- data[index,1]

  return(list(f_orig = f_orig, h_orig = h_orig, Fmax = Fmax, hmax = hmax, hmin = hmin))
}

################################################################################

CovMatrix <- function(uh, uF, hvec, Fvec) {
  # uh, uF - nejistoty hloubky a zatizeni
  # Fvec - vektor namerenych hodnot zatizeni
  ## vystup: kovariancni matice

  N <- length(Fvec)
  M0 <- matrix(0, nrow = N, ncol = N)
  M1 <- uh^2 * diag(hvec)
  M2 <- uF^2 * diag(Fvec)^2
  covM <- cbind(rbind(M1, M0), rbind(M0, M2))
  return(covM)
}

################################################################################

ConditionForIteration <- function(L0, L1, th) {
  ## Vystupem je hodnota TRUE, pokud je podminka pro vstup do while
  ## cyklu splnena, nebo FALSE jinak.
  ## L0 . . . list obsahujici hodnoty paramtru z minuleho kroku
  ## L1 . . . list obsahujici soucasne hodnoty parametru
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
  ## Funkce na vytvoreni matice B2 v algoritmu.
  ## LOF  . . . List vsech funkci zadanych uzivatelem. Derivace podle parametru
  ##            jsou na pozicich 3 a dale ve spravnem poradi.
  ## LP   . . . list obsahujici hodnoty parametrtu ve spravnem poradi
  ## xvec . . . vektor x-ovych hodnot
  ## m    . . . pocet parametru

  B2 <- sapply(1:m, function(i) {
    sapply(xvec, function(val, LPP){do.call(LOF[[2+i]], args=c(val, LPP))}, LP)
  })

  return(-B2)
}

################################################################################

IsListOK <- function(List) {
  ## Funkce kontrolujici zda se v seznamu "List" nevyskytuji
  ## hodnoty NaN, NA, Inf or -Inf. Muze jit jak o polozky Listu,
  ## tak i jednotliva cisla v polozkach Listu.

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
  ## Funkce pridava dalsi objekty do RData souboru. Zkopirovano (a upraveno) z:
  ## https://stackoverflow.com/questions/38364745/add-a-data-frame-to-an-existing-rdata-file
  .dummy <- NULL
  if (!file.exists(rda_file)) save(.dummy, file = rda_file)

  old_e <- new.env()
  new_e <- new.env()

  load(file = rda_file, envir = old_e)

  # name_obj <- deparse(substitute(obj))   # get the name of the object
  ## Bylo v puvodni internetove verzi

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
