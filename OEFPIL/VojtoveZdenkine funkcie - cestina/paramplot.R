library(ggplot2)
paramplot <- function(output.list){
  # input - list of 'oefpil' objects
  # output - ggplot graf odhadu parametru with erorr bars (plus minus standard deviation)
  
  if(!is.list(output.list)){
    stop("Input has to be a list of 'oefpil' objects.")
  } #kontrola vstupu
  
  if(!IsListOK(output.list)){
    stop("There are NA or NaN in estimated parameter values or in the covariance matrix.")
  } #kontrola NaN hodnot
  
  # vytvoreni datove struktury vhodne pro ggplot
  N <- length(output.list)
  
  if (is.vector(output.list[[1]])){
    coefnames <- output.list$contents$names.of.parameters
    k <- length(coefnames)
    data <- matrix(NA, nrow = k, ncol = 4)
    data <- as.data.frame(data)
    data[,1] <- rep(1, each = k)
    data[,2] <- coefnames
    data[,3] <- unlist(output.list[1:k])
    data[,4] <- sqrt(diag(output.list$cov.m_Est))
  } else {
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

  ggplot(data, aes(x = model, y = est, col = model)) +
    geom_pointrange(aes(ymin = est - sdest, ymax =  est + sdest))+
    labs(x = "", y = "") +
    facet_wrap(~ cf, scale = "free") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
}

## examples paramplot()
library(MASS)
steamdata <- steam
colnames(steamdata) <- c("x","y")
n <- nrow(steamdata)
CM1 <- diag(rep(0.1,2*n))
CM2 <- diag(c(rep(0.2^2,n), rep(0.1^2,n)))
st1 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 5, b2 = 8, b3 = 200),
              CM1, useNLS = F)
st2 <- OEFPIL(steamdata, y ~ b1 * 10^(b2 * x/ (b3 + x)), list(b1 = 5, b2 = 8, b3 = 200),
              CM2, useNLS = F)
paramplot(list(st1,st2))

# pripadne by asi bylo fajn uvest i priklad pro jeden samostatny model, kdy mi to proste vyhodi
# parametry toho konkretniho modelu s erorr bary
paramplot(st1)

# pokud bychom chteli explicitne uvadet i priklad, jak to dopadne pri pouziti vstupu, kde mam
# v oefpil NaN hodnoty, protoze algoritmus nezkonvergoval, tak tam jde pouzit ten priklad 
# se spatnymi start. parametry (b1 = 0.1), kdy to pak hlasi chybu, ale nemyslim, ze je to
# uplne nutne 
