##################################################################
### Inicializace startovacich hodnot
##################################################################
source("TACR_Rsource.R")

num <- 28
file_name <- paste("BV52_09_07_2019.txt_", num, ".unload.txt", sep="")

# 1, 2, 4, 5 ## problemove soubory pro data nize; nikde by ale algoritmus nemel havarovat
#file_name <- paste("./Data/UNHT-Cu_unload/kalibrace-hloubky-2014-nocontact-FD_", num, "_mn.unload.txt", sep="")

data <- read.delim(file = file_name, sep="", header = F, col.names = c("y", "x"))

id = num
unload.data = T
ucut = 0.98
lcut = 0.2
uh = 0.5
uF = 0.001
max.iter = 100
see.iter.val = F
th = 0.001
signif.level = 0.05

useNLS = T
# save.file.name = "Export_iter.Rdata"
## jmeno souboru zadavame jen v pripade, pokud chceme ulozit prubezna data z iterace



ICL <- FindIndentCurve(data, lcut = lcut, ucut = ucut, unload.data = unload.data)
## nemame-li k dispozici cista data, musime je vyseknout
## mame-li uz odtezovaci krivku, jiste oseknuti musime provest take

data.cut <- data.frame(ICL$f_orig, ICL$h_orig)
colnames(data.cut) <- c("y", "x")
## vytvoreni datoveho souboru, ktery ma pojmenovane sloupce - nutno provest

### kovarianci matice ###
vec.uni <- rep(1, length(data.cut$x))
CM <- CovMatrix(uh, uF, vec.uni, vec.uni)
## kovariancni matice vytvorena podle navrhu v algoritmu;
## v prubehu se nemeni

###########################################################################
### Vstupni list startovacich hodnot pro OEFPIL ###
###########################################################################

# f_orig <- ICL$f_orig # je v data.cut
# h_orig <- ICL$h_orig # je v data.cut
Fmax <- ICL$Fmax
hmax <- ICL$hmax
hmin <- ICL$hmin

# startovaci volba pro alpha, hp, m a podle navrhu v algoritmu
# (volba podle navrhu a korektni z fyzikalnho hlediska)
m.start <- 1.5
hp.start <- 0.9 * hmin
alpha.start <- Fmax / ((hmax - hp.start) ^ m.start)
start.val <- list(alpha=alpha.start, m=m.start, hp=hp.start)
# start.val <- c(alpha.start, m.start, hp.start)
names(start.val) <- c("alpha", "m", "hp")
form <- yy ~ alpha * (xx - hp) ^ m
###########################################################################
# OEFPIL ###
###########################################################################
typeof(data.cut)
data.cut1 <- as.matrix(data.cut)
colnames(data.cut1) <- c('[,1]','[,2]')

matrix(c(1,2,3,4), nrow = 2)

data.cut2 <- as.data.frame(data.cut)
names(data.cut2) <- c("yy","xx")

output.form <- OEFPIL(data.cut2, form, start.val, CM,  max.iter = max.iter,
                         see.iter.val = see.iter.val,
                         th = th, signif.level = signif.level, useNLS = F)
a <- summary(output.form)
output.form$alpha_start.val

xx <- seq(from = min(output.form$contents$x), to = max(output.form$contents$x), by = 0.1)
## vykreslovaci sekvence

CB.list <- plot(output.form, signif.level = signif.level,
                               xx = xx, add.plot = T)

################################################################################
# Priprava na vkladani dat primo do funkce OEFPIL ###
# data <- data.cut
# save.file.name = "Export_iter.Rdata"


# Priprava na vkladani dat primo do funkce OEFPILIter
# y0 <- y; x0 <- x; L <- lst.start.val
#  Konec pripravy dat ###
################################################################################


################################################################################
### NanoIndent ###
################################################################################
output.form.NI <- NanoIndent(data, unload.data = unload.data, ucut = ucut, lcut = lcut,
                             uh = uh, uF = uF, max.iter = max.iter,
                             see.iter.val = see.iter.val, th = th,
                             signif.level = signif.level, useNLS = useNLS)
summary(output.form.NI)

# Zadani NanoIndent s vlastnimi startovacimi hodnotami parametru
output.form.NI <- NanoIndent(data, alpha.start = alpha.start, m.start = m.start,
                             hp.start = hp.start, unload.data = unload.data, ucut = ucut, lcut = lcut,
                             uh = uh, uF = uF, max.iter = max.iter,
                             see.iter.val = see.iter.val, th = th,
                             signif.level = signif.level, useNLS = useNLS)
summary(output.form.NI)

xx <- seq(from = min(output.form.NI$contents$x), to = max(output.form.NI$contents$x), by = 0.1)
## vykreslovaci sekvence

CB.list <- plot(output.form.NI, signif.level = signif.level,
                               xx = xx, add.plot = T)
################################################################################


### vstupni list startovacich hodnot -- odpovida souboru 28 v BV52 ###
# alpha.start <- 0.002474752
# m.start <- 1.5
# hp.start <- 65.98112
#
# start.val <- list(alpha=alpha.start, m=m.start, hp=hp.start)
# start.val <- c(alpha.start, m.start, hp.start)
# names(start.val) <- c("alpha", "m", "hp")


# formule definovana uzivatelem
#form <- y ~ sin(exp(alpha*(x - hp)^m) + 1) + 5
# form <- y ~ alpha*(x - hp)^m

# Zkusebni zadani 1:
# form <- y ~ a * x * pi + b * sin (a ^ 2)
#  start.val <- list(b = 0.1, c = 1, a = 1, pi = 45)
#  names(start.val) <- c("b", "c", "a", "pi")

# Zkusebni zadani 2:
#form <- y ~ 24.5*x^2 + c1*x + c2*x^0.5 + c3*x^0.25 + c4*x^0.125 + c5*x^0.0625
#start.val <- c(1:5)
#names(start.val) <- c("c1", "c2", "c3", "c4", "c5")

# data.cut = spravna data se kterymi bychom meli pracovat
# data.p = spatna data na kterych zkousime podminky ve fci
#data.p <- data.cut
#data.p[,3] <- data.cut[,1]

##################################################################
### Kontrola vysledku se starsi verzi algoritmu
##################################################################

expected_m = c(0,2)
see_est = F

ReturnList <- ParEst(data, id = num, unload.data = unload.data, ucut = ucut, lcut = lcut,
                     CM, uh = uh, uF = uF, max_iter = max.iter,
                     expected_m = expected_m, see_est = see_est, th = th, signif_level = signif.level,
                     useNLS = useNLS )


##################################################################
### Testovani problemoveho souboru
##################################################################
# algoritmus nezkonverguje a v listu odhadu parametru se objevuji NA hodnoty
num <- 10
file_name <- paste("./DATA/Si_unload/WG80530_fddata.TXT_", num, ".unload.txt", sep="")
data <- read.delim(file = file_name, sep="", header = F, col.names = c("y", "x"))
##################################################################
### Pracovne-testovaci sekce
##################################################################

A1 <- list(a = 1, b = 2, c = 3)
A2 <- list(a = 1.2, b = 2.2, c = 3.2)
A3 <- list(a = 1, b = 2.2, c = 3.3, d = 4)

save(A1, file = "ulozeneListy.Rdata")
save(A2, file = "ulozeneListy.rda")
AddObjectToRdata(A3, "ulozeneListy.Rdata", overwrite = FALSE)

##################################################################
### Test vykresleni
##################################################################
# vstupni hodnoty
data.h <- data$h
data.f <- data$f
xx <- seq(from = min(data.h), to = max(data.h), by = 0.1)
## vybrane h podle nasi funkce
alpha_previous.step <- ReturnList$alpha_previous.step
m_previous.step <- ReturnList$m_previous.step
hp_previous.step <- ReturnList$hp_previous.step
cov_m <- ReturnList$cov_m
signif.level <- 0.05

PointwiseCB_kresli(data.h=data.h, data.f=data.f, xx=xx,
            alpha_previous.step=alpha_previous.step, m_previous.step=m_previous.step,
            hp_previous.step=hp_previous.step, cov_m=cov_m, signif.level=signif.level)

xx <- seq(from = min(data$h), to = max(data$h), length = 500)
yy <- ReturnList[[num]]$alpha * (xx - ReturnList[[num]]$hp) ^ ReturnList[[num]]$m
yy0 <- ReturnList[[num]]$alpha0 * (xx - ReturnList[[num]]$hp0) ^ ReturnList[[num]]$m0
lines(xx, yy, lwd = 2, col = "blue")
lines(xx, yy0, lwd = 2, col = "red")
