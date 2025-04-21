## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----step6v-------------------------------------------------------------------
library(fda)
library(puls)
# Build a simple fd object from already smoothed smoothed_arctic
data(smoothed_arctic)
NBASIS <- 300
NORDER <- 4
y <- t(as.matrix(smoothed_arctic[, -1]))
splinebasis <- create.bspline.basis(rangeval = c(1, 365),
                                    nbasis = NBASIS,
                                    norder = NORDER)
fdParobj <- fdPar(fdobj = splinebasis,
                  Lfdobj = 2,
                  # No need for any more smoothing
                  lambda = .000001)
yfd <- smooth.basis(argvals = 1:365, y = y, fdParobj = fdParobj)

## ----pulssplitv, fig.cap="Four cluster result of PULS on Arctic ice extent data for years 1979-1986, 1989-2013.", fig.height=5, fig.width = 6, echo=T, warning=F, tidy=T----
Jan <- c(1, 31); Feb <- c(31, 59); Mar <- c(59, 90)
Apr <- c(90, 120); May <- c(120, 151); Jun <- c(151, 181)
Jul <- c(181, 212); Aug <- c(212, 243); Sep <- c(243, 273)
Oct <- c(273, 304); Nov <- c(304, 334); Dec <- c(334, 365)

intervals <-
  rbind(Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec)

PULS4_pam <- PULS(toclust.fd = yfd$fd, intervals = intervals,
                  nclusters = 4, method = "pam")
plot(PULS4_pam)

## ----contrainedv, echo = T, tidy=T--------------------------------------------
constrained_PULS <- PULS(yfd$fd, intervals = intervals, nclusters = 4, 
                         spliton = 7:9, method = "pam")
print(constrained_PULS)

