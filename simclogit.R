library(drgee)

expit <- function(x) {
    return (1/(1+exp(-x)))
}

logit <- function(p) {
    return (log(p/(1-p)))
}

beta1 <- 1.5
beta2 <- 2
beta3 <- -0.5
## number of clusters
m <- 100000
## cluster sizes
m.size <- rpois(m,1.5)+3
## m.size <- rep(5, m)
lopnr <- rep(1:m,times=m.size)
n <- length(lopnr)
l1 <- rnorm(n, sd = 1)
l2 <- rnorm(n, sd = 1)
l3 <- rnorm(n, sd = 1)
beta0 <- rep(rnorm(m, mean=0, sd=0.5), m.size)
y <- rbinom(n, 1, p = expit(beta0 + l1 * beta1 + l2 * beta2 ) )

simdata <- data.frame(lopnr, y, l1, l2, l3)

setwd("c:/Users/johzet/Rwd/drgee")
write.table(simdata, "simdata.dat")


setwd("c:/Users/johzet/Rwd/drgee")
simdata <- read.table("simdata.dat")

####################################

library(drgee)

library(tictoc)

tic()
geefit <- gee(formula = y ~ l1 + l2, link = "logit",
              clusterid = "lopnr", data = simdata, cond = TRUE)
toc()

vcov(geefit)
## geefit <- gee(formula = y ~ l1 + l2 + l3, link = "logit", clusterid = "lopnr", cond = TRUE)

Rprof(gee(formula = y ~ l1 + l2, link = "logit",
              clusterid = "lopnr", data = simdata, cond = TRUE) )


tablex <- rbind(
    geefit$time0,
    geefit$time1,
    geefit$time2,
    geefit$time3)

rownames(tablex) <- c("clogit", "Cpp-res", "R-res all", "R-res disc")
tablex

library(xtable)

xtable(tablex[, 1:3], digits = 3)

## coef(clfit)
## coef(geefit)

## vcov(clfit)
## vcov(geefit)

#######################################################################################
#######################################################################################
#######################################################################################

## DR-estimation when
## the main model is
## E(Y|A,L1,L2)-E(Y|A=0,L1,L2)=beta0*A+beta1*A*L1
## and the outcome nuisance model is
## E(Y|A=0,L1,L2)=gamma0+gamma1*L1+gamma2*L2
## and the exposure nuisance model is
## E(A|Y=0,L1,L2)=expit(alpha0+alpha1*L1+alpha2*l2)

library(drgee)

expit<-function(x) exp(x)/(1+exp(x))

n<-5000

## nuisance
l1<-rnorm(n, mean = 0, sd = 1)
l2<-rnorm(n, mean = 0, sd = 1)

beta0<-1.5
beta1<-1
gamma0<--1
gamma1<--2
gamma2<-2
alpha0<-1
alpha1<-5
alpha2<-3

## Exposure generated from the exposure nuisance model
a<-rbinom(n,1,expit(alpha0 + alpha1*l1 + alpha2*l2))
## Outcome generated from the main model and the
## outcome nuisance model
y<-rnorm(n,
mean = beta0 * a + beta1 * a * l1 + gamma0 + gamma1 * l1 + gamma2 * l2,
sd = 1)

simdata<-data.frame(y,a,l1,l2)

## outcome nuisance model misspecified and
## exposure nuisance model correctly specified

## DR-estimation
dr.est <- drgee(oformula = formula(y~l1),
eformula = formula(a~l1+l2),
iaformula = formula(~l1),
olink = "identity", elink = "logit",
data = simdata, estimation.method = "dr")
summary(dr.est)

## O-estimation
o.est <- drgee(exposure = "a", oformula = formula(y~l1),
iaformula = formula(~l1), olink = "identity", data = simdata,
estimation.method = "o")
summary(o.est)

## E-estimation
e.est <- drgee(outcome = "y", eformula = formula(a~l1+l2),
iaformula = formula(~l1), elink="logit", data = simdata,
estimation.method = "e")
summary(e.est)

#######################################################################################
#######################################################################################
#######################################################################################

library(data.table)

beta1 <- 1.5
## number of clusters
m <- 5
## cluster sizes
m.size <- rep(4, m)
lopnr <- rep(1:m,times=m.size)
n <- length(lopnr)
l1 <- rnorm(n, sd = 1)
beta0 <- rep(rnorm(m, mean=0, sd=0.5), m.size)
y0 <- rbinom(n, 1, p = expit(beta0) )
y1 <- rbinom(n, 1, p = expit(beta0 + beta1 * l1) )

dt1 <- as.data.table(list(lopnr, t(cbind(y0, y1))))
dt1

setkey(dt1, lopnr)

agg_sums <- dt1[, j = list( sum(y0), sum(y1)), by = lopnr]
agg_sums <- dt1[, lapply(.SD, sum), by = lopnr]
agg_sums
as.matrix(agg_sums)
as.matrix( agg_sums )[, -1]
