set.seed(1)
setwd('~/Multi-species-model/JT/Link function/')

##
## Libraries and functions
##

source('makePlotLink.R')
source('simMsJTLink.R')
source('mcmcMsJTLink.R')
source('rMVN.R')
library(mvtnorm)

##
## Initialize simulation parameters
##

n <- 50 #25
N <- 500 #300
J <- 2 #4

## detection probability
beta <- c(-3, 3)
tau.beta <- length(beta)
U <- matrix(c(rnorm(n, 0, 1), rnorm(n, 0, 1)), ncol = n, byrow = TRUE)

## occurrence probability
eta <- c(3, 1)
tau.eta <- length(eta)
V <- matrix(c(rnorm(n, 0, 1), rnorm(n, 0, 1)), ncol = n, byrow = TRUE)

##
## Sample data
##

data <- makeMultiSpec(n, N, J, beta, U, eta, V)

hist(data$p, freq = FALSE)
hist(data$psi, freq = FALSE)


##
## Set up hyperpriors knowing the true values
##


n.aug <- 250

dim(data$Y)[2]
alpha.lambda <- 1.5
beta.lambda <- 5
curve(dbeta(x, alpha.lambda, beta.lambda))
abline(v = (N - dim(data$Y)[2]) / n.aug, col = 'red')

##
## Initialize MCMC parameters
##

n.mcmc <- 5000
n.burn <- floor(n.mcmc / 5)
mu.beta <- rep(0, tau.beta)
sigma.squared.beta <- 100
mu.eta <- rep(0, tau.eta)
sigma.squared.eta <- 100
Z.init <- 0.015

start <- Sys.time()
out <- mcmcMS(data$Y, n.aug, mu.beta, sigma.squared.beta, mu.eta, sigma.squared.eta, alpha.lambda, beta.lambda, beta.tune, eta.tune, n.mcmc, Z.init)
finish <- Sys.time() - start
finish

dim(data$Y)[2]
N

make.plot(out)
names(out)

