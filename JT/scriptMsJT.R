set.seed(1)
setwd('~/Multi-species-model/JT/')

##
## Libraries and functions
##

source('make.plot.R')
source('simMsJT.R')
source('mcmcMsJT.R')

##
## Initialize simulation parameters
##

n <- 25
N <- 300
J <- 4

# presence probability
alpha.psi <- 6
beta.psi <- 6
psi <- rbeta(N, alpha.psi, beta.psi) # alpha.psi = 1, beta.psi = 3
curve(dbeta(x, alpha.psi, beta.psi))

# detection probability
alpha.p <- 1.5 
beta.p <- 8
p <- rbeta(N,alpha.p, beta.p) # alpha.p = 1, beta.p = 3
curve(dbeta(x, alpha.p, beta.p))

##
## Sample data
##

data <- makeMultiSpec(n, N, J, alpha.psi, beta.psi, alpha.p, beta.p)

##
## Set up hyperpriors knowing the true values
##

#alpha.psi <- 1
sigma.squared.alpha.psi.tune <- 0.01
alpha.alpha.psi <- 4
beta.alpha.psi <- 4
curve(dgamma(x, alpha.alpha.psi, beta.alpha.psi), from = 0, to = 2)
abline(v = alpha.psi, col = 'red')

#beta.psi <- 3
sigma.squared.beta.psi.tune <- 0.01
alpha.beta.psi <- 12
beta.beta.psi <- 4
curve(dgamma(x, alpha.beta.psi, beta.beta.psi), from = 0, to = 10)
abline(v = beta.psi, col = 'red')

#alpha.p <- 2
sigma.squared.alpha.p.tune <- 0.01
alpha.alpha.p <- 6
beta.alpha.p <- 3
curve(dgamma(x, alpha.alpha.p, beta.alpha.p), from = 0, to = 10)
abline(v = alpha.p, col = 'red')

#beta.p <- 4
sigma.squared.beta.p.tune <- 0.01
alpha.beta.p <- 16
beta.beta.p <- 4
curve(dgamma(x, alpha.beta.p, beta.beta.p), from = 0, to = 10)
abline(v = beta.p, col = 'red')

alpha.lambda <- 2
beta.lambda <- 5
curve(dbeta(x, alpha.lambda, beta.lambda))

##
## Initialize MCMC parameters
##

n.aug <- 1000
n.mcmc <- 4000
alpha.p.tune <- 0.05
beta.p.tune <- 0.05
alpha.psi.tune <- 0.05
beta.psi.tune <- 0.05

start <- Sys.time()
out <- mcmcMS(data$Y, n.aug, alpha.alpha.p, beta.alpha.p, alpha.beta.p, beta.beta.p, alpha.alpha.psi, alpha.beta.psi, beta.alpha.psi, beta.beta.psi, alpha.lambda, beta.lambda, alpha.p.tune, beta.p,tune, alpha.psi.tune, beta.psi.tune, n.mcmc)
finish <- Sys.time() - start
finish

dim(data$Y)[2]

make.plot(out)