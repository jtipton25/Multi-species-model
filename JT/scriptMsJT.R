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

n <- 50 #25
N <- 500 #300
J <- 2 #4

# presence probability
alpha.psi <- 5
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

hist(data$p, freq = FALSE)
curve(dbeta(x, alpha.p, beta.p), add = TRUE, col = 'red')

hist(data$psi, freq = FALSE)
curve(dbeta(x, alpha.psi, beta.psi), add = TRUE, col = 'red')


##
## Set up hyperpriors knowing the true values
##

#alpha.psi <- 5
sigma.squared.alpha.psi.tune <- 0.01
alpha.alpha.psi <- 10
beta.alpha.psi <- 2
curve(dgamma(x, alpha.alpha.psi, beta.alpha.psi), from = 0, to = 10)
abline(v = alpha.psi, col = 'red')

#beta.psi <- 6
sigma.squared.beta.psi.tune <- 0.01
alpha.beta.psi <- 12
beta.beta.psi <- 2
curve(dgamma(x, alpha.beta.psi, beta.beta.psi), from = 0, to = 10)
abline(v = beta.psi, col = 'red')

#alpha.p <- 1.5
sigma.squared.alpha.p.tune <- 0.01
alpha.alpha.p <- 6
beta.alpha.p <- 3
curve(dgamma(x, alpha.alpha.p, beta.alpha.p), from = 0, to = 10)
abline(v = alpha.p, col = 'red')

#beta.p <- 8
sigma.squared.beta.p.tune <- 0.01
alpha.beta.p <- 16
beta.beta.p <- 2
curve(dgamma(x, alpha.beta.p, beta.beta.p), from = 0, to = 12)
abline(v = beta.p, col = 'red')


n.aug <- 250

dim(data$Y)[2]
alpha.lambda <- 1.5
beta.lambda <- 5
curve(dbeta(x, alpha.lambda, beta.lambda))
abline(v = (N - dim(data$Y)[2]) / n.aug, col = 'red')
##
## Initialize MCMC parameters
##

n.mcmc <- 200
n.burn <- floor(n.mcmc / 5)
alpha.p.tune <- 0.15
beta.p.tune <- 0.15
alpha.psi.tune <- 0.15
beta.psi.tune <- 0.15

start <- Sys.time()
out <- mcmcMS(data$Y, n.aug, alpha.alpha.p, beta.alpha.p, alpha.beta.p, beta.beta.p, alpha.alpha.psi, alpha.beta.psi, beta.alpha.psi, beta.beta.psi, alpha.lambda, beta.lambda, alpha.p.tune, beta.p,tune, alpha.psi.tune, beta.psi.tune, n.mcmc)
finish <- Sys.time() - start
finish

dim(data$Y)[2]
N

make.plot(out)
names(out)

