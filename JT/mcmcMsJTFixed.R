set.seed(1)
setwd('~/Multi-species-model/JT/')

##
## Libraries and functions
##

source('make.plot.R')
source('simMsJT.R')
source('mcmcMsJTNew.R')

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

n.aug <- 250

dim(data$Y)[2]
alpha.lambda <- 1.5
beta.lambda <- 5
curve(dbeta(x, alpha.lambda, beta.lambda))
abline(v = (N - dim(data$Y)[2]) / n.aug, col = 'red')
##
## Initialize MCMC parameters
##

n.mcmc <- 20000
n.burn <- floor(n.mcmc / 5)
alpha.p.tune <- 0.05
beta.p.tune <- 0.05
alpha.psi.tune <- 0.15
beta.psi.tune <- 0.15
Z.init <- 0.0075

start <- Sys.time()
out <- mcmcMS(data$Y, n.aug, alpha.p, beta.p, alpha.psi, beta.psi, alpha.lambda, beta.lambda, n.mcmc, Z.init)
finish <- Sys.time() - start
finish

dim(data$Y)[2]
N
plot(out$N.save, type = 'l')
n.burn <- floor(n.mcmc / 5) + 1
layout(matrix(1:4, nrow = 2))
plot(out$N.save[n.burn:n.mcmc], type = 'l')
hist(out$N.save[n.burn:n.mcmc])
abline(v = N, col = 'red')

hist(apply(out$p.save[, n.burn:n.mcmc], 1, mean), freq = FALSE)
curve(dbeta(x, alpha.p, beta.p), add = TRUE) ## note the lack of identifiability
hist(apply(out$psi.save[, n.burn:n.mcmc], 1, mean), freq = FALSE)
curve(dbeta(x, alpha.psi, beta.psi), add = TRUE) ## note the lack of identifiability

