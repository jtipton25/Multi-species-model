####################

############Script for multi-species trial 1

K=25
psi <- runif(K,0.4,0.8)
p <- runif(K,0.4,0.6)
n=250
J=5

source("sim_data_MS1.R")
output1 <- sim_data_ms1(K,psi,p,n,J)
summary(output1)

Y <- output1$Y
dim(Y)

n.mcmc <- 100

k.aug <- 50

alpha.psi <- .5
sigma.squared.alpha.psi.tune <- 0.01
alpha.alpha.psi <- 0.5
beta.alpha.psi <- 0.5

beta.psi <- .5
sigma.squared.beta.psi.tune <- 0.01
alpha.beta.psi <- 0.5
beta.beta.psi <- 0.5

alpha.p <- .5
sigma.squared.alpha.p.tune <- 0.01
alpha.alpha.p <- 0.5
beta.alpha.p <- 0.5

beta.p <- .5
sigma.squared.beta.p.tune <- 0.01
alpha.beta.p <- 0.5
beta.beta.p <- 0.5

alpha.lambda <- 2
beta.lambda <- 2

source("MS.occ.mcmc_1yr.R")
output1 <- MS.occ.mcmc_1yr(Y,k.aug,alpha.psi,sigma.squared.alpha.psi.tune,alpha.alpha.psi,beta.alpha.psi,
                              beta.psi,sigma.squared.beta.psi.tune,alpha.beta.psi,beta.beta.psi,
                              alpha.p,sigma.squared.alpha.p.tune,alpha.alpha.p,beta.alpha.p,
                              beta.p,sigma.squared.beta.p.tune,alpha.beta.p,beta.beta.p,
                              alpha.lambda, beta.lambda, n.mcmc)

