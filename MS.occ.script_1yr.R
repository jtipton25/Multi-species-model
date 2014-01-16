####################

############Script for multi-species trial 1

K=150 #change to N
alpha.psi <- 1 
beta.psi <- 3
psi <- rbeta(K,alpha.psi, beta.psi) # alpha.psi = 1, beta.psi = 3
curve(dbeta(x, alpha.psi, beta.psi)

alpha.p <- 2 
beta.p <- 4
p <- rbeta(K,alpha.p, beta.p) # alpha.p = 1, beta.p = 3
curve(dbeta(x, alpha.p, beta.p))

n=250
J=5

source("sim_data_MS1.R")
output1 <- sim_data_ms1(K,psi,p,n,J)
summary(output1)

Y <- output1$Y
dim(Y)

n.mcmc <- 100

#alpha.psi <- 1
sigma.squared.alpha.psi.tune <- 0.01
alpha.alpha.psi <- 4
beta.alpha.psi <- 4
curve(dgamma(x, alpha.alpha.psi, beta.alpha.psi), from = 0, to = 2)
mean(rgamma(1000, alpha.alpha.psi, beta.alpha.psi))	

#beta.psi <- 3
sigma.squared.beta.psi.tune <- 0.01
alpha.beta.psi <- 12
beta.beta.psi <- 4
curve(dgamma(x, alpha.beta.psi, beta.beta.psi), from = 0, to = 10)
mean(rgamma(1000, alpha.beta.psi, beta.beta.psi))

#alpha.p <- 2
sigma.squared.alpha.p.tune <- 0.01
alpha.alpha.p <- 6
beta.alpha.p <- 3
curve(dgamma(x, alpha.alpha.p, beta.alpha.p), from = 0, to = 10)
mean(rgamma(1000, alpha.alpha.p, beta.alpha.p))
	
#beta.p <- 4
sigma.squared.beta.p.tune <- 0.01
alpha.beta.p <- 16
beta.beta.p <- 4
curve(dgamma(x, alpha.beta.p, beta.beta.p), from = 0, to = 10)
mean(rgamma(1000, alpha.beta.p, beta.beta.p))
	
alpha.lambda <- 5
beta.lambda <- 3

k.aug <- 500
		
source("MS.occ.mcmc_1yr.R")
output1 <- MS.occ.mcmc_1yr(Y,k.aug,alpha.psi,sigma.squared.alpha.psi.tune,alpha.alpha.psi,beta.alpha.psi,
                              beta.psi,sigma.squared.beta.psi.tune,alpha.beta.psi,beta.beta.psi,
                              alpha.p,sigma.squared.alpha.p.tune,alpha.alpha.p,beta.alpha.p,
                              beta.p,sigma.squared.beta.p.tune,alpha.beta.p,beta.beta.p,
                              alpha.lambda, beta.lambda, n.mcmc)

