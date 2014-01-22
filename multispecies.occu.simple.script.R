
# Multi-species Occupancy model: Simulation and model test.


setwd("/Users/broms/Dropbox/HootenClass_Fall2013/ZIP_and_Occupancy/MultispeciesModels")

###
### Simulate the data.
###

J <- 4    # number of surveys
N <- 200  # number of species
n <- 25  # number of sites

Omega = 300  # number of "possible species"
#lambda <- N / Omega  # probability of being a "real" species

curve( dbeta(x, 2, 4) )
curve( dbeta(x, 1, 3) )

alpha.p <- 1.5
beta.p <- 8
alpha.psi <- 3
beta.psi <- 3
curve( dbeta(x, 1.5,8) )
curve( dbeta(x, 3, 3) )

# ***** Assuming that detection and occupancy probabilities are independent. ******
det.prob <- rbeta(N, alpha.p, beta.p)
occu.prob <- rbeta(N, alpha.psi, beta.psi)

z <- matrix(0, nrow=n, ncol=Omega)
tmp.y <- z
for(k in 1:N){
  z[, k] <- rbinom(n, 1, occu.prob[k])
  tmp.y[, k] <- rbinom(n, J, det.prob[k]) * z[, k]
}
(realizedN <- length( which(apply(z, 2, sum) > 0) ) )  # N = 200 of the 200 possible birds actually occupy a site
y <- tmp.y[ , which(apply(tmp.y, 2, sum) > 0)]
dim(y)   # 181 birds ever observed

###
### Run the model on the data.
###
tune <- 0.2
source("multispp.occu.simple.mcmc.R")
out <- multispp.occu.simple.mcmc(y, J, n.mcmc=200,
			alpha.psi.tune=tune, beta.psi.tune=tune,
			alpha.p.tune=tune, beta.p.tune=tune)
outB <- more.mcmc(out, 200)
#
#  With vaguer priors:
tune=0.5
source("multispp.occu.simple.mcmc.R")
out2 <- multispp.occu.simple.mcmc(y, J, n.mcmc=2000, Omega=250,
 prior.shape.alpha.p=1, prior.rate.alpha.p=0.1,
 prior.shape.beta.p=1, prior.rate.beta.p=0.1,
 prior.shape.alpha.psi=1, prior.rate.alpha.psi=0.1,
 prior.shape.beta.psi=1, prior.rate.beta.psi=0.1,
			alpha.psi.tune=tune, beta.psi.tune=tune,
			alpha.p.tune=tune, beta.p.tune=tune)

out2B <- more.mcmc(out2, 3000)
#




###
### test on real data
### BBS data from Chapter 12 of R and D Book
###


#  They estimates a spp richness of N = 138 (SD=19)
getwd()
dat <- read.csv("MultiSppData_RandDbook_Ch12.csv")
head(dat)
real.y <- t(dat[, -1])
dim(real.y)  # n=50 sites, k=99 observed species.
J <- 11  # 11 days of sampling.
Omega <- 500

tune=0.5
source("multispp.occu.simple.mcmc.R")
out2 <- multispp.occu.simple.mcmc(real.y, J, n.mcmc=10000,
 prior.shape.alpha.p=1, prior.rate.alpha.p=0.01,
 prior.shape.beta.p=1, prior.rate.beta.p=0.01,
 prior.shape.alpha.psi=1, prior.rate.alpha.psi=0.01,
 prior.shape.beta.psi=1, prior.rate.beta.psi=0.01,
			alpha.psi.tune=tune, beta.psi.tune=tune,
			alpha.p.tune=tune, beta.p.tune=tune)
			
out2B <- more.mcmc(out2, 5000)



lambda.tilde.num <- (out2$lambda * out2$psi ^ apply(out2$z, 2, sum, na.rm=T) * (1-out2$psi) ^ (apply(1-out2$z, 2, sum, na.rm=T)) )
summary(lambda.tilde.num)

###
### test on real data- Swiss BBS
###

# analysis described Chapter 12 of R and D Book
#  N = 170 (151-195), while K = 134
getwd()
dat <- read.csv("")
head(dat)
# n=254 sites, k=99 observed species.
J <- 3  # 2 or 3 sampling occasions
## Detection requires a date!?!
Omega <- 500

tune=0.2
source("multispp.occu.simple.mcmc.R")
out3 <- multispp.occu.simple.mcmc(real.y, J, n.mcmc=200,
 prior.shape.alpha.p=1, prior.rate.alpha.p=0.1,
 prior.shape.beta.p=1, prior.rate.beta.p=0.1,
 prior.shape.alpha.psi=1, prior.rate.alpha.psi=0.1,
 prior.shape.beta.psi=1, prior.rate.beta.psi=0.1,
			alpha.psi.tune=tune, beta.psi.tune=tune,
			alpha.p.tune=tune, beta.p.tune=tune)
out3B <- more.mcmc(out3, 5000)
