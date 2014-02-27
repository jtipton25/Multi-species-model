
# Multi-species Occupancy model: Simulation and model test.


setwd("~/Dropbox/HootenClass_Fall2013/ZIP_and_Occupancy/MultispeciesModels")

###
### Simulate the data.
###

K <- 3    # number of surveys; previously called J
N <- 300  # true number of species
J <- 25  # number of sites; previously called n

Omega = 400  # number of "possible species"
#lambda <- N / Omega  # probability of being a "real" species

curve( dbeta(x, 2, 4) )
curve( dbeta(x, 1, 3) )

ALPHA.p <- 1.5
BETA.p <-  8
ALPHA.psi <- 0.4  #3
BETA.psi <- 0.4  #3
curve( dbeta(x, 1.5, 8) )
curve( dbeta(x, 0.4, 0.4) )

# ***** Assuming that detection and occupancy probabilities are independent. ******
det.prob <- rbeta(N, ALPHA.p, BETA.p)
occu.prob <- rbeta(N, ALPHA.psi, BETA.psi)

z <- matrix(NA, nrow=Omega, ncol=J)
tmp.y <- z
for(k in 1:N){
  z[k, ] <- rbinom(J, 1, occu.prob[k])
  tmp.y[k, ] <- rbinom(J, K, det.prob[k]) * z[k, ]
}
(realizedN <- length( which(apply(z, 1, sum) > 0) ) )  # N = 200 of the 200 possible birds actually occupy a site
y <- tmp.y[which(apply(tmp.y, 1, sum) > 0), ]
dim(y)   # 274 of 300 birds ever observed
n <- dim(y)[1]
in.dat <- y

###
### Run the model on the data.
###
tune <- 0.1
source("multispp.occu.simple.mcmc_LessIndexing.R")
out <- multispp.occu.simple.mcmc(y, K, N=N, N.MCMC=5000,
			alpha.psi.tune=tune, beta.psi.tune=tune,
			alpha.p.tune=tune, beta.p.tune=tune)
#outB <- more.mcmc(out, 1000)


#
#  With vaguer priors:
tune=0.8
tune.psi=0.2
source("multispp.occu.simple.mcmc_LessIndexing.R")
out2 <- multispp.occu.simple.mcmc(y, K,  N=N, N.MCMC=2000, 
  alpha.lambda=1, beta.lambda=1,
 prior.shape.alpha.p=1, prior.rate.alpha.p=0.1,
 prior.shape.beta.p=1, prior.rate.beta.p=0.1,
 prior.shape.alpha.psi=1, prior.rate.alpha.psi=0.1,
 prior.shape.beta.psi=1, prior.rate.beta.psi=0.1,
			alpha.psi.tune=tune.psi, beta.psi.tune=tune.psi,
			alpha.p.tune=tune, beta.p.tune=tune)

out2B <- more.mcmc(out2, 3000)
#



###
### test on real data- Dorazio et al. Ecology paper
###

bbDat <- read.table("breedingBirdData.txt", header=T, sep=",")
head(bbDat)  # they estimate 93 species total
butterflies <- read.table("butterflyData.txt", header=T, sep=",")
head(butterflies)  # they estimate 28 species total

in.dat <- butterflies  # bbDat
y <- as.matrix(in.dat)
n <- dim(in.dat)[1]
J <- dim(in.dat)[2]
K <- 11
nzeroes <- 100

tune <- 0.3
tune.psi <- 0.1
source("multispp.occu.simple.mcmc_LessIndexing.R")
out <- multispp.occu.simple.mcmc(y, K, N=93, N.MCMC=500,
			alpha.psi.tune=tune, beta.psi.tune=tune,
			alpha.p.tune=tune, beta.p.tune=tune)
out2 <- multispp.occu.simple.mcmc(y, K,  N=93, N.MCMC=5000, 
  alpha.lambda=1, beta.lambda=1,
 prior.shape.alpha.p=2, prior.rate.alpha.p=0.5,
 prior.shape.beta.p=2, prior.rate.beta.p=0.5,
 prior.shape.alpha.psi=2, prior.rate.alpha.psi=0.5,
 prior.shape.beta.psi=2, prior.rate.beta.psi=0.5,
			alpha.psi.tune=tune.psi, beta.psi.tune=tune.psi,
			alpha.p.tune=tune, beta.p.tune=tune)



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
