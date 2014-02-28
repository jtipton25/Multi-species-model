
# Dorazio multispeices model
??jags
library(coda)
library(R2jags)

?bugs2jags

getwd()
# setwd("/Users/broms/Dropbox/HootenClass_Fall2013/ZIP_and_Occupancy/MultispeciesModels")
setwd('~/Multi-species-model/KB/')

bbDat <- read.table("breedingBirdData.txt", header=T, sep=",")
head(bbDat)
butterflies <- read.table("butterflyData.txt", header=T, sep=",")
head(butterflies)

in.dat <- bbDat


n <- dim(in.dat)[1]
J <- dim(in.dat)[2]
K <- 3  #11
nzeroes <- 100
X <- rbind(as.matrix(in.dat), matrix(0, nrow=nzeroes, ncol=J) )
Z <- ifelse(X>0, 1, NA)
w <- c( rep(1, n), rep(NA, nzeroes) )
data <- list("X", "Z", "w", "n", "nzeroes", "J", "K")  
      
params <- c("omega", "psi.mean", "theta.mean", "sigma.u",
  "sigma.v", "rho", "N")

# Initial values
inits <- function(){ list(
	omega = runif(1),
	psi.mean = runif(1),
	theta.mean = runif(1),
	tau.u = rgamma(1, 0.1, 0.1),
	tau.v = rgamma(1, 0.1, 0.1), 
	rho = runif(1, -1, 1)
) }


DorazioModel <- function() {	
	omega ~ dunif(0, 1)
	
	psi.mean ~ dunif(0, 1)
	alpha <- log(psi.mean) - log(1 - psi.mean)
	
	theta.mean ~ dunif(0, 1)
	beta <- log(theta.mean) - log(1 - theta.mean)
	
	tau.u ~ dgamma(0.1, 0.1)
	tau.v ~ dgamma(0.1, 0.1)
	rho ~ dunif(-1, 1)
	var.eta <- tau.v / (1 - pow(rho, 2) )
	
	sigma.u <- 1 / sqrt(tau.u)
	sigma.v <- 1 / sqrt(tau.v)
	
	# LIKELIHOOD
	for (i in 1:(n+nzeroes)) {
		w[i] ~ dbern(omega)
		phi[i] ~ dnorm(alpha, tau.u)
		
		mu.eta[i] <- beta + (rho * sigma.v/sigma.u) * (psi[i] - alpha)
		eta[i] ~ dnorm(mu.eta[i], var.eta)
		
		logit(psi[i]) <- phi[i]
		logit(theta[i]) <- eta[i]
		
		mu.psi[i] <- psi[i]*w[i]
		for (j in 1:J) {
			Z[i, j] ~ dbern(mu.psi[i])
			mu.theta[i, j] <- theta[i] * Z[i, j]
			X[i, j] ~ dbin(mu.theta[i,j], K)
		}
	}
	
	n0 <- sum( w[(n+1):(n+nzeroes)])
	N <- n + n0
}

system.time( out <- jags(data, inits, params, model.file=DorazioModel, n.chains = 2,
   n.thin = 1, n.iter = 40, n.burnin = 5)  )


system.time( out <- jags(data, inits, params, model.file=DorazioModel, n.chains = 2,
   n.thin = 2, n.iter = 2000, n.burnin = 100)  )
print(out, 2)
traceplot(out)

out.update <- update(out, n.iter=10000)
print(out.update, 2)



###
### what patterns of probabilities are produced with the Dorazio model?
###

library(boot)

sigma.u <- 1.2a
( p.bar <- 0.8 ) # runif(1, 0, 1) )
( alpha <- logit(p.bar) )
( sigma.u <- runif(1, 0.5, 5) )
u <- rnorm(1000, 0, sigma.u)
hist( inv.logit(u + alpha), breaks=30 )