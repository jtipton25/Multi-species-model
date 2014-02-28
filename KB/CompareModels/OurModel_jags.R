
# Dorazio multispeices model
library(coda)
library(R2jags)

#setwd("~/Dropbox/HootenClass_Fall2013/ZIP_and_Occupancy/MultispeciesModels")

#bbDat <- read.table("breedingBirdData.txt", header=T, sep=",")
#head(bbDat)
#butterflies <- read.table("butterflyData.txt", header=T, sep=",")
#head(butterflies)

#y <- as.matrix(bbDat)

#n <- dim(y)[1]  # number of observed species
#J <- dim(y)[2]  # number of sites
#K <- 11  # is number of surveys per site

#nzeroes <- 200
#X <- rbind(y, matrix(0, nrow=nzeroes, ncol=J) )
#Z <- ifelse(X>0, 1, NA)
#w <- c( rep(1, n), rep(NA, nzeroes) )

data <- list("X", "Z", "w", "n", "nzeroes", "J", "K")        
params <- c("lambda", "N", "psi.alpha", "psi.beta", "p.alpha", "p.beta")
params_fixed <- c("lambda", "N") ## "psi.alpha", "psi.beta", "p.alpha", "p.beta", "N")

# Initial values
inits <- function(){ list(
	lambda = n/ (n + nzeroes),
	psi.alpha = rgamma(1, 2, 1),
	psi.beta = rgamma(1, 2, 1),
	p.alpha = rgamma(1, 2, 1),
	p.beta = rgamma(1, 2, 1)
) }
inits_fixed <- function(){ list(
	lambda = n/ (n + nzeroes)
) }



OurModel <- function() {	
	lambda ~ dunif(0, 1)
	
	psi.alpha ~ dgamma(2, 1)
	psi.beta ~ dgamma(2, 1)
	p.alpha ~ dgamma(2, 1)
	p.beta ~ dgamma(2, 1)
	
	# LIKELIHOOD
	for (i in 1:(n+nzeroes)) {
		w[i] ~ dbern(lambda)
		p[i] ~  dbeta(p.alpha, p.beta)  # dbeta(1.5, 8)
		psi[i] ~ dbeta(psi.alpha, psi.beta)  #dbeta(1, 1)  # 
		mu.psi[i] <- psi[i] * w[i]
							
		for (j in 1:J) {

			Z[i, j] ~ dbern(mu.psi[i])
			mu.theta[i, j] <- p[i] * Z[i, j]
			X[i, j] ~ dbin(mu.theta[i,j], K)
		}
	}
	
	n0 <- sum( w[(n+1):(n+nzeroes)])
	N <- n + n0
}

#system.time( out2 <- jags(data, inits, params, model.file=OurModel, n.chains = 2,
#   n.thin = 2, n.iter = 40, n.burnin = 5)  )
#print(out2, 2)
#traceplot(out2)
#out.update2<- update(out2, n.iter=1000)
#print(out.update2, 2)

#traceplot(out.update2)
#curve(dbeta(x, 0.14, 0.78))

 #jagsfit.mcmc <- as.mcmc.list(out.update2)
 #   traceplot(jagsfit.mcmc)
 #   xyplot(jagsfit.mcmc)
 #   densityplot(jagsfit.mcmc)


OurModel_fixed <- function() {	
	lambda ~ dunif(0, 1)
	
	# LIKELIHOOD
	for (i in 1:(n+nzeroes)) {
		w[i] ~ dbern(lambda)
		p[i] ~  dbeta(1.5, 8)  #dbeta(p.alpha, p.beta)  # dbeta(1.5, 8)
		psi[i] ~ dbeta(1, 1)  # dbeta(psi.alpha, psi.beta)  #dbeta(1, 1)  # 
		mu.psi[i] <- psi[i] * w[i]		
				
		for (j in 1:J) {
			Z[i, j] ~ dbern(mu.psi[i])
			mu.theta[i, j] <- p[i] * Z[i, j]
			X[i, j] ~ dbin(mu.theta[i,j], K)
		}
	}
	
	n0 <- sum( w[(n+1):(n+nzeroes)])
	N <- n + n0
}
