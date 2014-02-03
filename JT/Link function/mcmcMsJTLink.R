##
## simple multispecies occupancy model
##
## Created 1.17.2014 John Tipton
##

mcmcMS <- function(Y, U, V, n.aug, mu.beta, sigma.squared.beta, mu.eta, sigma.squared.eta, alpha.lambda, beta.lambda, beta.tune, eta.tune, n.mcmc, Z.init){ # Y is a matrix of binomial counts
	
	##
	## libraries and functions
	##
	
# 	makeSumYZ <-function(s, Yaug, Z){
# 		Yaug[, s] %*% Z[, s]
# 	}
# 	
# 	makeSumJYZ <-function(s, Yaug, Z, J){
# 	  ((J - Yaug[, s]) %*% Z[, s])
# 	}
	
	##
  ## set up for MCMC
	##
  
  n <- dim(Y)[1]
	K <- dim(Y)[2]
	## Augment the data matrix				 
	Yaug <- cbind(Y, matrix(rep(rep(0, n), n.aug), ncol = n.aug))
	Omega <- K + n.aug
	## setup indexes for MCMC loop
  Y0 <- apply(Yaug, 2, sum) == 0
  Y1 <- !Y0	
#   (K+1):Omega
  Yaug0 <- Yaug == 0
  nYaug0 <- !Yaug0

  ##
	## Initialize variables
	##
	tau.beta <- length(mu.beta)
  beta <- rmvnorm(1, mu.beta, sigma.squared.beta * diag(tau.beta))
  tau.eta <- length(mu.eta)  
  eta <- rnmnorm(1, mut.eta, sigma.squared.eta * diag(tau.eta))
  lambda <- rbeta(1, alpha.lambda, beta.lambda)
	p <- expit(U)
	psi <- rbeta(Omega, alpha.psi, beta.psi)
	W <- vector(length = Omega)
	W[Y1] <- 1
	W[Y0] <- rbinom(n.aug, 1, lambda)	
	Z <- matrix(nrow = n, ncol = Omega)
	Z[nYaug0] <- 1 
	Z[Yaug0] <- rbinom(sum(Yaug0), 1, Z.init)
	sumZ <- apply(Z, 2, sum)  
	
	##
	## Save variables
	##
	
	p.save <- matrix(nrow = Omega, ncol = n.mcmc)
	psi.save <- matrix(nrow = Omega, ncol = n.mcmc)
	N.save <- vector(length = n.mcmc)
	beta.save <- vector(length = n.mcmc)
	eta.save <- vector(length = n.mcmc)
	lambda.save <- vector(length = n.mcmc)
	beta.accept <- 0
	eta.accept <- 0
	
	for(l in 1:n.mcmc){
		if(l %% 100 == 0){
			cat(l, " ")
		} 
	
		##
		## Sample W
		##  
		
		lambda.tilde <- rep(((1 - psi)^n * lambda) / ((1 - psi)^n * lambda + (1 - lambda)), n)
		W[Y0 & (sumZ > 0)] <- 1
		W[Y0 & (sumZ == 0)] <- rbinom(sum(Y0 & (sumZ == 0)), 1, lambda.tilde[Y0 & (sumZ == 0)])
		W1 <- W == 1
		sumW <- sum(W)
		W.mat <- matrix(rep(W1, n), nrow = n, byrow = TRUE)
		    
		##
		## Sample Z
		##
		
		psi.tilde <- matrix(rep(((1 - p)^J * psi) / ((1 - p)^J * psi + (1 - psi)), n), nrow = n, ncol = Omega, byrow = TRUE)[Yaug0 & W.mat]
		Z[Yaug0 & !W.mat] <- 0
		Z[Yaug0 & W.mat] <- rbinom(sum(Yaug0 & W.mat), 1, psi.tilde)
		#     psi.tilde <- matrix(rep(((1 - p)^J * psi) / ((1 - p)^J * psi + (1 - psi)), n), nrow = n, ncol = Omega, byrow = TRUE)[Yaug0 & !W.mat]
		#     Z[Yaug0 & W.mat] <- 0
		#     Z[Yaug0 & !W.mat] <- rbinom(sum(Yaug0 & !W.mat), 1, psi.tilde)
		sumZ <- apply(Z, 2, sum)
		
		## 
		## Sample p
		##
		
		sumYZ <- apply(Yaug * Z, 2, sum)
		sumJYZ <- apply((J - Yaug) * Z, 2, sum)
		p <- rbeta(Omega, alpha.p + sumYZ, beta.p + sumJYZ) 

		##
		## Sample Psi
		##
		
		psi <- rbeta(Omega, alpha.psi + sumZ * W, beta.psi + (n - sumZ) * W)
		
    ##
	  ## Sample beta
  	##
		
    
    ## need to think about how to use covariates
	  beta.star <- rnorm(2, beta, beta.tune)
      p.star <- expit(beta.star %*% U)
		  mh1.beta <- sum(dbinom(Yaug, J, p.star, log = TRUE)) + dgamma(alpha.p.star, alpha.alpha.p, beta.alpha.p, log = TRUE) + dgamma(beta.p.star, alpha.beta.p, beta.beta.p, log = TRUE)
	  	mh2.p <- sum(dbeta(p, alpha.p, beta.p, log = TRUE)) + dgamma(alpha.p, alpha.alpha.p, beta.alpha.p, log = TRUE) + dgamma(beta.p, alpha.beta.p, beta.beta.p, log = TRUE)
  		mh.p <- exp(mh1.p - mh2.p)
		
		  if(mh.p > runif(1)){
	  		alpha.p <- alpha.p.star
  			beta.p <- beta.p.star
			  p.accept <- p.accept + 1 / n.mcmc
		  }
    rm(alpha.p.star)
    rm(beta.p.star)	

	  ##
	  ## Sample alpha.psi and beta.psi
  	##
		
	  alpha.psi.star <- rnorm(1, alpha.psi, alpha.psi.tune)
	  beta.psi.star <- rnorm(1, beta.psi, beta.psi.tune)
  	if(alpha.psi.star > 0 & beta.psi.star > 0){
		  mh1.psi <- sum(dbeta(psi, alpha.psi.star, beta.psi.star, log = TRUE)) + dgamma(alpha.psi.star, alpha.alpha.psi, beta.alpha.psi, log = TRUE) + dgamma(beta.psi.star, alpha.beta.psi, beta.beta.psi, log = TRUE)
#   	  mh1.psi <- sum(dbeta(psi[W1], alpha.psi.star, beta.psi.star, log = TRUE)) + dgamma(alpha.psi.star, alpha.alpha.psi, beta.alpha.psi, log = TRUE) + dgamma(beta.psi.star, alpha.beta.psi, beta.beta.psi, log = TRUE)
  	  mh2.psi <- sum(dbeta(psi, alpha.psi, beta.psi, log = TRUE)) + dgamma(alpha.psi, alpha.alpha.psi, beta.alpha.psi, log = TRUE) + dgamma(beta.psi, alpha.beta.psi, beta.beta.psi, log = TRUE)
#       mh2.psi <- sum(dbeta(psi[W1], alpha.psi, beta.psi, log = TRUE)) + dgamma(alpha.psi, alpha.alpha.psi, beta.alpha.psi, log = TRUE) + dgamma(beta.psi, alpha.beta.psi, beta.beta.psi, log = TRUE)
  		mh.psi <- exp(mh1.psi - mh2.psi)
      
		  if(mh.psi > runif(1)){
	  		alpha.psi <- alpha.psi.star
  			beta.psi <- beta.psi.star
			  psi.accept <- psi.accept + 1 / n.mcmc
		  }
	  }
    rm(alpha.psi.star)
    rm(beta.psi.star)
	
  	##
    ## Sample lambda	
    ##				
							
    lambda <- rbeta(1, alpha.lambda + sumW, beta.lambda + Omega - sumW)
		
	  ##
  	## Save samples
	  ##

  	p.save[, l] <- p # matrix(nrow = Omega, ncol = n.mcmc)
	  psi.save[, l] <- psi # matrix(nrow = Omega, ncol = n.mcmc)
  	N.save[l] <- sumW # vector(length = n.mcmc)
  	beta.save[l] <- beta # vector(length = n.mcmc)
  	eta.save[l] <- eta # vector(length = n.mcmc)	
		lambda.save[l] <- lambda
  }
	list(p.save = p.save, psi.save = psi.save, N.save = N.save, beta.save = beta.save, eta.save = eta.save, lambda.save = lambda.save, beta.accept = beta.accept, eta.accept = eta.accept, Z.save = Z.save)
}
