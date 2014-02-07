##
## simple multispecies occupancy model
##
## Created 1.17.2014 John Tipton
##

mcmcMS <- function(Y, J, n.aug, alpha.alpha.p, beta.alpha.p, alpha.beta.p, beta.beta.p, alpha.alpha.psi, alpha.beta.psi, beta.alpha.psi, beta.beta.psi, alpha.lambda, beta.lambda, alpha.p.tune, beta.p,tune, alpha.psi.tune, beta.psi.tune, n.mcmc, Z.init){ # Y is a matrix of binomial counts
	
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
#   Yaug <- cbind(Y, matrix(rep(rep(0, n), ncol = n.aug), nrow = n.aug))
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
	
	alpha.p <- rgamma(1, alpha.alpha.p, beta.alpha.p)
	beta.p <- rgamma(1, alpha.beta.p, beta.beta.p)
	alpha.psi <- rgamma(1, alpha.alpha.psi, beta.alpha.psi)
	beta.psi <- rgamma(1, alpha.beta.psi, beta.beta.psi)
	lambda <- rbeta(1, alpha.lambda, beta.lambda)
	p <- rbeta(Omega, alpha.p, beta.p)
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
	alpha.p.save <- vector(length = n.mcmc)
	beta.p.save <- vector(length = n.mcmc)
	alpha.psi.save <- vector(length = n.mcmc)
	beta.psi.save <- vector(length = n.mcmc)
	lambda.save <- vector(length = n.mcmc)
	p.accept <- 0
	psi.accept <- 0
	
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
	  ## Sample alpha.p and beta.p
  	##
		
	  alpha.p.star <- rnorm(1, alpha.p, alpha.p.tune)
	  beta.p.star <- rnorm(1, beta.p, beta.p.tune)
  	if(alpha.p.star > 0 & beta.p.star > 0){
		  mh1.p <- sum(dbeta(p, alpha.p.star, beta.p.star, log = TRUE)) + dgamma(alpha.p.star, alpha.alpha.p, beta.alpha.p, log = TRUE) + dgamma(beta.p.star, alpha.beta.p, beta.beta.p, log = TRUE)
# 		  mh1.p <- sum(dbeta(p[W1], alpha.p.star, beta.p.star, log = TRUE)) + dgamma(alpha.p.star, alpha.alpha.p, beta.alpha.p, log = TRUE) + dgamma(beta.p.star, alpha.beta.p, beta.beta.p, log = TRUE)
	  	mh2.p <- sum(dbeta(p, alpha.p, beta.p, log = TRUE)) + dgamma(alpha.p, alpha.alpha.p, beta.alpha.p, log = TRUE) + dgamma(beta.p, alpha.beta.p, beta.beta.p, log = TRUE)
# 		  mh2.p <- sum(dbeta(p[W1], alpha.p, beta.p, log = TRUE)) + dgamma(alpha.p, alpha.alpha.p, beta.alpha.p, log = TRUE) + dgamma(beta.p, alpha.beta.p, beta.beta.p, log = TRUE)
  		mh.p <- exp(mh1.p - mh2.p)
		
		  if(mh.p > runif(1)){
	  		alpha.p <- alpha.p.star
  			beta.p <- beta.p.star
			  p.accept <- p.accept + 1 / n.mcmc
		  }
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
	  alpha.p.save[l] <- alpha.p # vector(length = n.mcmc)
  	beta.p.save[l] <- beta.p # vector(length = n.mcmc)
	  alpha.psi.save[l] <- alpha.psi # vector(length = n.mcmc)
  	beta.psi.save[l] <- beta.psi # vector(length = n.mcmc)	
		lambda.save[l] <- lambda
  }
	list(p.save = p.save, psi.save = psi.save, N.save = N.save, alpha.p.save = alpha.p.save, beta.p.save = beta.p.save, alpha.psi.save = alpha.psi.save, beta.psi.save = beta.psi.save, lambda.save = lambda.save, p.accept = p.accept, psi.accept = psi.accept)
}