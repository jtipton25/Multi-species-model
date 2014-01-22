##
## simple multispecies occupancy model
##
## Created 1.17.2014 John Tipton
##

mcmcMS <- function(Y, n.aug, alpha.alpha.p, beta.alpha.p, alpha.beta.p, beta.beta.p, alpha.alpha.psi, alpha.beta.psi, beta.alpha.psi, beta.beta.psi, alpha.lambda, beta.lambda, alpha.p.tune, beta.p,tune, alpha.psi.tune, beta.psi.tune, n.mcmc){ # Y is a matrix of binomial counts
	
	##
	## libraries and functions
	##
	
	makeSumYZ <-function(s, Y, Z){
		sum(Y[, s] * Z[, s])
	}
	
	makeSumJYZ <-function(s, Y, Z, J){
		sum((J - Y[, s]) * Z[, s])
	}
	
	##
	## Initialize variables
	##
	
	alpha.p <- rgamma(1, alpha.alpha.p, beta.alpha.p)
	beta.p <- rgamma(1, alpha.beta.p, beta.beta.p)
	alpha.psi <- rgamma(1, alpha.alpha.psi, beta.alpha.psi)
	beta.psi <- rgamma(1, alpha.beta.psi, beta.beta.psi)
	lambda <- rbeta(1, alpha.lambda, beta.lambda)

	n <- dim(Y)[1]
	K <- dim(Y)[2]
					 
	Y.aug <- cbind(Y, matrix(rep(rep(0, n), n.aug), ncol = n.aug))
	Omega <- K + n.aug
	Y1 <- 1:K
  Y0 <- (K+1):Omega
	Y.0.list <- vector('list', length = Omega)
	for(s in 1:Omega){
		Y.0.list[[s]] <- Y.aug[, s] == 0
	}
	p <- rbeta(Omega, alpha.p, beta.p)
	psi <- rbeta(Omega, alpha.psi, beta.psi)
	W <- vector(length = Omega)
	W[ - Y0] <- 1
  W[Y0] <- rbinom(n.aug, 1, lambda)	
  Z <- matrix(nrow = n, ncol = Omega)
  Z[Y.aug > 0] <- 1 
	Z[Y.aug == 0] <- 0
	p.accept <- 0
	psi.accept <- 0
	
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
		
	for(l in 1:n.mcmc){
		if(l %% 100 == 0){
			cat(l, " ")
		} 
	
	  ##
	  ## Sample W
	  ##	
	
	  sum.Z <- apply(Z, 2, sum)	
	  lambda.tilde <- ((psi^sum.Z * (1 - psi) ^ (1 - sum.Z) * lambda) / ((psi^sum.Z * (1 - psi) ^ (1 - sum.Z) * lambda) + 1 - lambda))[Y0]
	  W[Y0] <-	rbinom(n.aug, 1, lambda.tilde)
	  W.1 <- W == 1
	  sumW <- sum(W)
	
	  ##
    ## Sample Z
	  ##
	
		##
		## work on speeding this up !!!
		##
		
	  psi.tilde <- (1 - p)^J * psi /((1 - p)^J * psi + 1 - psi) 
		for(k in 1:Omega){
			for(i in 1:n){
	  	#if(Y.aug[i, k] == 0 && W[k] == 1){
		    if(Y.0.list[[k]][i]){
				  if(W.1[k]){
   	    Z[i, k] <- rbinom(1, 1, psi.tilde[k])
				  }
		   }
		  }
	  }
		
	  sumZ <- apply(Z, 2, sum)
		
	  ## 
	  ## Sample p
	  ##
		
	  #sumYZ <- rep(0, Omega)
	  #sumYZ[1:K] <- sapply(1:K, makeSumYZ, Y = Y, Z = Z)
	  #sumJYZ <- rep(0, Omega)
	  #sumJYZ[1:K] <- sapply(1:K, makeSumJYZ, Y = Y, Z = Z, J = J)
		
	  #for(k in 1:K){
		#  p[k] <- rbeta(1, alpha.p + sumYZ[k], beta.p + sumJYZ[k])
	  #}
		
		sumYZ <- sapply(1:K, makeSumYZ, Y = Y, Z = Z)
		sumJYZ <- sapply(1:K, makeSumJYZ, Y = Y, Z = Z, J = J)
		p[1:K] <- rbeta(1, alpha.p + sumYZ, beta.p + sumJYZ)
		p[Y0] <- rbeta(n.aug, alpha.p, beta.p)
		
	  ##
	  ## Sample Psi
  	##
	
#	  for(k in 1:Omega){
#	  	if(W.1[k]){
#    	psi[k] <- rbeta(1, alpha.psi + sumZ[k], beta.psi + n - sumZ[k])
#		  }
#	  }	
		psi[W.1] <- rbeta(1, alpha.psi + sumZ, beta.psi + n - sumZ)

	  ##
	  ## Sample alpha.p and beta.p
  	##
		
	  alpha.p.star <- rnorm(1, alpha.p, alpha.p.tune)
	  beta.p.star <- rnorm(1, beta.p, beta.p.tune)
  	if(alpha.p.star > 0 && beta.p.star > 0){
		  mh1.p <- sum(dbeta(p, alpha.p.star, beta.p.star, log = TRUE)) + dgamma(alpha.p.star, alpha.alpha.p, beta.alpha.p, log = TRUE) + dgamma(beta.p.star, alpha.beta.p, beta.beta.p, log = TRUE)
	  	mh2.p <- sum(dbeta(p, alpha.p, beta.p, log = TRUE)) + dgamma(alpha.p, alpha.alpha.p, beta.alpha.p, log = TRUE) + dgamma(beta.p, alpha.beta.p, beta.beta.p, log = TRUE)
  		mh.p <- exp(mh1.p - mh2.p)
		
		  if(mh.p > runif(1)){
	  		alpha.p <- alpha.p.star
  			beta.p <- beta.p.star
			  p.accept <- p.accept + 1 / n.mcmc
		  }
	  }
		
	  ##
	  ## Sample alpha.p and beta.p
  	##
		
	  alpha.psi.star <- rnorm(1, alpha.psi, alpha.psi.tune)
	  beta.psi.star <- rnorm(1, beta.psi, beta.psi.tune)
  	if(alpha.psi.star > 0 && beta.psi.star > 0){
		  mh1.psi <- sum(dbeta(p, alpha.psi.star, beta.psi.star, log = TRUE)) + dgamma(alpha.psi.star, alpha.alpha.psi, beta.alpha.psi, log = TRUE) + dgamma(beta.psi.star, alpha.beta.psi, beta.beta.psi, log = TRUE)
	  	mh2.psi <- sum(dbeta(p, alpha.psi, beta.psi, log = TRUE)) + dgamma(alpha.psi, alpha.alpha.psi, beta.alpha.psi, log = TRUE) + dgamma(beta.psi, alpha.beta.psi, beta.beta.psi, log = TRUE)
  		mh.psi <- exp(mh1.psi - mh2.psi)
			
		  if(mh.psi > runif(1)){
	  		alpha.psi <- alpha.psi.star
  			beta.psi <- beta.psi.star
			  psi.accept <- psi.accept + 1 / n.mcmc
		  }
	  }
	
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
	list(p.save = p.save, psi.save = psi.save, N.save = N.save, alpha.p.save = alpha.p.save, beta.p.save = beta.p.save, alpha.psi.save = alpha.psi.save, beta.psi.save = beta.psi.save, lambda = lambda, p.accept = p.accept, psi.accept = psi.accept)
}