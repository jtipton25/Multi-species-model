multispp.occu.simple.mcmc <- function(y=y, J=J, Omega=500,
									  alpha.lambda=3, beta.lambda=5, 
                                      prior.shape.alpha.p=6, prior.rate.alpha.p=3,
                                      prior.shape.beta.p=16, prior.rate.beta.p=4,
                                      prior.shape.alpha.psi=12, prior.rate.alpha.psi=4,
                                      prior.shape.beta.psi=4, prior.rate.beta.psi=4,
                                      alpha.psi.tune=0.1, beta.psi.tune=0.1,
                                      alpha.p.tune=0.1, beta.p.tune=0.1,
                                      drawPlots=T, n.mcmc=1000){
  
  #
  #  Kristin Broms (20140103), Last Modified:  20140103
  #  Multispecies Ocucpancy Model
  #  Model goals: 1.  estimate spp richness.  
  #     2.  Determing detection and occu prob. distributions
  #  Model assumes data has already been augmented with 0s.
  
  #####
  #####  Libraries and Subroutines
  #####
  
  library(pscl)  # for inverse gamma functions
  
  #####
  #####  Setup Variables 
  #####
  
  n <- dim(y)[[1]]
  psi <- rep(NA, Omega)
  p <- rep(NA, Omega)
  mat0 <- matrix(rep(0, n*(Omega-dim(y)[[2]]) ), nrow=n)
  Y <- cbind(y, mat0)  # augmented data matrix
  psi.acc <- 0
  p.acc <- 0
  
  psi.alpha.save <- rep(0, n.mcmc)
  psi.beta.save <- rep(0, n.mcmc)
  p.alpha.save <- rep(0, n.mcmc)
  p.beta.save <- rep(0, n.mcmc)
  N.save <- rep(0, n.mcmc)
  
  n.burn <- n.mcmc/10
  
  #####
  #####  Priors and Starting Values 
  #####

  z <- ifelse(Y == 0, NA, 1)

  spp.onsite <-  (apply(z, 2, sum, na.rm=T) > 0)
  K <- dim(y)[[2]]
  w <- c(rep(1, K), rep(0, Omega - K) )
  
  alpha.psi <- rgamma(1, prior.shape.alpha.psi, prior.rate.alpha.psi)
  beta.psi <- rgamma(1, prior.shape.beta.psi, prior.rate.beta.psi)
  alpha.p <- rgamma(1, prior.shape.alpha.p, prior.rate.alpha.p)
  beta.p <- rgamma(1, prior.shape.beta.p, prior.rate.beta.p)
  
  psi[(w==1)] <- rbeta(sum(w==1), alpha.psi, beta.psi)
  p[(spp.onsite)] <- rbeta(sum(spp.onsite), alpha.p, beta.p)
  lambda <- rbeta(1, alpha.lambda, beta.lambda)

  
  
  #####
  #####  Begin MCMC Loop 
  #####  
  for(n.mcmc in 1:n.mcmc){
    if(n.mcmc %% 100 == 0) cat(n.mcmc, " "); flush.console()
    

 	#####
	##### Sample w_k's
	#####
	
	y.k <- apply(Y, 2, sum)  # was a species detected at any site?
	lambda.tilde.num <- (lambda * psi ^ apply(z, 2, sum, na.rm=T) * (1-psi) ^ (apply(1-z, 2, sum, na.rm=T)) )
	lambda.tilde <- lambda.tilde.num / (lambda.tilde.num + 1 - lambda)
	w[y.k==0]  <- rbinom(sum(y.k==0), 1, lambda.tilde[y.k==0])


    #####
    #####  Sample z 
    #####
 
    psi.tilde <- ( psi * (1-p) ^ J ) / ( psi * (1-p) ^ J + 1-psi )     
    for(k in 1:Omega){  # Loops are sooooooo sloooooooooooooooow!
      if(w[k]==1){
      	for(i in 1:n){
      		if(Y[i,k] == 0)   z[i, k] <- rbinom( 1, 1, psi.tilde[k] )     		
 #   		z[i, w[k]==1 && Y[i, k]==0] <- rbinom( sum(w[k]==1 && Y[i, k]==0), 1, psi.tilde[w[k]==1 & Y[i, k]==0] )     		
      	}
      }
    }
	  
    #####
    ##### Sample p
    #####
    
    p[(w==1)] <- rbeta(sum(w==1), apply(Y*z, 2, sum, na.rm=T) + alpha.p, apply(z*(J-Y), 2, sum, na.rm=T) + beta.p)
 
 
    #####
    ##### Sample psi
    #####
    
    psi[w==1] <- rbeta(sum(w==1), apply(z, 2, sum, na.rm=T) + alpha.psi, apply(1-z, 2, sum, na.rm=T) + beta.psi)
        
     
    #####
    #####  Sample p parameters
    #####

    alpha.p.star <- rnorm(1, alpha.p, alpha.p.tune)
    beta.p.star <- rnorm(1, beta.p, beta.p.tune)
    
    if(alpha.p.star > 0 && beta.p.star > 0){
      mh1.p <- sum( dbeta(p[w==1], alpha.p.star, beta.p.star, log=TRUE) ) + 
        log(dgamma(alpha.p.star, prior.shape.alpha.p, prior.rate.alpha.p) ) + 
        log(dgamma(beta.p.star, prior.shape.beta.p, prior.rate.beta.p) )
      mh2.p <- sum( dbeta(p[w==1], alpha.p, beta.p, log=TRUE) ) + 
        log(dgamma(alpha.p, prior.shape.alpha.p, prior.rate.alpha.p) ) + 
        log(dgamma(beta.p, prior.shape.beta.p, prior.rate.beta.p) )

      mh.ratio.p <- exp( mh1.p - mh2.p ) 
#      cat("\n", alpha.p, "\t", alpha.p.star, "\t the betas:  ", beta.p, "\t", beta.p.star)
#      cat("\n", mh1.p, "\t", mh2.p, "\t", mh.ratio.p)
      #print(mh.ratio.p)  
      
      if(mh.ratio.p > runif(1)){      
        alpha.p <- alpha.p.star      
        beta.p <- beta.p.star      
        p.acc <- p.acc + 1    
      }       
    }
 
    #####
    #####  Sample psi parameters
    #####
    
    alpha.psi.star <- rnorm(1, alpha.psi, alpha.psi.tune)
    beta.psi.star <- rnorm(1, beta.psi, beta.psi.tune)
    
    if(alpha.psi.star > 0 & beta.psi.star > 0){
      mh1.psi <- sum( dbeta(psi[w==1], alpha.psi.star, beta.psi.star, log=TRUE) ) + 
        log(dgamma(alpha.psi.star, prior.shape.alpha.psi, prior.rate.alpha.psi) ) + 
        log(dgamma(beta.psi.star, prior.shape.beta.psi, prior.rate.beta.psi) )
      mh2.psi <- sum( dbeta(psi[w==1], alpha.psi, beta.psi, log=TRUE) ) + 
        log(dgamma(alpha.psi, prior.shape.alpha.psi, prior.rate.alpha.psi) ) + 
        log(dgamma(beta.psi, prior.shape.beta.psi, prior.rate.beta.psi) )
   
      mh.ratio.psi <- exp( mh1.psi - mh2.psi )    
 #     cat("\n", mh1.psi, "\t", mh2.psi, "\t", mh.ratio.psi)
      #print(mh.ratio.psi)  

      if(mh.ratio.psi > runif(1)){      
        alpha.psi <- alpha.psi.star      
        beta.psi <- beta.psi.star      
        psi.acc <- psi.acc + 1    
      }       
    } 
 
 
    #####
    #####  Sample lambda
    #####
   
   lambda <- rbeta(1, alpha.lambda + sum(w), beta.lambda + sum(1-w) )
   
   
    #####  
    #####  Save Samples   
    #####
 
#    cat("\n", alpha.p, "\t", alpha.p.star, "\t the betas:  ", beta.p, "\t", beta.p.star)
#    cat("\n", alpha.psi, "\t", alpha.psi.star, "\t the betas:  ", beta.psi, "\t", beta.psi.star)
   
    psi.alpha.save[n.mcmc] <- alpha.psi  
    psi.beta.save[n.mcmc] <- beta.psi  
    p.alpha.save[n.mcmc] <- alpha.p  
    p.beta.save[n.mcmc] <- beta.p  
    N.save[n.mcmc] <- sum(w)

  }
  cat("\n")



  if(drawPlots==T){
	####
	####  Make Plots
	####
    
	layout( matrix(c(1:4, 5, 5, 6, 6), 4, 2 ) )
	plot( p.alpha.save[-(1:n.burn)], type="l", lty=1, main="Detection prob. Alpha" )
	abline(h=2)
	plot( p.beta.save[-(1:n.burn)], type="l", lty=1 , main="Detection prob. Beta" )
	abline(h=4)
	plot( psi.alpha.save[-(1:n.burn)], type="l", lty=1, main="Occu Prob alpha" )
	abline(h=1)
	plot( psi.beta.save[-(1:n.burn)], type="l", lty=1, main="Occu Prob beta")
    abline(h=3)
	
  	curve(dbeta(x, mean(p.alpha.save[-(1:n.burn)]), mean(p.beta.save[-(1:n.burn)]) ) )
  	curve(dbeta(x, mean(psi.alpha.save[-(1:n.burn)]), mean(psi.beta.save[-(1:n.burn)]) ) )
    }	
	
	#####
	##### Print estimates and 95% CIs to screen
	#####

    cat( "\n Detection prob.:  ")
    cat("\t", round(mean(p.alpha.save[-(1:n.burn)]), 2), "\t", 
              round(mean(p.beta.save[-(1:n.burn)]), 2 ) )

    cat( "\n Occupancy prob.:  ")
    cat("\t", round(mean(psi.alpha.save[-(1:n.burn)]), 2), "\t", 
              round(mean(psi.beta.save[-(1:n.burn)]), 2), "\n")

     cat( "\n As medians, Detection prob.:  ")
    cat("\t", round(median(p.alpha.save[-(1:n.burn)]), 2), "\t", 
              round(median(p.beta.save[-(1:n.burn)]), 2) )

    cat( "\n As medians, Occupancy prob.:  ")
    cat("\t", round(median(psi.alpha.save[-(1:n.burn)]), 2), "\t", 
              round(median(psi.beta.save[-(1:n.burn)]), 2), "\n" )
    
    cat("\n Species richness:  ", quantile(N.save[-(1:n.burn)], c(0.025, 0.5, 0.975) ) )



  #####
  #####  Write Output 
  #####
  
  cat("\n Psi acceptance rate: \t", psi.acc/n.mcmc)
  cat("\n Det. prob acceptance rate: \t", p.acc/n.mcmc, "\n")

  list(psi.alpha.save=psi.alpha.save, psi.beta.save=psi.beta.save,
     p.alpha.save=p.alpha.save, p.beta.save=p.beta.save,
     N.save=N.save, n.mcmc=n.mcmc, n.burn=n.burn, psi.acc=psi.acc, 
     p=p, psi=psi, z=z, Y=Y, w=w, J=J, n=n, Omega=Omega,
     alpha.lambda=alpha.lambda, beta.lambda=beta.lambda, lambda=lambda,
     prior.shape.alpha.p=prior.shape.alpha.p, prior.rate.alpha.p=prior.rate.alpha.p,
     prior.shape.beta.p=prior.shape.beta.p, prior.rate.beta.p=prior.rate.beta.p,
     prior.shape.alpha.psi=prior.shape.alpha.psi, prior.rate.alpha.psi=prior.rate.alpha.psi,
     prior.shape.beta.psi=prior.shape.beta.psi, prior.rate.beta.psi=prior.rate.beta.psi,
     alpha.psi.tune=alpha.psi.tune, beta.psi.tune=beta.psi.tune,
     alpha.p.tune=alpha.p.tune, beta.p.tune=beta.p.tune, 
     lambda.tilde=lambda.tilde)
 
}



more.mcmc <- function(out, n.mcmc, drawPlots=T){
  
  #  Run more iterations after a model has been initiated
    
  #####
  #####  Setup Variables 
  #####
  
  n <- dim(out$Y)[[1]]
  psi.acc <- 0
  p.acc <- 0
  
  psi.alpha.save <- rep(0, n.mcmc)
  psi.beta.save <- rep(0, n.mcmc)
  p.alpha.save <- rep(0, n.mcmc)
  p.beta.save <- rep(0, n.mcmc)
  N.save <- rep(0, n.mcmc)
  
  
  #####
  #####  Priors and Starting Values 
  #####

  psi <- out$psi
  p <- out$p
  Y <- out$Y
  z <- out$z
  w <- out$w
  J <- out$J
  Omega <- out$Omega
  n <- out$n
  lambda <- out$lambda
  
  alpha.psi <- tail(out$psi.alpha.save, 1)
  beta.psi <- tail(out$psi.alpha.save, 1)
  alpha.p <- tail(out$psi.alpha.save, 1)
  beta.p <- tail(out$psi.alpha.save, 1)
  
    
  
  #####
  #####  Begin MCMC Loop 
  #####  
  for(n.mcmc in 1:n.mcmc){
    if(n.mcmc %% 100 == 0) cat(n.mcmc, " "); flush.console()
    
    #####
    #####  Sample z 
    #####
 
    psi.tilde <- ( psi * (1-p) ^ J ) / ( psi * (1-p) ^ J + 1-psi )     
    for(k in 1:Omega){  # Loops are sooooooo sloooooooooooooooow!
      if(w[k]==1){
      	for(i in 1:n){
      		if(Y[i,k] == 0)   z[i, k] <- rbinom( 1, 1, psi.tilde[k] )     		    		
      	}
      }
    }

 	#####
	##### Sample w_k's
	#####
	
	y.k <- apply(Y, 2, sum)  # was a species detected at any site?
	lambda.tilde.num <- (lambda * psi ^ apply(z, 2, sum, na.rm=T) * (1-psi) ^ (apply(1-z, 2, sum, na.rm=T)) )
	lambda.tilde <- lambda.tilde.num / (lambda.tilde.num + 1 - lambda)
	w[y.k==0]  <- rbinom(sum(y.k==0), 1, lambda.tilde[y.k==0])

	  
    #####
    ##### Sample p
    #####
    
    p[(w==1)] <- rbeta(sum(w==1), apply(Y*z, 2, sum, na.rm=T) + alpha.p, 
    								   apply(z*(J-Y), 2, sum, na.rm=T) + beta.p)
 
 
    #####
    ##### Sample psi
    #####
    
    psi[w==1] <- rbeta(sum(w==1), apply(z, 2, sum, na.rm=T) + alpha.psi, 
    								   apply(1-z, 2, sum, na.rm=T) + beta.psi)
        
     
    #####
    #####  Sample p parameters
    #####

    alpha.p.star <- rnorm(1, alpha.p, out$alpha.p.tune)
    beta.p.star <- rnorm(1, beta.p, out$beta.p.tune)
    
    if(alpha.p.star > 0 && beta.p.star > 0){
      mh1.p <- sum( dbeta(p[w==1], alpha.p.star, beta.p.star, log=TRUE) ) + 
        log(dgamma(alpha.p.star, out$prior.shape.alpha.p, out$prior.rate.alpha.p) ) + 
        log(dgamma(beta.p.star, out$prior.shape.beta.p, out$prior.rate.beta.p) )
      mh2.p <- sum( dbeta(p[w==1], alpha.p, beta.p, log=TRUE) ) + 
        log(dgamma(alpha.p, out$prior.shape.alpha.p, out$prior.rate.alpha.p) ) + 
        log(dgamma(beta.p, out$prior.shape.beta.p, out$prior.rate.beta.p) )

      mh.ratio.p <- exp( mh1.p - mh2.p )       
      if(mh.ratio.p > runif(1)){      
        alpha.p <- alpha.p.star      
        beta.p <- beta.p.star      
        p.acc <- p.acc + 1    
      }       
    }
 
    #####
    #####  Sample psi parameters
    #####
    
    alpha.psi.star <- rnorm(1, alpha.psi, out$alpha.psi.tune)
    beta.psi.star <- rnorm(1, beta.psi, out$beta.psi.tune)
    
    if(alpha.psi.star > 0 & beta.psi.star > 0){
      mh1.psi <- sum( dbeta(psi[w==1], alpha.psi.star, beta.psi.star, log=TRUE) ) + 
        log(dgamma(alpha.psi.star, out$prior.shape.alpha.psi, out$prior.rate.alpha.psi) ) + 
        log(dgamma(beta.psi.star, out$prior.shape.beta.psi, out$prior.rate.beta.psi) )
      mh2.psi <- sum( dbeta(psi[w==1], alpha.psi, beta.psi, log=TRUE) ) + 
        log(dgamma(alpha.psi, out$prior.shape.alpha.psi, out$prior.rate.alpha.psi) ) + 
        log(dgamma(beta.psi, out$prior.shape.beta.psi, out$prior.rate.beta.psi) )
   
      mh.ratio.psi <- exp( mh1.psi - mh2.psi )    
      if(mh.ratio.psi > runif(1)){      
        alpha.psi <- alpha.psi.star      
        beta.psi <- beta.psi.star      
        psi.acc <- psi.acc + 1    
      }       
    } 
 
 
    #####
    #####  Sample lambda
    #####
   
    lambda <- rbeta(1, out$alpha.lambda + sum(w), out$beta.lambda + sum(1-w) )
   
      
    #####  
    #####  Save Samples   
    #####
   
    psi.alpha.save[n.mcmc] <- alpha.psi  
    psi.beta.save[n.mcmc] <- beta.psi  
    p.alpha.save[n.mcmc] <- alpha.p  
    p.beta.save[n.mcmc] <- beta.p  
    N.save[n.mcmc] <- sum(w)

  }
  cat("\n")


  if(drawPlots==T){
	####
	####  Make Plots
	####
    
	layout( matrix(c(1:4, 5, 5, 6, 6), 4, 2 ) )
	plot( p.alpha.save, type="l", lty=1, main="Detection prob. Alpha")
	abline(h=2)
	plot( p.beta.save, type="l", lty=1 , main="Detection prob. Beta" )
	abline(h=4)
	plot( psi.alpha.save, type="l", lty=1, main="Occu Prob alpha")
	abline(h=1)
	plot( psi.beta.save, type="l", lty=1, main="Occu Prob beta")
    abline(h=3)
	
  	curve(dbeta(x, mean(p.alpha.save), mean(p.beta.save) ) )
  	curve(dbeta(x, mean(psi.alpha.save), mean(psi.beta.save) ) )
    }	
	
	#####
	##### Print estimates and 95% CIs to screen
	#####

    cat( "\n Detection prob.:  ")
    cat("\t", round(mean(p.alpha.save), 2), "\t", 
              round(mean(p.beta.save), 2 ) )
    cat( "\n Occupancy prob.:  ")
    cat("\t", round(mean(psi.alpha.save), 2), "\t", 
              round(mean(psi.beta.save), 2), "\n" )
    cat("\n Species richness:  ", quantile(N.save, c(0.025, 0.5, 0.975) ) )


  #####
  #####  Write Output 
  #####
  
  cat("\n Psi acceptance rate: \t", psi.acc/n.mcmc)
  cat("\n Det. prob acceptance rate: \t", p.acc/n.mcmc, "\n")

  list(psi.alpha.save=c(out$psi.alpha.save, psi.alpha.save), 
       psi.beta.save=c(out$psi.beta.save, psi.beta.save),
       p.alpha.save=c(out$p.alpha.save, p.alpha.save), 
       p.beta.save=c(out$p.beta.save, p.beta.save),
       N.save=c(out$N.save, N.save), 
       n.mcmc=n.mcmc, n.burn=out$n.burn, psi.acc=psi.acc, 
     p=p, psi=psi, z=z, Y=Y, w=w, J=J, n=n, Omega=Omega,
     alpha.lambda=alpha.lambda, beta.lambda=beta.lambda, lambda=lambda,
     prior.shape.alpha.p=prior.shape.alpha.p, prior.rate.alpha.p=prior.rate.alpha.p,
     prior.shape.beta.p=prior.shape.beta.p, prior.rate.beta.p=prior.rate.beta.p,
     prior.shape.alpha.psi=prior.shape.alpha.psi, prior.rate.alpha.psi=prior.rate.alpha.psi,
     prior.shape.beta.psi=prior.shape.beta.psi, prior.rate.beta.psi=prior.rate.beta.psi,
     alpha.psi.tune=alpha.psi.tune, beta.psi.tune=beta.psi.tune,
     alpha.p.tune=alpha.p.tune, beta.p.tune=beta.p.tune)
 
}