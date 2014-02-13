multispp.occu.simple.mcmc <- function(y, K=3, nzeroes=100, N=N,
									  alpha.lambda=3, beta.lambda=5, 
                                      prior.shape.alpha.p=6, prior.rate.alpha.p=3,
                                      prior.shape.beta.p=16, prior.rate.beta.p=4,
                                      prior.shape.alpha.psi=12, prior.rate.alpha.psi=4,
                                      prior.shape.beta.psi=4, prior.rate.beta.psi=4,
                                      alpha.psi.tune=0.1, beta.psi.tune=0.1,
                                      alpha.p.tune=0.1, beta.p.tune=0.1,
                                      drawPlots=T, N.MCMC=1000){
  
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
  
  J <- dim(y)[2]  # previously called n
  n <- dim(y)[1]  # previously K
  Omega <- n + nzeroes
  psi <- rep(NA, Omega)
  p <- rep(NA, Omega)
  mat0 <- matrix(rep(0, J*nzeroes), ncol=J)
  Y <- rbind(y, mat0)  # augmented data matrix
  psi.acc <- 0
  p.acc <- 0
  
  psi.alpha.save <- rep(0, N.MCMC)
  psi.beta.save <- rep(0, N.MCMC)
  p.alpha.save <- rep(0, N.MCMC)
  p.beta.save <- rep(0, N.MCMC)
  N.save <- rep(0, N.MCMC)

  n.burn <- N.MCMC/10
  
  #####
  #####  Priors and Starting Values 
  #####

  y.k <- apply(Y, 1, sum)  # was a species detected at any site?
  z <- ifelse(Y > 0, 1, NA)
  w <- c(rep(1, n), rep(0, nzeroes) )
   
   # set priors equal to true values.
  alpha.psi <- rgamma(1, prior.shape.alpha.psi, prior.rate.alpha.psi) #ALPHA.psi  # 
  beta.psi <- rgamma(1, prior.shape.beta.psi, prior.rate.beta.psi) # BETA.psi  # 
  alpha.p <- rgamma(1, prior.shape.alpha.p, prior.rate.alpha.p)  # ALPHA.p  # 
  beta.p <- rgamma(1, prior.shape.beta.p, prior.rate.beta.p)  # BETA.p  # 
  
  psi <- rbeta(Omega, alpha.psi, beta.psi)
  p <- rbeta(Omega, alpha.p, beta.p)
  lambda <- rbeta(1, alpha.lambda, beta.lambda)

 
  
  #####
  #####  Begin MCMC Loop 
  #####  
  for(n.mcmc in 1:N.MCMC){
    if(n.mcmc %% 100 == 0) cat(n.mcmc, " "); flush.console()


 	#####
	##### Sample w_k's
	#####
	
	z.k <- apply(z, 1, sum, na.rm=T)
	w[z.k==1] <- 1
	
	lambda.tilde.num <- ( lambda * (1-psi) ^ J )
	lambda.tilde <- lambda.tilde.num / (lambda.tilde.num + 1 - lambda)
	w[z.k==0]  <- rbinom(sum(z.k==0), 1, lambda.tilde[z.k==0])


    #####
    #####  Sample z 
    #####
 
    psi.tilde <- ( psi * (1-p) ^ K ) / ( psi * (1-p) ^ K + 1-psi )     
    for (k in 1:Omega) {  # Loops are sooooooo sloooooooooooooooow!
      if (w[k] == 1) {
      	for (i in 1:J) {
      		if (Y[k, i] == 0)   z[k, i] <- rbinom( 1, 1, psi.tilde[k] )     				
      	}
      }
      if (w[k] == 0) {
      	z[k, ] <- 0  # ( 1 - p[k] ) ^ J
      }
    }
    
 	  
    #####
    ##### Sample p
    #####
    
    p <- rbeta(Omega, alpha.p, beta.p)  #rbeta(Omega, apply(Y*z, 1, sum, na.rm=T) + alpha.p, apply(z*(K-Y), 1, sum, na.rm=T) + beta.p)

 
    #####
    ##### Sample psi
    #####
    
    psi[w==1] <- rbeta(sum(w==1), apply(z, 1, sum, na.rm=T)[(w==1)] + alpha.psi, apply(1-z, 1, sum, na.rm=T)[(w==1)] + beta.psi)
    psi[w==0] <- rbeta(sum(w==0), alpha.psi, beta.psi)    
  #    psi <- rbeta(Omega, alpha.psi, beta.psi)    
 
    
    #####
    #####  Sample p parameters
    #####

    alpha.p.star <- rnorm(1, alpha.p, alpha.p.tune)
    beta.p.star <- rnorm(1, beta.p, beta.p.tune)
    
    if(alpha.p.star > 0 && beta.p.star > 0){
      mh1.p <- sum( dbeta(p, alpha.p.star, beta.p.star, log=TRUE) ) + 
        log(dgamma(alpha.p.star, prior.shape.alpha.p, prior.rate.alpha.p) ) + 
        log(dgamma(beta.p.star, prior.shape.beta.p, prior.rate.beta.p) )
      mh2.p <- sum( dbeta(p, alpha.p, beta.p, log=TRUE) ) + 
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
      mh1.psi <- sum( dbeta(psi, alpha.psi.star, beta.psi.star, log=TRUE) ) + 
        dgamma(alpha.psi.star, prior.shape.alpha.psi, prior.rate.alpha.psi, log=T) + 
        dgamma(beta.psi.star, prior.shape.beta.psi, prior.rate.beta.psi, log=T)
      mh2.psi <- sum( dbeta(psi, alpha.psi, beta.psi, log=TRUE) ) + 
        dgamma(alpha.psi, prior.shape.alpha.psi, prior.rate.alpha.psi, log=T) + 
        dgamma(beta.psi, prior.shape.beta.psi, prior.rate.beta.psi, log=T )
   
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
#    cat("\t", sum(w), "\t")

  }
  cat("\n")



  if(drawPlots==T){
	####
	####  Make Plots
	####
    
	layout( matrix(c(1:5, 6, 6, 7, 7, 8), 5, 2 ) )
	par(mar=c(4, 3, 3, 1)+0.1)
	plot( p.alpha.save[-(1:n.burn)], type="l", lty=1, main="Detection prob. Alpha", ylab="")
	if (exists("ALPHA.p")==T)   abline(h= ALPHA.p)
	plot( p.beta.save[-(1:n.burn)], type="l", lty=1 , main="Detection prob. Beta", ylab="")
	if (exists("BETA.p")==T)   abline(h= BETA.p)
	plot( psi.alpha.save[-(1:n.burn)], type="l", lty=1, main="Occu Prob alpha", ylab="")
	if (exists("ALPHA.psi")==T)  abline(h= ALPHA.psi)
	plot( psi.beta.save[-(1:n.burn)], type="l", lty=1, main="Occu Prob beta", ylab="")
    if (exists("BETA.psi")==T)  abline(h= BETA.psi)
   plot( N.save, type="l", lty=1, main="Species richness", ylab="")
    abline(h=N, col="red", lty=2)
	
  	curve(dbeta(x, mean(p.alpha.save[-(1:n.burn)]), mean(p.beta.save[-(1:n.burn)]) ), ylab="" )
  	if (exists("ALPHA.p")==T)   curve(dbeta(x, ALPHA.p, BETA.p), add=T, col="red", lty=2)
  	curve(dbeta(x, mean(psi.alpha.save[-(1:n.burn)]), mean(psi.beta.save[-(1:n.burn)]) ), ylab="")
  	if (exists("ALPHA.psi")==T)   curve(dbeta(x, ALPHA.psi, BETA.psi), add=T, col="red", lty=2)
  	
  	hist(N.save)
  	abline(v=N, col="red", lty=2)
  	
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
    
    cat("\n Species richness:  ", quantile(N.save[-(1:n.burn)], c(0.025, 0.5, 0.975) ) )
    cat("\n Mean:  ", mean(N.save[-(1:n.burn)]) )



  #####
  #####  Write Output 
  #####
  
  cat("\n Psi acceptance rate: \t", psi.acc/N.MCMC)
  cat("\n Det. prob acceptance rate: \t", p.acc/N.MCMC, "\n")

  list(psi.alpha.save=psi.alpha.save, psi.beta.save=psi.beta.save,
     p.alpha.save=p.alpha.save, p.beta.save=p.beta.save,
     N.save=N.save, N.MCMC=N.MCMC, n.burn=n.burn, psi.acc=psi.acc, 
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
  y.k <- apply(Y, 2, sum)  # was a species detected at any site?

  alpha.p.tune <- out$alpha.p.tune
  beta.p.tune <- out$beta.p.tune
  alpha.psi.tune <- out$alpha.psi.tune
  beta.psi.tune <- out$beta.psi.tune
  prior.shape.alpha.p <- out$prior.shape.alpha.p 
  prior.rate.alpha.p  <- out$prior.rate.alpha.p
  prior.shape.beta.p  <- out$prior.shape.beta.p
  prior.rate.beta.p <- out$prior.rate.beta.p
  prior.shape.alpha.psi <- out$prior.shape.alpha.psi
  prior.rate.alpha.psi <- out$prior.rate.alpha.psi
  prior.shape.beta.psi <- our$prior.shape.beta.psi
  prior.rate.beta.psi <- out$prior.rate.beta.psi
  
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
	##### Sample w_k's
	#####
	
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
      	}
      }
      if(w[k]==0){
      	z[, k] <- ( 1 - p[k] ) ^ J
      }
    }
    
 	  
    #####
    ##### Sample p
    #####
    
    p<- rbeta(Omega, apply(Y*z, 2, sum, na.rm=T) + alpha.p, apply(z*(J-Y), 2, sum, na.rm=T) + beta.p)

 
    #####
    ##### Sample psi
    #####
    
    psi[w==1] <- rbeta(sum(w==1), apply(z, 2, sum, na.rm=T)[(w==1)] + alpha.psi, apply(1-z, 2, sum, na.rm=T)[(w==1)] + beta.psi)
    psi[w==0] <- rbeta(sum(w==0), alpha.psi, beta.psi)    
    
    
    #####
    #####  Sample p parameters
    #####

    alpha.p.star <- rnorm(1, alpha.p, alpha.p.tune)
    beta.p.star <- rnorm(1, beta.p, beta.p.tune)
    
    if(alpha.p.star > 0 && beta.p.star > 0){
      mh1.p <- sum( dbeta(p, alpha.p.star, beta.p.star, log=TRUE) ) + 
        log(dgamma(alpha.p.star, prior.shape.alpha.p, prior.rate.alpha.p) ) + 
        log(dgamma(beta.p.star, prior.shape.beta.p, prior.rate.beta.p) )
      mh2.p <- sum( dbeta(p, alpha.p, beta.p, log=TRUE) ) + 
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
      mh1.psi <- sum( dbeta(psi, alpha.psi.star, beta.psi.star, log=TRUE) ) + 
        log(dgamma(alpha.psi.star, prior.shape.alpha.psi, prior.rate.alpha.psi) ) + 
        log(dgamma(beta.psi.star, prior.shape.beta.psi, prior.rate.beta.psi) )
      mh2.psi <- sum( dbeta(psi, alpha.psi, beta.psi, log=TRUE) ) + 
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
    
	layout( matrix(c(1:5, 6, 6, 7, 7, 8), 5, 2 ) )
	plot( p.alpha.save, type="l", lty=1, main="Detection prob. Alpha")
	abline(h=alpha.p)
	plot( p.beta.save, type="l", lty=1 , main="Detection prob. Beta" )
	abline(h=beta.p)
	plot( psi.alpha.save, type="l", lty=1, main="Occu Prob alpha")
	abline(h=alpha.psi)
	plot( psi.beta.save, type="l", lty=1, main="Occu Prob beta")
    abline(h=beta.psi)
    plot( N.save, type="l", lty=1, main="Species richness")
    abline(h=dim(y))
	
  	curve(dbeta(x, mean(p.alpha.save), mean(p.beta.save) ) )
  	curve(dbeta(x, alpha.p, beta.p), add=T, col="red")
  	curve(dbeta(x, mean(psi.alpha.save), mean(psi.beta.save) ) )
  	curve(dbeta(x, alpha.psi, beta.psi), add=T, col="red")
  	
  	hist(N.save)
  	abline(v=dim(y))

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
  
  cat("\n \n Psi acceptance rate: \t", psi.acc/n.mcmc)
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