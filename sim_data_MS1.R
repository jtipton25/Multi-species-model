sim_data_ms1 <- function(K,psi,p,n,J){ 

#   K=25
#   psi <- runif(K,0.4,0.8)
#   p <- runif(K,0.4,0.6)
#   n=250
#   J=5
  
  library(Rlab)
  library(mvtnorm)
  library(boot)
  
  ###setting up the data for each species
  Y.all <- NULL;
  for(i in 1:K){ 
    zi.temp <- rbern(n,round(psi[i],2))
    p.temp <- round(p[i],2)
    Y.temp <- matrix(NA,ncol=J,nrow=n)#J days and n sites
    head(Y.temp)
      for(j in 1:n){
      Y.temp[j,]  <- rbinom(J, size = 1, prob = p.temp*(zi.temp[j]))
      }
    Y.all <- array( c(Y.temp,Y.all),dim = c( n , J , i ))
  }
dim(Y.all)
  
  list(Y.all=Y.all)

}