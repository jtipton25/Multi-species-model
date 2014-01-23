MS.occ.mcmc_1yr <- function(Y,k.aug,alpha.psi,sigma.squared.alpha.psi.tune,alpha.alpha.psi,beta.alpha.psi,
                              beta.psi,sigma.squared.beta.psi.tune,alpha.beta.psi,beta.beta.psi,
                              alpha.p,sigma.squared.alpha.p.tune,alpha.alpha.p,beta.alpha.p,
                              beta.p,sigma.squared.beta.p.tune,alpha.beta.p,beta.beta.p,
                              alpha.lambda, beta.lambda, n.mcmc){

  
# Ys <- matrix of a sum of detections for each site, for each species
  
  
#####
#####Subroutines and libraries
#####


library(mvtnorm)
library(Rlab)
library(boot)


#####


##get the number of sites and sum of how many times each site was visited
  
Ys <- apply(Y, 3, rowSums)

n <- dim(Y)[1]
  
J <- dim(Y)[2]

K <- dim(Y)[3]

m <-  dim(Y)[1]


##save variables

#for each mcmc iteration
N.save <- rep(0,n.mcmc)  ##save sum(W) to get species richness
Lambda.save <- rep(0,n.mcmc)

#for each species, for each mcmc iteration

psi.save <- array(0,c(n.mcmc,n,K)) ##saving Psi for each species
p.save <- array(0,c(n.mcmc,n,K)) ##saving p for each species


##Set up starting values for each species, and access them through the K loop

##setting up the zi vector, and defining where we know that it is = 1, and the site was occupied
zi <- matrix(0,nrow=n,ncol=K)

zi[Ys >0] <- 1

Y.nz.idx <- zi==1

z.sum <- colSums(zi)


#####
#####Starting values and priors
#####

##starting values for psi and p
##should we use informative priors here?

# p <- rbeta(K, alpha.p, beta.p)
# psi <- rbeta(K, alpha.psi, beta.psi)

p <- rep(0,K) 
psi <- rep(0,K) 

for(i in 1:K){ 
p[i] <- mean(Ys[Ys[,i]>0,i]/J)
psi[i] <- z.sum[i]/n
}

n.burn <- floor(n.mcmc / 10)

##starting values for w

w.aug <- rbern(k.aug,p=lambda) 

W <- c(rep(1,K),w.aug)

w <- length(W)

w.sum <- sum(W)

#####
#####Star MCMC loop
#####

for(k in 1:n.mcmc){
  if(k %% 100 == 0){
    cat(k, '')
  }    
  
    ##
    ## Sample lambda
    ##
    
    lambda <- rbeta(1,alpha.lambda + w.sum, beta.lambda + (w - w.sum))
    
    ##
    ## Sample W
    ##
    
    w.aug <- rbern(k.aug,p=lambda)
    
    W <- c(rep(1,K),w.aug)
    
    w <- length(W)
    
    w.sum <- sum(W)
  
  
  
      for(i in 1:K){ ##we only want to know psi and p for the species we know are there, but borrowing 
                     ##information about p and psi from the species of other detected species
                     ##this is why estimates are succeptable to the info from species in their group
                     ##but how can you estimate psi and p for the species that were never detected?

          ##
          ## Sample Z when zi[,i]==0
          ##
          
          psi.tilde <- ((1 - p[i])^J * psi[i]) / ((1 - p[i])^J * psi[i] + (1 - psi[i]))
          zi[!Y.nz.idx[,i],i] <- rbern(sum(!Y.nz.idx[,i]),psi.tilde) 
          z.sum.k <- sum(zi[,i])
          
          ##
          ## Sample psi
          ##
          
          psi <- rbeta(1, alpha.psi + z.sum.k, beta.psi + (m - z.sum.k))##is this right?      
          
          ####
          ####  Sample alpha.psi and beta.psi 
          ####
          
          alpha.psi.star <- rnorm(1, alpha.psi, sigma.squared.alpha.psi.tune)
          beta.psi.star <- rnorm(1, beta.psi, sigma.squared.beta.psi.tune)
          
          
          mh1.psi <- sum(dbeta(psi[i],alpha.psi.star,beta.psi.star))+ dgamma(alpha.psi.star,alpha.alpha.psi,beta.alpha.psi)
                     + dgamma(beta.psi.star,alpha.beta.psi,beta.beta.psi)
          mh2.psi <- sum(dbeta(psi[i],alpha.psi,beta.psi))+ dgamma(alpha.psi,alpha.alpha.psi,beta.alpha.psi)
                      + dgamma(beta.psi,alpha.beta.psi,beta.beta.psi)
          mh.psi <- exp(mh1.psi-mh2.psi)
          
          if(mh.psi > runif(1)){
            alpha.psi=alpha.psi.star
            beta.psi=beta.psi.star 
          }
          
          
          ####
          ####  Sample p when Ys > 0, or Y.nz.idx
          ####
          
          p[i] <- rbeta(1, alpha.p + sum(Ys[,i]), beta.p + (sum(Y.nz.idx[,i]*J) - sum(Ys[,i])))
          
          ####
          ####  Sample alpha.p and beta.p
          ####
          
          alpha.p.star <- rnorm(1, alpha.p, sigma.squared.alpha.p.tune)
          beta.p.star <- rnorm(1, beta.p, sigma.squared.beta.p.tune)

          mh1.p <- sum(dbeta(p[i],alpha.p.star,beta.p.star))+ dgamma(alpha.p.star,alpha.alpha.p,beta.alpha.p)
          + dgamma(beta.p.star,alpha.beta.p,beta.beta.p)
          mh2.p <- sum(dbeta(p[i],alpha.p,beta.p))+ dgamma(alpha.p,alpha.alpha.p,beta.alpha.p)
          + dgamma(beta.p,alpha.beta.p,beta.beta.p)
          mh.p <- exp(mh1.p-mh2.p)
          
          if(mh.p > runif(1)){
            alpha.p=alpha.p.star
            beta.p=beta.p.star 
          }
          
          
          #####
          #####Save p, psi for all k iterations
          #####
  
          psi.save[,,i] <- psi
          p.save[,,i] <- p        
        
      }
      
      
      N.save[,k] <- sum(W)
      Lambda.save[,k] <- lambda
      
    }
      

list(psi.save=psi.save,Lambda.save=Lambda.save,N.save = N.save,p.save=p.save)


}