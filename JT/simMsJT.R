#
# Libraries and functions
##

##
## need to figure out how to use covariates
##

makeMultiSpec <- function(n, N, J, alpha.p, beta.p, alpha.psi, beta.psi){
  ## try simulating from the hyper prior distribution
  
  Y.full <- matrix(nrow = n, ncol = N)
  Z <- matrix(nrow = n, ncol = N)
  p <- rbeta(N, alpha.p, beta.p)
  psi <- rbeta(N, alpha.psi, beta.psi)
  
  
  ##
  ## Sample Z[i, k]
  ##
  
  for(i in 1:n){
    for(k in 1:N){
      Z[i, k] <- rbinom(1, 1, psi[k])
    }
  }
  
  ##
  ## Sample y[i, k] conditional on Z[i, k]
  ##
  
  for(i in 1:n){
    for(k in 1:N){
      Y.full[i, k] <- if(Z[i, k] == 0){
        0
      } else {
        sum(rbinom(J, 1, p[k]))
      }
    }
  }
  
  Y <- Y.full[, which(apply(Y.full, 2, sum) != 0)]
  
  ##
  ## return data
  ##
  
  list(Y = Y, Z = Z, Y.full = Y.full, psi = psi, p = p)
}

