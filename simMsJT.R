<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
set.seed(1)
>>>>>>> d7e1c3221d5f628d3888fcd74648c5f4eb7daff4
>>>>>>> ed83440c2d77ad34d15e73f302852531a43dd8f4

##
## Libraries and functions
##

<<<<<<< HEAD
=======
<<<<<<< HEAD
=======

##
## Initialize simulation parameters
##

n <- 10
N <- 500
J <- 2
Omega <- 2000

# presence probability
alpha.psi <- 1 
beta.psi <- 3
psi <- rbeta(N, alpha.psi, beta.psi) # alpha.psi = 1, beta.psi = 3
			
# detection probability
alpha.p <- 2 
beta.p <- 4
p <- rbeta(N,alpha.p, beta.p) # alpha.p = 1, beta.p = 3
lambda <- 0.25


>>>>>>> d7e1c3221d5f628d3888fcd74648c5f4eb7daff4
>>>>>>> ed83440c2d77ad34d15e73f302852531a43dd8f4
makeMultiSpec <- function(n, N, J, alpha.psi, beta.psi, alpha.p, beta.p){
	psi <- rbeta(N, alpha.psi, beta.psi) # alpha.psi = 1, beta.psi = 3
	p <- rbeta(N,alpha.p, beta.p) # alpha.p = 1, beta.p = 3
	Y.full <- matrix(nrow = n, ncol = N)
	Z <- matrix(nrow = n, ncol = N)
	
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
	
	list(Y = Y, Z = Z, Y.full = Y.full)
}

