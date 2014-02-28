
#  Code to learn and practice rstan

library(rstan)
#set_cppo("fast")
 set_cppo("debug")

#setwd("~/Dropbox/STAN")

multispp_code <- '
  data {
  	int<lower=1> J;  // number of sites
    int<lower=1> K;  // number of sampling occassions per site
  	int<lower=0> nzeroes;  // n + nzeroes
  	int<lower=1> n;  // number of observed species
  	int<lower=n> Omega;  // number of observed species
  	int<lower=0, upper=K> Y[Omega, J];  // matrix of the detections for each species and site
  	int<lower=0, upper=K> Y_ind[Omega, J];  // Y converted to 0s and 1s
  
    int<lower=0, upper=1> w_obs[n];  // Need to break up W and Z matrices into observed/unobserved
    int<lower=0, upper=1> Z_obs[n];  // Need to break up W and Z matrices into observed/unobserved
  	int<lower=0> length_z1;
  	int<lower=0> length_z0;

    int<lower=0> Site_ct0[n];
    int<lower=0> Site_ct1[n];
  	int<lower=0, upper=1> Z_y1[length_z1];

  }
  
  parameters {
  	real<lower=0> psi_alpha;
  	real<lower=0> psi_beta;
  	real<lower=0> p_alpha;
  	real<lower=0> p_beta;
  	real<lower=0, upper=1> lambda;
  	
   	real psi[Omega];
  	real det_prob[J];
  }
  model {
    int w_miss[nzeroes];
    int Z_miss[nzeroes, J];
    int Z_y0[length_z0];

  	lambda ~ beta(1, 1);
  	
  	psi_alpha ~ gamma(2, 1);
  	psi_beta ~ gamma(2, 1);
  	p_alpha ~ gamma(2, 1);
  	p_beta ~ gamma(2, 1);

  	psi ~ beta(psi_alpha, psi_beta);
 	det_prob ~ beta(p_alpha, p_beta);
	
	w_obs ~ bernoulli(lambda);
	w_miss ~ bernoulli(lambda);
	
	for (i in 1:nzeroes) {
		if (w_miss[i] == 1) {
			Z_miss[i, ] ~ bernoulli(psi[i]);
		} else{
			Z_miss[i, ] <- 0;
		}
	}
  	for(i in 1:n) { 
  		ct1_tmp <- 0
  		ct0_tmp <- 0
  		for (k in 1:(n-1)) {
  			Z_y1[ (1+Site_ct1[k]):(ct1_tmp + Site_ct1[k+1] ] ~ bernoulli(psi[k]);
  			ct1_tmp <- ct1_tmp + Site_ct1[k+1];

  			Z_y0[ (1+Site_ct0[k]):(ct0_tmp + Site_ct0[k+1] ] ~ bernoulli(psi[k]);
  			ct0_tmp <- ct0_tmp + Site_ct0[k+1];

  		}	
  		Z_obs[i] ~ binomial(K, psi[i]);	
  		for(j in 1:J) {
  			Y[i, j] ~ binomial(K, det_prob[i]*Z[i,j]);
  		}
  	}
  }
'

fit <- stan(model_code = multispp_code, data=multi_dat, iter = 10, chains=2)



Y_ind <- ifelse(y>0, 1, 0)
Site_ct1 <- apply(Y_ind, 1, sum)
Site_ct0 <- J - Site_ct1
z_y1 <- rep(1, sum(Site_ct1) )
# z_y0  ## a transformed parameter (transformed because we don't care about it.)

Site_ct1 <- c(0, Site_ct1)
Site_ct0 <- c(0, Site_ct0)
length_z1 <- sum(Site_ct1)
length_z0 <- sum(Site_ct0)
nzeroes <- 100
n <- dim(y)[1]
multi_dat <- list(J = dim(y)[2],  K = K, nzeroes=nzeroes, n=dim(y)[1], Omega=n+nzeroes,
                  y = rbind(y, matrix(0, nrow=nzeroes, ncol=J)),
                  w_obs = rep(1, n) )
                    
fit <- stan(model_code = multispp_code, data=multi_dat, iter = 10, chains=2)

fit1 <- stan(file = '8schools.stan', data=schools_dat, iter=5000, chains=2)
fit2 <- stan(fit = fit1, data=schools_dat, iter=5000, chains=2)
# fit2 is an S4 object of class "stanfit"
print(fit2)
plot(fit2)

la <- extract(fit2, permuted = TRUE)  # each parameter is in a separate list.
mu <- la$mu 
mean(mu)
str(la)

a <- extract(fit2, permuted = FALSE)  # an n.iter X n.chain X n.parameter array of the output
str(a)

a2 <- as.array(fit2)  # same result as a
str(a2)
?as.array.stanfit

m <- as.matrix(fit2) # an (n.iter*n.chain) X n.parameter matrix
str(m)



###
### alternate way to input data
###

schools_dat <- read_rdump("8schools.rdump")
fit <- stan(model_code = schools_code, data=schools_dat, iter = 1000, chains=2)

#  instead of above two lines
source("8schools.rdump")  # 
fit <- stan(file="8schools.stan", data=c("J", "y", "sigma"), iter=1000, chains=2)


#####################################################################
######################################################################

###
### RATS
###

y <- read.table("rats.txt", header=TRUE)
head(y)
x <- c(8, 15, 22, 29, 36)
rats_dat <- list(N = nrow(y), T = ncol(y),
                 x = x, y = y, xbar = mean(x))
rats_fit <- stan(file = "rats.stan", data=rats_dat, verbose=F)


###
### Running parallel chains
###

example(sflist2stanfit)
scode <- "
  parameters {
  	real y[2];
  }
  model {
  	y[1] ~ normal(0, 1);
  	y[2] ~ double_exponential(0, 2);
  }
"
seed <- 123
f1 <- stan(model_code = scode, chains=1, seed=seed, chain_id=1)
f2 <- stan(fit=f1, chains=2, seed=seed, chain_id=2:3)
f12 <- sflist2stanfit( list(f1, f2) )
print(f12)

library(parallel)
# Run 4 chains in parallel after the model has been compiled
sflist1 <- mclapply(1:4, mc.cores=4, 
                    function(i)  stan(fit=f1, seed = seed, chains=1, chain_id=i, refresh=-1) )
f3 <- sflist2stanfit(sflist1)
print(f3)