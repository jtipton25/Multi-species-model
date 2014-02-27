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

