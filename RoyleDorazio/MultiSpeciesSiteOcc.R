MultiSpeciesSiteOcc = function(nrepls, X) {

  start.time = Sys.time()
  
  # augment data matrix with an arbitrarily large number of zero row vectors
  nzeroes = 100
  n = dim(X)[1]
  nsites = dim(X)[2]
  Xaug = rbind(X, matrix(0, nrow=nzeroes, ncol=nsites))

  # create arguments for bugs()
  sp.data = list(n=n, nzeroes=nzeroes, J=nsites, K=nrepls, X=Xaug)

  sp.params = list('alpha', 'beta', 'rho', 'sigma.u', 'sigma.v', 'omega', 'N')

  sp.inits = function() {
    omegaGuess = runif(1, n/(n+nzeroes), 1)
    psi.meanGuess = runif(1, .25,1)
    theta.meanGuess = runif(1, .25,1)
    rhoGuess = runif(1, 0,1)
    sigma.uGuess = 1
    sigma.vGuess = 1
    list(omega=omegaGuess, psi.mean=psi.meanGuess, theta.mean=theta.meanGuess, tau.u=1/(sigma.uGuess^2), tau.v=1/(sigma.vGuess^2), rho=rhoGuess,
         w=c(rep(1, n), rbinom(nzeroes, size=1, prob=omegaGuess)),
         phi=rnorm(n+nzeroes, log(psi.meanGuess/(1.-psi.meanGuess)), sigma.uGuess),
         eta=rnorm(n+nzeroes, log(theta.meanGuess/(1.-theta.meanGuess)), sigma.vGuess),
         Z = matrix(rbinom((n+nzeroes)*nsites, size=1, prob=psi.meanGuess), nrow=(n+nzeroes))
         )
  }

  # fit model to data using WinBUGS code
  library(R2WinBUGS)
  fit = bugs(sp.data, sp.inits, sp.params,
    model.file='MultiSpeciesSiteOccModel.txt',
    debug=F, n.chains=4, n.iter=55000, n.burnin=5000, n.thin=50)
  
  end.time = Sys.time()
  elapsed.time = difftime(end.time, start.time, units='mins')
  cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes\n', sep=''))
  
  list(fit=fit, data=sp.data, X=X)
}
