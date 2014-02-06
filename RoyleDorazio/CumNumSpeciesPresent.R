CumNumSpeciesPresent = function(nsites, alpha, sigmaU, N) {

  # Computes a sample of the posterior-predictive distribution of the (cumulative) number of species present at nsites.

  # compute posterior predictions of species occurrence probabilities
  ndraws = length(alpha)
  Nmax = max(N)
  logitPsi = matrix(NA, nrow=ndraws, ncol=Nmax)
  psi = logitPsi
  for (i in 1:ndraws) {
    logitPsi[i,1:N[i]] = rnorm(N[i], mean=alpha[i], sd=sigmaU[i])
    psi[i, 1:N[i]] = 1/(1+exp(-logitPsi[i,1:N[i]]))
  }
  
  # compute posterior predictions of species presence at each site
  z = array(NA, dim=c(ndraws, Nmax, nsites))
  for (i in 1:ndraws) {
    for (j in 1:N[i]) {
      z[i,j, ] = rbinom(nsites, size=1, prob=psi[i,j])
    }
  }

  # compute posterior predictions of cumulative number of species present
  M = matrix(NA, nrow=ndraws, ncol=nsites)
  for (i in 1:ndraws) {
    for (j in 1:nsites) {
      zsum = rep(NA, N[i])
      if (j>1) {
        zsum = apply(z[i, 1:N[i], 1:j], 1, sum)
      }
      else {
        zsum = z[i, 1:N[i], 1]
      }
      M[i,j] = sum(zsum>0)
    }
  }

  # compute summary stats for plotting
  nSpeciesPresent = matrix(NA, nrow=3, ncol=nsites)
  for (j in 1:nsites) {
    x = M[,j]
    nSpeciesPresent[1, j] = mean(x)
    nSpeciesPresent[2:3, j] = quantile(x, probs=c(.05, .95))
  }

  # plot results
  ylimits = c(min(nSpeciesPresent[2,]), max(nSpeciesPresent[3,]))
  plot(1:nsites, nSpeciesPresent[1,], pch=16, ylim=ylimits, type='b',
       xlab='Number of sample locations', ylab='Number of species', las=1, cex.axis=1.2, cex.lab=1.5, cex=1.5)
  segments(1:nsites, nSpeciesPresent[2,], 1:nsites, nSpeciesPresent[3,])

  list(meanAndquantiles=nSpeciesPresent, summaryStats=summary(M))
}

