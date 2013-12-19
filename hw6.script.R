setwd('/Volumes/Data Drive 1/Google Drive/Stats/3 - Fall 2013/Bayesian Ecological Modelling/biweek 6/homework 6/')
setwd('~/Google Drive/Stats/3 - Fall 2013/Bayesian Ecological Modelling/biweek 6/homework 6/')

data <- read.table(file = 'hw6.txt', header = TRUE, sep = '')

source('hw6.mcmc.R')

J <- 4
observed <- apply(data[, 1:J], 1, sum)
y <- c(observed, rep(0,1000))
x <- c(data$spp, rep(0, 1000))
mu.alpha.0 <- -2
mu.alpha.1 <- -1
sigma.squared.alpha.0 <- .5
sigma.squared.alpha.1 <- 1
alpha.psi <- .5
beta.psi <- 1
alpha.pi <- 1.5
beta.pi <- 1
n.mcmc <- 10000

sigma.squared.alpha.0.tune <- 0.15
sigma.squared.alpha.1.tune <- 0.20
p.x.tune <- 0.5

mcmc.out <- mcmc.fit(y, x, J, mu.alpha.0, mu.alpha.1, sigma.squared.alpha.0, sigma.squared.alpha.1, alpha.psi, beta.psi, alpha.pi, beta.pi, sigma.squared.alpha.0.tune, sigma.squared.alpha.1.tune, p.x.tune, n.mcmc)
mcmc.out$accept.x

layout(matrix(1:4, nrow = 2))
  plot(mcmc.out$alpha.0.save, type = 'l', main =paste('alpha.0 acceptance rate', round(mcmc.out$accept.alpha.0, digits = 3)))
  plot(mcmc.out$alpha.1.save, type = 'l', main =paste('alpha.1 acceptance rate', round(mcmc.out$accept.alpha.1, digits = 3)))
  plot(mcmc.out$psi.save, type = 'l', main =paste('plot of psi'))
  hist(mcmc.out$N.save, freq = FALSE)
  

  plot(mcmc.out$N.save, type = 'l')

apply(mcmc.out$x.save, 1, mean)
