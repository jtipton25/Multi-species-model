##
## HW 6 abundance model
##
mcmc.fit <- function(y, x, J, mu.alpha.0, mu.alpha.1, sigma.squared.alpha.0, sigma.squared.alpha.1, alpha.psi, beta.psi, alpha.pi, beta.pi, sigma.squared.alpha.0.tune, sigma.squared.alpha.1.tune, p.x.tune, n.mcmc){

  ##
  ## libraries and subroutines
  ##

  logit <- function(expit){
    log(expit) / (1 - log(expit))
  }

  expit <- function(logit){
    exp(logit) / (1 + exp(logit))
  }

  ##
  ## initialize variables
  ##

  pi.x <- 0.5 ## assume salamander species is equally distributed a priori
  m <- length(y)
  y.0 <- y == 0
  y.1 <- y != 0
  augmented <- length(y[y.0])
  z <- rep(0, length(y))
  z[y.1] <- rep(1, length(y[y.1]))
  x[y.0] <- rbinom(augmented, 1, pi.x)  
  alpha.0 <- rnorm(1, mu.alpha.0, sqrt(sigma.squared.alpha.0))
  alpha.1 <- rnorm(1, mu.alpha.1, sqrt(sigma.squared.alpha.1))
  p <- expit(alpha.0 + alpha.1 * x)
  psi <- rbeta(1, alpha.psi, beta.psi)
  n.burn <- floor(n.mcmc / 10)

  ##
  ## setup save variables
  ##

  z.save <- matrix(0, nrow = length(x), ncol = n.mcmc) - n.burn
  p.save <- matrix(0, nrow = length(x), ncol = n.mcmc - n.burn)
  psi.save <- rep(0, n.mcmc - n.burn)
  alpha.0.save <- rep(0, n.mcmc - n.burn)
  alpha.1.save <- rep(0, n.mcmc - n.burn)
  x.save <- matrix(0, nrow = length(x), ncol = n.mcmc - n.burn)
  pi.save <- rep(0, n.mcmc - n.burn)
  N.save <- rep(0, n.mcmc - n.burn)
  accept.alpha.0 <- 0
  accept.alpha.1 <- 0
  accept.x <- 0


  for(k in 1:n.mcmc){
    if(k %% 100 == 0){
      cat(k, '')
    }    

    ##
    ## Sample Z
    ##

    psi.tilde <- ((1 - p)^J * psi) / ((1 - p)^J * psi + 1 - psi)
    z[y.0] <- rbinom(augmented, 1, psi.tilde)
    z.sum <- sum(z)
    z.1 <- z == 1
    z.0 <- z == 0

    ##
    ## Sample psi
    ##

    psi <- rbeta(1, alpha.psi + z.sum, beta.psi + m - z.sum)

    ##
    ## sample alpha.0
    ## 

    alpha.0.star <- rnorm(1, alpha.0, sigma.squared.alpha.0.tune)
    p.star <- expit(alpha.0.star + alpha.1 * x)
    mh1.alpha.0 <- sum(dbinom(y[z.1], J, prob = p.star[z.1], log = TRUE)) + dnorm(alpha.0.star, mu.alpha.0, sqrt(sigma.squared.alpha.0), log = TRUE)
    mh2.alpha.0 <- sum(dbinom(y[z.1], J, prob = p[z.1], log = TRUE)) + dnorm(alpha.0, mu.alpha.0, sqrt(sigma.squared.alpha.0), log = TRUE)
    mh.alpha.0 <- exp(mh1.alpha.0 - mh2.alpha.0)

    if(mh.alpha.0 > runif(1)){
      alpha.0 <- alpha.0.star
      accept.alpha.0 <- accept.alpha.0 + 1 / n.mcmc
      p <- p.star      
    }    

    ##
    ## sample alpha.1
    ## 

    alpha.1.star <- rnorm(1, alpha.1, sigma.squared.alpha.1.tune)
    p.star <- expit(alpha.0 + alpha.1.star * x)
    mh1.alpha.1 <- sum(dbinom(y[z.1], J, prob = p.star[z.1], log = TRUE)) + dnorm(alpha.1.star, mu.alpha.1, sqrt(sigma.squared.alpha.1), log = TRUE)
    mh2.alpha.1 <- sum(dbinom(y[z.1], J, prob = p[z.1], log = TRUE)) + dnorm(alpha.1, mu.alpha.1, sqrt(sigma.squared.alpha.1), log = TRUE)
    mh.alpha.1 <- exp(mh1.alpha.1 - mh2.alpha.1)

    if(mh.alpha.1 > runif(1)){
      alpha.1 <- alpha.1.star
      accept.alpha.1 <- accept.alpha.1 + 1 / n.mcmc
      p <- p.star    
    }    

    ##
    ## sample x
    ##

    ## For y = 0 and z = 1
    
     x.star <- rbinom(length(x[y.0 & z.1]), 1, p.x.tune)
     p.star <- expit(alpha.0 + alpha.1 * x.star)
     mh.x.1 <- dbinom(y[y.0 & z.1], J, p.star, log = TRUE) + dbinom(x.star, 1, pi.x, log = TRUE)
     mh.x.2 <- dbinom(y[y.0 & z.1], J, p[y.0 & z.1], log = TRUE) + dbinom(x[y.0 & z.1], 1, pi.x, log = TRUE)
     mh.x <- exp(mh.x.1 - mh.x.2)
     idx <- which(mh.x > runif(length(x[y.0 & z.1])))
     x[y.0 & z.1][idx] <- x.star[idx]
     accept.x <- accept.x + (length(x.star[idx]) / length(x.star)) / n.mcmc
#    x.sum <- sum(x[y.0 & z.1])
#    x.total <- length(x[y.0 & z.1])
#    x.sum.star <- x.sum + sample(c(-1,1), 1)
#    if(x.sum.star >= 0 & x.sum.star <= x.total){
#      x.star <- c(rep(1, x.sum.star), rep(0, x.total - x.sum.star))
#      p.star <- expit(alpha.0 + alpha.1 * x.star)
#      mh.x.1 <- sum(dbinom(y[y.0 & z.1], J, p.star, log = TRUE) + dbinom(x.star, 1, pi.x, log = TRUE))
#      mh.x.2 <- sum(dbinom(y[y.0 & z.1], J, p[y.0 & z.1], log = TRUE) + dbinom(x[y.0 & z.1], 1, pi.x, log = TRUE))
#      mh.x <- exp(mh.x.1 - mh.x.2)

#      if(mh.x > runif(1)){
#        x[y.0 & z.1] <- x.star
#        accept.x <- accept.x + 1 / (n.mcmc)
#      }
#    }
  
    ## for y = 0 and z = 0
    x[y.0 & z.0] <- rbinom(length(x[y.0 & z.0]), 1, pi.x)
    
    ## update p
    p <- expit(alpha.0 + alpha.1 * x)

    ##
    ## sample pi.x
    ##

#    pi.x <- rbeta(1, alpha.pi + sum(x[y.0]), beta.pi + augmented  - sum(x[y.0]))
    pi.x <- rbeta(1, alpha.pi + sum(x), beta.pi + m - sum(x))
    
    ##
    ## save variables
    ##
    if(k > n.burn){
      z.save[, k - n.burn] <- z
      p.save[, k - n.burn] <- p
      psi.save[k - n.burn] <- psi
      alpha.0.save[k - n.burn] <- alpha.0
      alpha.1.save[k - n.burn] <- alpha.1
      x.save[, k - n.burn] <- x
      pi.save[k - n.burn] <- pi.x
      N.save[k - n.burn] <- sum(z)
      }
    }

  ##
  ## write output
  ##

  list(z.save = z.save, p.save = p.save, psi.save = psi.save, alpha.0.save = alpha.0.save, alpha.1.save = alpha.1.save, x.save = x.save, pi.save = pi.save, accept.alpha.0 = accept.alpha.0, accept.alpha.1 = accept.alpha.1, accept.x = accept.x, N.save = N.save)
}