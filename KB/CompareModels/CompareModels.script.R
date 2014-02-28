
#  Read in data/ rcreate simulated data, 
#  then use to compare our model to the Dorazio model

setwd("~/Multi-species-model/KB/CompareModels")

###
### Option 1:  Use Dorazio data
###

bbDat <- read.table("~/Multi-species-model/KB/breedingBirdData.txt", header=T, sep=",")
head(bbDat)
butterflies <- read.table("~/Multi-species-model/KB/butterflyData.txt", header=T, sep=",")
head(butterflies)

y <- as.matrix(butterflies)
n <- dim(y)[1]
J <- dim(y)[2]
K <- 18  # 11 sampling occassions for BBS data, 18 for the butterflies
nzeroes <- 100
X <- rbind(as.matrix(y), matrix(0, nrow=nzeroes, ncol=J) )
Z <- ifelse(X>0, 1, NA)
w <- c( rep(1, n), rep(NA, nzeroes) )


###
### Or simulate data ----------------------
###

FOLDER <- "Sim3_25sites"

K <- 3    # number of surveys
N <- 100  # true number of species
J <- 25 # 10  # number of sites
nzeroes <- 150

ALPHA.p <- 1
BETA.p <-  8
ALPHA.psi <- 6  # 2  #0.4  #3
BETA.psi <- 2  #0.4  #3

#***** Assuming that detection and occupancy probabilities are independent.
det.prob <- rbeta(N, ALPHA.p, BETA.p)
occu.prob <- rbeta(N, ALPHA.psi, BETA.psi)

z <- matrix(NA, nrow=N, ncol=J)
tmp.y <- z
for(k in 1:N){
  z[k, ] <- rbinom(J, 1, occu.prob[k])
  tmp.y[k, ] <- rbinom(J, K, det.prob[k]) * z[k, ]
}
(realizedN <- length( which(apply(z, 1, sum) > 0) ) )  # N = 200 of the 200 possible birds actually occupy a site
y <- tmp.y[which(apply(tmp.y, 1, sum) > 0), ]
dim(y)   # 274 of 300 birds ever observed

n <- dim(y)[1]
X <- rbind(as.matrix(y), matrix(0, nrow=nzeroes, ncol=J) )
Z <- ifelse(X>0, 1, NA)
w <- c( rep(1, n), rep(NA, nzeroes) )



###
###  Run All models
###

out <- list()
Nhat <- matrix(NA, nrow=7, ncol=3)
times <- Nhat
nt <- 2
niter <- 10000
nb <- 1000
nc <- 2

###
### Run Dorazio Model --------
###

source("~/Multi-species-model/KB/DorazioModel/DorazioModel.R")
tmp <- system.time( out$D <- jags(data, inits, params, model.file=DorazioModel, n.chains = nc,
   n.thin = nt, n.iter = niter, n.burnin = nb)  )
print(tmp)
times[1, ] <- tmp[1:3]
print(out$D, 2)
#par(mfrow=c(3, 3))
#traceplot(out$D)
#plot(out$D)

(Nhat[1, ] <- quantile(out$D$BUGSoutput$sims.list$N, c(0.025, 0.5, 0.975) ))

#out.update <- update(out$D, n.iter=10000)
#print(out.update, 2)

jagsfit.mcmc <- as.mcmc.list(out$D$BUGSoutput)
quartz()
xyplot(jagsfit.mcmc, main="Dorazio")
densityplot(jagsfit.mcmc, main="Dorazio")

jpeg(paste(FOLDER, "/DorazioPlots_trace.jpg", sep="") )
xyplot(jagsfit.mcmc, main="Dorazio")
dev.off()

jpeg(paste(FOLDER, "/DorazioPlots_posterior.jpg", sep="") )
densityplot(jagsfit.mcmc, main="Dorazio")
dev.off()

###
### Run our model in JAGS ------
###

source("OurModel_jags.R")
tmp <- system.time( out$JAGS <- jags(data, inits, params, model.file=OurModel, n.chains = nc,
   n.thin = nt, n.iter = niter, n.burnin = nb)  )
print(tmp)
times[2, ] <- tmp[1:3]
print(out$JAGS, 2)
(Nhat[2, ] <- quantile(out$JAGS$BUGSoutput$sims.list$N, c(0.025, 0.5, 0.975) ))

jagsfit.mcmc <- as.mcmc.list(out$JAGS$BUGSoutput)
quartz()
xyplot(jagsfit.mcmc, main="Our JAGS")
densityplot(jagsfit.mcmc, main="Our JAGS")

jpeg(paste(FOLDER, "/OurJAGSPlots_trace.jpg", sep="") )
xyplot(jagsfit.mcmc, main="Our JAGS")
dev.off()

jpeg(paste(FOLDER, "/OurJAGSPlots_posterior.jpg", sep="") )
densityplot(jagsfit.mcmc, main="Our JAGS")
dev.off()
###
    
###
### Run our model in JAGS--fixed effects-----
###

source("OurModel_jags.R")
tmp <- system.time( out$JAGSfixed <- jags(data, inits_fixed, params_fixed, model.file=OurModel_fixed, n.chains = nc,
   n.thin = nt, n.iter = niter, n.burnin = nb)  )
# print(tmp)
times[3, ] <- tmp[1:3]
print(out$JAGSfixed, 2)
#par(mfrow=c(3, 3))
#traceplot(out$D)
#plot(out$JAGS)
(Nhat[3, ] <- quantile(out$JAGSfixed$BUGSoutput$sims.list$N, c(0.025, 0.5, 0.975) ))

jagsfit.mcmc <- as.mcmc.list(out$JAGSfixed$BUGSoutput)
quartz()
xyplot(jagsfit.mcmc, main="Our JAGS, fixed p")
densityplot(jagsfit.mcmc, main="Our JAGS, fixed p")

jpeg(paste(FOLDER, "/OurJAGSfixedpPlots_trace.jpg", sep="") )
xyplot(jagsfit.mcmc, main="Our JAGS, fixed p")
dev.off()

jpeg(paste(FOLDER, "/OurJAGSfixedpPlots_posterior.jpg", sep="") )
densityplot(jagsfit.mcmc, main="Our JAGS, fixed p")
dev.off()
    
###
###  Run my mcmc code-----
###

sink(paste(FOLDER, "/overall_results.txt", sep="") )
cat("True species is ", N, ". \n But \'Realized N\' is ", realizedN, ".  \n")
cat("Observed species is ", dim(y)[1], ". \n")
cat("ALPHA.p=  ", ALPHA.p, ",  BETA.p= ", BETA.p, "\n")
cat("ALPHA.psi=  ", ALPHA.psi, ",  BETA.psi= ", BETA.psi, "\n")
cat("Broms model, all random effects.  \n \n")  
source("~/Multi-species-model/KB/multispp.occu.simple.mcmc_LessIndexing.R")
tune <- 0.1
tmp <- system.time( out$me <- multispp.occu.simple.mcmc(y, K, N=median(out$D$BUGSoutput$sims.list$N), 
            N.MCMC=2*niter,
			alpha.psi.tune=tune, beta.psi.tune=tune,
			alpha.p.tune=tune, beta.p.tune=tune, drawPlots=T) )
# print(tmp)
times[4, ] <- tmp[1:3]
(Nhat[4, ] <- quantile(out$me$N.save[-(1:out$me$n.burn)], c(0.025, 0.5, 0.975) ) )

  
###
###  Run my mcmc code - p fixed effects with uniform prior------
###

quartz()
cat("\n \n \n Broms model, fixed p with uniform prior.  \n \n")  
source("~/Multi-species-model/KB/multispp.occu.simple.mcmc_LessIndexing.R")
tune <- 0.1
tmp <- system.time( out$pFlat <- multispp.occu.fixedp.mcmc(y, K, N=median(out$D$BUGSoutput$sims.list$N), 
            N.MCMC=2*niter,
			alpha.psi.tune=tune, beta.psi.tune=tune,
			drawPlots=T) )
# print(tmp)
times[5, ] <- tmp[1:3]
(Nhat[5, ] <- quantile(out$pFlat$N.save[-(1:out$pFlat$n.burn)], c(0.025, 0.5, 0.975) ) )

###
###  Run my mcmc code - p fixed effects with informative prior-----
###

cat("\n \n \n Broms model, fixed p with informative prior.  \n \n")  
source("~/Multi-species-model/KB/multispp.occu.simple.mcmc_LessIndexing.R")
tune <- 0.1
tmp <- system.time( out$pInform <- multispp.occu.fixedp.mcmc(y, K, N=median(out$D$BUGSoutput$sims.list$N), 
            N.MCMC=2*niter,
			alpha.psi.tune=tune, beta.psi.tune=tune,
			alpha.p=1.5, beta.p=8,
			drawPlots=T) )
# print(tmp)
times[6, ] <- tmp[1:3]
(Nhat[6, ] <- quantile(out$pInform$N.save[-(1:out$pInform$n.burn)], c(0.025, 0.5, 0.975) ) )


rownames(times) <- c("Dorazio", "JAGS", "FixedJAGS", "Broms",
  "Broms_uniformP", "Broms_informedP", "Tipton")
rownames(Nhat) <- c("Dorazio", "JAGS", "FixedJAGS", "Broms",
  "Broms_uniformP", "Broms_informedP", "Tipton")
#Save outputs!!!!
# sink("BBS_Comparison/BBS_results.txt")
# sink("Butterfly_Comparison/butterfly_results.txt")

#sink("Sim1/overall_results.txt")
cat("\n Dorazio Model:  \n")
	print(out$D, 2)
    cat("\n\n\n\n")
cat( "\n JAGS Model:  \n")
print(out$JAGS, 2)# cat("ERROR, didn't run. \n")
    cat("\n\n\n\n")    
cat( "\n JAGS Model, fixed p:  \n")
print(out$JAGSfixed, 2)
    cat("\n\n\n\n")
cat("\n\n\n\n\n  Run times:  \n")
times
cat("\n\n\n\n\n  Spp Richness with CI's:  \n")
Nhat
sink()


###
###  Run Tipton's code ------
###

#rm(n); rm(J); rm(Z); rm(w)
source("~/Multi-species-model/JT/mcmcMsJT.R")
tune <- 0.1
#tmp <- system.time( out$Tipton <- mcmcMS(Y=t(y), J=K, n.aug=100,
#  alpha.alpha.p=2, beta.alpha.p=1, alpha.beta.p=2, beta.beta.p=1, 
#  alpha.alpha.psi=2, alpha.beta.psi=1, beta.alpha.psi=2, beta.beta.psi=1, 
#  alpha.lambda=1, beta.lambda=1, 
#  alpha.p.tune=tune, beta.p.tune=tune, 
#  alpha.psi.tune=tune, beta.psi.tune=tune, 
#  n.mcmc=2*niter, Z.init=0.1) )

#print(tmp)
#times[7, ] <- tmp[1:3]
#(Nhat[7, ] <- quantile(out$Tipton$N.save[-(1:out$Tipton$n.burn)], c(0.025, 0.5, 0.975) ) )

