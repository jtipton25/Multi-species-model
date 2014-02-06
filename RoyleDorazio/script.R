setwd('~/Multi-species-model/RoyleDorazio/')
source('CumNumSpeciesPresent.R')
source('MultiSpeciesSiteOcc.R')

X = as.matrix(read.csv("butterflyData.txt"))
nrepls = 18
butterfly = MultiSpeciesSiteOcc(nrepls, X)

summary(butterfly$fit$sims.matrix)

alpha.post = butterfly$fit$sims.matrix[,"alpha"]
sigmaU.post = butterfly$fit$sims.matrix[,"sigma.u"]
N.post = butterfly$fit$sims.matrix[,"N"]

nsites = 8
CumNumSpeciesPresent(nsites, alpha.post, sigmaU.post, N.post)
$meanAndquantiles

$summaryStats

  