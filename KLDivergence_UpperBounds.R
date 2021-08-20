
################################################################################
################################################################################
#
# PURPOSE OF THIS SCRIPT
# This script calculates one-sided upper tolerance bounds for Kullback-Leibler
# divergences between a multivariate normal distribution, representing the
# true phenotypic distribution of a species, and estimates of that distribution
# based on random samples of various sizes. Each one-sided upper tolerance bound
# contains a given percentile of Kullback-Leibler divergences, with a given
# confidence level. The one-side upper bounds are "distribution free," meaning
# no assumptions about the distribution of the Kullback-Leibler divergences are
# used to calculate the one-side upper bounds (Meeker et al. 2017. Statistical
# intervals: A guide for practitioners and researchers. John Wiley & Sons.). In
# this script, one-side upper bounds are calculated for several sample sizes and
# number of dimensions of the phenotypic distribution of species, following the
# computational method in page 80 of Meeker et al. (2017). The particular sample
# sizes, number of dimensions and parameters defining the multivariate normal
# distributions match those in computer simulation experiments designed to
# examine the effect of sample size (number of specimens) on the performance of
# normal mixture models in species delimitation based on phenotypic data.     
#
# FILES REQUIRED TO RUN THE SCRIPT
# None.
#
# CONTENTS OF THIS SCRIPT
# 1. Preliminaries.
# 2. Definition of sample sizes.
# 3. Calculation of upper tolerance bounds for Kullback-Leibler divergence.
# 4. Compare upper bounds across sample sizes and dimensionalities.
# 5. Write file with upper bounds for all sample sizes and dimensionalities.
#
################################################################################
################################################################################


################################################################################
# 1. Preliminaries.
################################################################################

#Load packages
library(mvtnorm)
library(clusterGeneration)

#Define directory to save file with upper bounds for all sample sizes and 
#dimensionalities
dir.to.save.file <- "replace with path to directory to save file"

#Function for calculating Kullback-Leibler divergence between two normal 
#distributions. Obtained from package "rags2ridges"; original code can be
#found at https://rdrr.io/cran/rags2ridges/src/R/rags2ridges.R#sym-KLdiv
KLdiv <- function(Mest, Sest, Mtrue, Strue)
	{
	KLd <- (sum(diag(solve(Sest)%*%Strue)) 
	+ t(Mest - Mtrue)%*%solve(Sest)%*%(Mest - Mtrue)
	- nrow(Strue) - log(det(Strue)) + log(det(Sest)))/2

	return(as.numeric(KLd))
	}


################################################################################
# 2. Definition of sample sizes.
################################################################################

#Define the set of total sample sizes, which are the total number of specimens
#obtained from a species group
n <- seq(50, 1450, 100)

#Declare the vector that will contain the sample sizes per species, which are
#the number of specimens obtained from a single species in a species group 
specimens.per.species <- vector(mode="list", length=length(n))
names(specimens.per.species) <- n

#Run loop to obtain all unique sample sizes per species
for(i in 1:length(n)){

	specimens <- vector(mode="list", length=9)

	#sampling.proportion = 0.5 and s = 8
	specimens[[1]] <- n[i] * c(0.14, 0.14, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12)
	#sampling.proportion = 0.5 and s = 4 
	specimens[[2]] <- n[i] * c(0.26, 0.26, 0.24, 0.24)
	#sampling.proportion = 0.5 and s = 2 
	specimens[[3]] <- n[i] * c(0.5, 0.5)

	#sampling.proportion = 0.7 and s = 8
	specimens[[4]] <- n[i] * c(0.36, 0.16, 0.10, 0.08, 0.08, 0.08, 0.08, 0.06)
	#sampling.proportion = 0.7 and s = 4 
	specimens[[5]] <- n[i] * c(0.50, 0.20, 0.16, 0.14)
	#sampling.proportion = 0.7 and s = 2 
	specimens[[6]] <- n[i] * c(0.7, 0.3)

	#sampling.proportion = 0.9 and s = 8
	specimens[[7]] <- n[i] * c(0.54, 0.20, 0.10, 0.06, 0.04, 0.02, 0.02, 0.02)
	#sampling.proportion = 0.9 and s = 4 
	specimens[[8]] <- n[i] * c(0.72, 0.18, 0.06, 0.04)
	#sampling.proportion = 0.9 and s = 2 
	specimens[[9]] <- n[i] * c(0.9, 0.1)

	specimens.per.species[[i]] <- sort(unique(round(unlist(specimens))))
}

#Obtain and sort unique values of per species sample size 
specimens.per.species <- sort(unique(round(unlist(specimens.per.species))))

#Add the value 16, which is the smallest per species sample size at which the
#Kullback-Leibler divergence can be calculated when dimensionality is p = 15;
#the upper boundary of Kullback-Leibler divergence at per species sample size
#16 will serve as the standard for per species sample sizes < 16. 
specimens.per.species <- sort(c(16, specimens.per.species))


################################################################################
# 3. Calculation of upper tolerance bounds for 
# legend("topright", c("p = 5", "p = 10", "p = 15"), pch=21, lty=1, 
#        col=c("blue", "green", "red"))
################################################################################

################################################################################
# 3.1. Multivariate normal distributions with p = 5 dimensions.

#Define number of dimensions
p <- 5 
#Define number of species to be simulated for each sample size
species <- 15000
#Create matrix that will contain Kullback-Leibler divergence for different
#per species sample sizes
KLdivergence.p5 <- matrix(NA, nrow=length(specimens.per.species), ncol=species)

start.time <- Sys.time() #Start timing loop

for(i in 1:length(specimens.per.species)){
	for(j in 1:species){
		#Simulate the species mean vector
		M.A <- runif(p, -1, 1)
		#Simulate the species variance-covariance, note that variances are sampled
		#from a uniform distribution between 1e-12 and 1, this is the full
		#range of variance values used in the computer simulation experiment
		#designed to examine the effect of sample size on the performance of
		#normal mixture models in species delimitation based on phenotypic data.
		VCV.A <- genPositiveDefMat(covMethod="unifcorrmat", rangeVar=c(1e-12,1), dim=p)$Sigma
		#Sample species
		A.s <- rmvnorm(specimens.per.species[i], mean = M.A, sigma = VCV.A)
		em.A <- colMeans(A.s)
		evcv.A <- var(A.s)	
		KLdivergence.p5[i,j] <- tryCatch(KLdiv(Mest=em.A, Sest=evcv.A, Mtrue=M.A, Strue=VCV.A), error=function(e) NA)
	}
}

Sys.time() - start.time #Finish timing loop

#Examine results
KLdivergence.p5[,1:5]
plot(as.vector(KLdivergence.p5) ~ rep(specimens.per.species, times=species))
plot(log(as.vector(KLdivergence.p5)) ~ rep(specimens.per.species, times=species))

#Calculate upper tolerance bounds for each per species sample size,
#using the computational method in page 80 of Meeker et al. (2017) 
Xperc <- 0.95 #the percentile for which the upper tolerance bounds are calculated
conf <- 0.999 #the confidence level of the upper tolerance bounds
u.index <- qbinom(conf, species, Xperc) + 1
u.l <- function(x) sort(x)[u.index]
KL.p5.ul <- apply(KLdivergence.p5, 1, u.l)

#Examine upper tolerance bounds for different per species sample sizes
plot(specimens.per.species, KL.p5.ul, type="o")
plot(specimens.per.species, log(KL.p5.ul), type="o")
plot(specimens.per.species, log(KL.p5.ul), type="o", xlim=c(0,20), xaxt="n")
axis(1, specimens.per.species)
abline(v=p+1, lty=3)
abline(h=log(0.1), lty=3)

################################################################################
# 3.2. Multivariate normal distributions with p = 10 dimensions.

#Define number of dimensions
p <- 10 
#Define number of species to be simulated for each per species sample size
species <- 15000
#Create matrix that will contain Kullback-Leibler divergence for different
#per species sample sizes
KLdivergence.p10 <- matrix(NA, nrow=length(specimens.per.species), ncol=species)

start.time <- Sys.time() #Start timing loop

for(i in 1:length(specimens.per.species)){
	for(j in 1:species){
		#Simulate the species mean vector
		M.A <- runif(p, -1, 1)
		#Simulate the species variance-covariance, note that variances are sampled
		#from a uniform distribution between 1e-12 and 1 (see comment for p = 5, above)
		VCV.A <- genPositiveDefMat(covMethod="unifcorrmat", rangeVar=c(1e-12,1), dim=p)$Sigma 
		A.s <- rmvnorm(specimens.per.species[i], mean = M.A, sigma = VCV.A)
		em.A <- colMeans(A.s)
		evcv.A <- var(A.s)	
		KLdivergence.p10[i,j] <- tryCatch(KLdiv(Mest=em.A, Sest=evcv.A, Mtrue=M.A, Strue=VCV.A), error=function(e) NA)
	}
}

Sys.time() - start.time #Finish timing loop

#Examine results
KLdivergence.p10[,1:5]
plot(as.vector(KLdivergence.p10) ~ rep(specimens.per.species, times=species))
plot(log(as.vector(KLdivergence.p10)) ~ rep(specimens.per.species, times=species))

#Calculate upper tolerance bounds for each per species sample size,
#using the computational method in page 80 of Meeker et al. (2017) 
Xperc <- 0.95 #the percentile for which the upper tolerance bounds are calculated
conf <- 0.999 #the confidence level of the upper tolerance bounds
u.index <- qbinom(conf, species, Xperc) + 1
u.l <- function(x) sort(x)[u.index]
KL.p10.ul <- apply(KLdivergence.p10, 1, u.l)

#Examine upper tolerance bounds for different per species sample sizes
plot(specimens.per.species, KL.p10.ul, type="o")
plot(specimens.per.species, log(KL.p10.ul), type="o")
plot(specimens.per.species, log(KL.p10.ul), type="o", xlim=c(0,20), xaxt="n")
axis(1, specimens.per.species)
abline(v=p+1, lty=3)
abline(h=log(0.1), lty=3)

################################################################################
# 3.3. Multivariate normal distributions with p = 15 dimensions.

#Define number of dimensions
p <- 15 
#Define number of species to be simulated for each per species sample size
species <- 15000
#Create matrix that will contain Kullback-Leibler divergence for different
#per species sample sizes
KLdivergence.p15 <- matrix(NA, nrow=length(specimens.per.species), ncol=species)

start.time <- Sys.time() #Start timing loop

for(i in 1:length(specimens.per.species)){
	for(j in 1:species){
		#Simulate the species mean vector
		M.A <- runif(p, -1, 1)
		#Simulate the species variance-covariance, note that variances are sampled
		#from a uniform distribution between 1e-12 and 1 (see comment for p = 5, above)
		VCV.A <- genPositiveDefMat(covMethod="unifcorrmat", rangeVar=c(1e-12,1), dim=p)$Sigma 
		A.s <- rmvnorm(specimens.per.species[i], mean = M.A, sigma = VCV.A)
		em.A <- colMeans(A.s)
		evcv.A <- var(A.s)	
		KLdivergence.p15[i,j] <- tryCatch(KLdiv(Mest=em.A, Sest=evcv.A, Mtrue=M.A, Strue=VCV.A), error=function(e) NA)
	}
}

Sys.time() - start.time #Finish timing loop

#Examine results
KLdivergence.p15[,1:5]
plot(as.vector(KLdivergence.p15) ~ rep(specimens.per.species, times=species))
plot(log(as.vector(KLdivergence.p15)) ~ rep(specimens.per.species, times=species))

#Calculate upper tolerance bounds for each per species sample size,
#using the computational method in page 80 of Meeker et al. (2017)
Xperc <- 0.95 #the percentile for which the upper tolerance bounds are calculated
conf <- 0.999 #the confidence level of the upper tolerance bounds
u.index <- qbinom(conf, species, Xperc) + 1
u.l <- function(x) sort(x)[u.index]
KL.p15.ul <- apply(KLdivergence.p15, 1, u.l)

#Examine upper tolerance bounds for different per species sample sizes
plot(specimens.per.species, KL.p15.ul, type="o")
plot(specimens.per.species, log(KL.p15.ul), type="o")
plot(specimens.per.species, log(KL.p15.ul), type="o", xlim=c(0,20), xaxt="n")
axis(1, specimens.per.species)
abline(v=p+1, lty=3)
abline(h=log(0.1), lty=3)


################################################################################
# 4. Compare upper bounds across sample sizes and dimensionalities.
################################################################################

plot(specimens.per.species, log(KL.p5.ul), ylim=log(range(KL.p5.ul, KL.p10.ul, KL.p15.ul, na.rm=T)), 
	type="o", col="blue", ylab="Kullback-Leibler divergence") #p5
points(specimens.per.species, log(KL.p10.ul), type="o", col="green") #p10
points(specimens.per.species, log(KL.p15.ul), type="o", col="red") #p15
#Absolute threshold
abline(h = log(0.1), lty=3)
#Show minimum per species sample size for which the upper bound is below the absolute threshold of 0.1 
abline(v=specimens.per.species[which(KL.p5.ul<=0.1)[1]], lty=3, col="blue") #p5
specimens.per.species[which(KL.p5.ul<=0.1)[1]] 
specimens.per.species[which(KL.p5.ul<=0.1)[1]]/( 5 + ((5*5)-5)/2 ) #specimens per parameter
abline(v=specimens.per.species[which(KL.p10.ul<=0.1)[1]], lty=3, col="green") #p10
specimens.per.species[which(KL.p10.ul<=0.1)[1]]
specimens.per.species[which(KL.p10.ul<=0.1)[1]]/( 10 + ((10*10)-10)/2 ) #specimens per parameter
abline(v=specimens.per.species[which(KL.p15.ul<=0.1)[1]], lty=3, col="red") #p15
specimens.per.species[which(KL.p15.ul<=0.1)[1]]
specimens.per.species[which(KL.p15.ul<=0.1)[1]]/( 15 + ((15*15)-15)/2 ) #specimens per parameter
legend("topright", c("p = 5", "p = 10", "p = 15"), pch=21, lty=1, col=c("blue", "green", "red"))

plot(specimens.per.species, log(KL.p5.ul), xlim=c(1,20), ylim=log(range(KL.p5.ul, KL.p10.ul, KL.p15.ul, na.rm=T)),
	type="o", col="blue", ylab="Kullback-Leibler divergence") #p5
points(specimens.per.species, log(KL.p10.ul), type="o", col="green") #p10
points(specimens.per.species, log(KL.p15.ul), type="o", col="red") #p15
#Show minimum per species sample size that allows estimation of Kullback-Leibler divergence
abline(v=c(6, 11, 16), lty=3, col=c("blue", "green", "red"))
#Absolute threshold
abline(h = log(0.1), lty=3)
legend("topleft", c("p = 5", "p = 10", "p = 15"), pch=21, lty=1, col=c("blue", "green", "red"))


################################################################################
# 5. Write file with upper bounds for all sample sizes and dimensionalities.
################################################################################

KLdivUB <- cbind(specimens.per.species, KL.p5.ul, KL.p10.ul, KL.p15.ul)
#Set directory to save file
setwd(dir.to.save.file)
write.table(KLdivUB, file=paste("KLdivUB_", format(Sys.time(), "%d%b%Y_%H%M%S"), ".txt", sep=""), sep=",")
