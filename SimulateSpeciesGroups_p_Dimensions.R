################################################################################
################################################################################
#
# PURPOSE OF THIS SCRIPT
# This script simulates the phenotypes of groups of species as mixtures of
# multivariate normal distributions (McLachlan and Peel. 2000. Finite mixture
# models. Hoboken, NJ: John Wiley and Sons [Series in probability and
# statistics]). In particular, each species in a group is represented by a
# normal distribution of p phenotypic traits, with means for each trait
# sampled from the uniform distribution between -1 and 1, variances for each
# trait sampled from the uniform distribution between 1e-10 and 1 (depending
# on the parameters, the variance of some traits is shrinked when multiplied
# by a value sampled from the uniform distribution between 0.1 and 0.5), and
# correlation coefficients between phenotypic traits sampled from a uniform
# distribution between -1 and 1. The phenotypic overlap between every pair of
# simulated species in a group is measured by "true.beta.star", defined as the
# proportion of individuals of each species contained within the most inclusive
# pair of highest density regions that are i) non-overlapping and ii) equal in
# coverage. For a useful discussion of highest density regions see:
# Hyndman (1996. Computing and graphing highest density regions. The American
# Statistician 50: 120-126).
#
# DATA FILES REQUIRED TO RUN THIS SCRIPT
# None.
#
# USAGE NOTES
# 
# 
# CONTENTS OF THIS SCRIPT
# 1. Preliminaries.
# 2. Define simulation conditions.
# 3. Simulate species groups.
# 4. Calculate the time taken to complete the simulation.
#
################################################################################
################################################################################

{ #Use brackets around code body so that stop() stops script execution immediately

################################################################################
# 1. Preliminaries
################################################################################

################################################################################
# 1.1. Load packages and define working directory

#Load packages
library(ellipse)
library(mvtnorm)
library(clusterGeneration)

#Define the working directory where all simulated species groups will be saved
#note: this script will, if needed, create subdirectories within root.dir to
#store species groups from which the same number of species will later be sampled.
#Within these subdirectories, this script will, if needed, create subdirectories
#to store species groups assigned to the same experimental treatment level.
root.dir <- "replace with path to directory to write all species objects"

################################################################################
# 1.2. Define the number of species in each species group

group.size <- 8

################################################################################
# 1.3. Define experimental treatment levels and number of replicates for each.
# Each treatment level is defined by five variables:
# 	s - number of species that will be sampled from the species group
#	p - number of phenotypic traits characterizing each species
#	LL.tbs - lower limit of acceptable true.beta.star values
#	sampling.proportion - value used to determine the proportion of specimens
#		that will be sampled from each species in the species group
#	n - total number of specimens that will be sampled from the species group
# Note: the values of s, sampling.proportion, and n do not affect the
#	simulation of species groups by this script. However, they affect the
#	sampling of species groups in scripts used subsequently. These values
#	should be defined here so that they can be included in the file name of
#	each simulated species group. That way, species groups are assigned to
#	different experimental treatment levels from the time they are created.

#For each experimental variable, define all possible levels.
#These values describe all possible experimental treatment levels and should
#be the same for every run of this script.
s.options <- c(2, 4, 8)
p.options <- c(5, 10, 15)
LL.tbs.options <- c(0, 0.1, 0.9)
sampling.proportion.options <- c(0.5, 0.7, 0.9)
n.options <- seq(50, 1450, 100)

#Check that the possible values for s (number of species to be sampled from a
#species group in a subsequent script), p (number of phenotypic traits), LL.tbs
#(lower limit of acceptable true.beta.star values), sampling proportion (proportion
#indicating how specimens are to be distributed when sampling from a species group
#in a subsequent script), and n (number of specimens to be sampled from a species
#group in a subsequent script) are those for which the code has been callibrated.
if(!identical(s.options, c(2, 4, 8)))
	stop("The code in the sampling script, which will be used subsequently, is not calibrated to work with s values other than 2, 4 and 8.")
if(!identical(p.options, c(5, 10, 15)))
	stop("The code is not calibrated to work with p values other than 5, 10 and 15.")
if(!identical(LL.tbs.options, c(0, 0.1, 0.9)))
	stop("The code is not calibrated to work with LL.tbs values other than 0, 0.1 and 0.9.")
if(!identical(sampling.proportion.options, c(0.5, 0.7, 0.9)))
	stop("The code in the sampling script, which will be used subsequently, is not calibrated to work with sampling.proportion values other than 0.5, 0.7, and 0.9.")
if(any(n.options%%50 != 0 | n.options < 0))
	stop("The code in the sampling script, which will be used subsequently, is not calibrated to work with n values other than positive multiples of 50.")

#The following facilitates running Sections 2, 3, and 4 of this script at
#specific treatment levels represented by all combinations of user-defined 
#values for each variable:
# Each vector defined below should be a subset of its corresponding levels
# defined above.
all.s <- c(8)
all.p <- c(5, 10)
all.LL.tbs <- c(0)
all.sampling.proportion <- c(0.5)
all.n <- c(50, 150, 250)

#Create dataframe to hold all treatment levels
num.tls <- prod(lengths(list(all.s, all.p, all.LL.tbs, all.sampling.proportion, all.n)))
all.tls <- data.frame(s=numeric(num.tls), p=numeric(num.tls), LL.tbs=numeric(num.tls), sampling.proportion=numeric(num.tls), n=numeric(num.tls), tl.id=character(num.tls), time.in.minutes=numeric(num.tls), stringsAsFactors=FALSE)
if(num.tls > 0)
	all.tls[,1:5] <- expand.grid(all.s, all.p, all.LL.tbs, all.sampling.proportion, all.n)

#Edit selected treatment levels/add additional treatment levels manually
all.tls <- edit(all.tls)

#Update num.tls
num.tls <- nrow(all.tls)

#Check that all supplied values are appropriate
if(!all(all.tls$s %in% s.options))
	stop("At least one of the added treatment levels contains an inappropriate value for s.")
if(!all(all.tls$p %in% p.options))
	stop("At least one of the added treatment levels contains an inappropriate value for p.")
if(!all(all.tls$LL.tbs %in% LL.tbs.options))
	stop("At least one of the added treatment levels contains an inappropriate value for LL.tbs.")
if(!all(all.tls$sampling.proportion %in% sampling.proportion.options))
	stop("At least one of the added treatment levels contains an inappropriate value for sampling.proportion.")
if(!all(all.tls$n %in% n.options))
	stop("At least one of the added treatment levels contains an inappropriate value for n.")

#Add treatment ID (a string representing the treatment level) to a column
#of all.tls
for(i in 1:num.tls)
{
	s.indicator <- paste0("s", all.tls$s[i])
	p.indicator <- paste0("p", paste0(rep(0, times = 2-nchar(as.character(all.tls$p[i]))), collapse=""), all.tls$p[i])
	if(all.tls$LL.tbs[i] == 0){ tbs.indicator <- "B0001"
	} else if(all.tls$LL.tbs[i] == 0.1){ tbs.indicator <- "B0109"
	} else if(all.tls$LL.tbs[i] == 0.9) tbs.indicator <- "B0910"
	sampling.proportion.indicator <- paste0("mp0", round(all.tls$sampling.proportion[i]*10))
	n.indicator <- paste0("n", paste0(rep(0, times = 4-nchar(as.character(all.tls$n[i]))), collapse=""), all.tls$n[i])
	
	#Store string with format "s#p##B####mp##n####"
	all.tls$tl.id[i] <- paste0(s.indicator, p.indicator, tbs.indicator, sampling.proportion.indicator, n.indicator)
}

#Define number of replicates for each treatment level
number.of.species.groups <- 181

################################################################################
# 1.4. Define the maximum number of searches to be done to add a new species
# to a particular species group. Once that number of unsuccessful searches is
# reached, the particular species group that was being built is discarded, and
# the simulation of a new species group begins.

counter.max <- 600

################################################################################
# 1.5. Choose how to display graphically the progress of the simulation 

#Determine if the simulation parameters should be printed
print.simulation.parameters <- TRUE

#Determine if the script will produce figures (produce.figures <- TRUE) 
#or not (produce.figures <- FALSE) 
produce.figures <- FALSE

#Select colors to display different species in phenotypic space
if(produce.figures) species.color <- c("red", "blue", "yellow", "green", "deeppink", "orange", "plum1", "cyan")[1:group.size] #only for group.size < 9
#Check color selection using a pie chart
if(produce.figures) pie(rep(1, group.size), labels=paste("species", 1:group.size), col=species.color, border=species.color, radius=1, clockwise=T)
#Check color selection using a figure legend
if(produce.figures) plot(0:1, 0:1, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n")
if(produce.figures) legend("topright", paste("species", 1:group.size), col=species.color, pch=19, lty=1)

#Select two phenotypic traits (out of p pheotypic traits) that will be used to
#display the phenotypes of species in each group; make sure these values do not
#exceed p; the choice of these two phenotypic traits does not influence
#the simulations, it only influences how the phenotypes of species groups are
#shown in figures:

if(produce.figures==F)
{
	T1 <- NULL
	T2 <- NULL
}

if(produce.figures==T)
{
	T1 <- 1
	T2 <- 2
}

################################################################################
# 1.6. Define a function to estimate true.beta.star between two species

#The function is named "tbs.ridgeline" and estimates true.beta.star between two
#species, say species A and B, with the following arguments:
#A.M, the vector of means of species A
#B.M, the vector of means of species B
#A.VCV, the variance-covariance matrix of species A 
#B.VCV, the variance-covariance matrix of species B

tbs.ridgeline <- function(A.M, B.M, A.VCV, B.VCV)
{
	#Check validity of argument values.
	if(!identical(length(A.M), dim(A.VCV)[1], dim(A.VCV)[2])) stop("A.M and A.VCV are non-conformable arguments")
	if(!identical(length(B.M), dim(B.VCV)[1], dim(B.VCV)[2])) stop("B.M and B.VCV are non-conformable arguments")
	if(!identical(length(A.M), length(B.M))) stop("A.M and B.M have different lengths")

	#Calculate ridgeline manifold.
	#The code for the ridgeline manifold is based on equation 4 in page 2045 of
	#Ray and Lindsay, 2005. "The topography of multivariate normal mixtures" The
	#Annals of Statistics 33: 2042-2065. the ridgeline manifold is evaluated at
	#various values of alpha and the resulting coordinates (in phenotypic space)
	#are captured in matrix ridgeline.coor.
	alpha <- seq(0,1,0.001)
	p.dimensions <- length(A.M)
	ridgeline.coor <- matrix(NA, nrow=length(alpha), ncol=p.dimensions)
	for (i in 1:length(alpha))
	{
		a <- solve((1-alpha[i])*solve(A.VCV)  +  alpha[i]*solve(B.VCV))
		b <- (1-alpha[i])*solve(A.VCV)%*%A.M  +  alpha[i]*solve(B.VCV)%*%B.M
		d <- a %*% b
		ridgeline.coor[i,] <- d
	}

	#Calculate true.beta.star
	md.A <- mahalanobis(ridgeline.coor, A.M, A.VCV)
	md.B <- mahalanobis(ridgeline.coor, B.M, B.VCV)
	true.beta.vector.A <- pchisq(md.A, p.dimensions, lower.tail = TRUE, log.p = FALSE)
	true.beta.vector.B <- pchisq(md.B, p.dimensions, lower.tail = TRUE, log.p = FALSE)
	beta.dif <- abs(true.beta.vector.A - true.beta.vector.B)
	o.beta.dif <- order(beta.dif)
	true.beta.star <- mean(c(true.beta.vector.A[o.beta.dif[1]], true.beta.vector.B[o.beta.dif[1]]))

	#Return a list with the estimate of true.beta.star and a matrix with alpha
	#values, the respective phenotypic coordinates of the ridgeline manifold,
	#and the respective proportion of the distribution of each species contained
	#within (hyper)elliptical regions.
	ridgeline <- cbind(alpha, ridgeline.coor, true.beta.vector.A, true.beta.vector.B)
	colnames(ridgeline) <- c("alpha", paste("ridgeline.coor", 1:p.dimensions, sep="."), "true.beta.vector.A", "true.beta.vector.B")	
	return(list(true.beta.star = true.beta.star, ridgeline = ridgeline))
}

################################################################################
# 1.7. Create file to update running times for simulations at each treatment
# level

save.running.times <- TRUE

if(save.running.times)
{
	setwd(root.dir)
	running.times.file.name <- paste0("simulation.running.times_", 
		number.of.species.groups, "replicates.per.level_", 
		format(Sys.time(), "%d%b%Y_%H%M%S"), "_ran_", 
		round(runif(1, min=1, max=1000000)), ".csv")
	write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE)
}


################################################################################
# 2. Define simulation conditions.
################################################################################

#Loop allows Sections 2, 3, and 4 of this script to execute once for each
#treatment level defined in Section 1
for(tl in 1:num.tls)
{

################################################################################
# 2.1. Create subdirectory to store species groups simulated at the same
# treatment level

#Obtain treatment level ID for this treatment level
tl.id <- all.tls$tl.id[tl]

#Set working directory to the folder where subfolders of new species groups will
#be created
setwd(root.dir)

#Create and set working directory to the folder which will store all simulated
#species groups with the current treatment level's level of s
s.dir <- paste0("species_s", all.tls$s[tl])
if(!dir.exists(s.dir))
	dir.create(s.dir)
setwd(s.dir)

#Create and set working directory to the folder which will store simulated
#species groups for this treatment level
tl.dir <- paste("species", tl.id, sep="_")
if(!dir.exists(tl.dir))
	dir.create(tl.dir)
setwd(tl.dir)

################################################################################
# 2.2. Obtain values for each experimental variable at this treatment level 

#Define the number of species that will be sampled from each species group,
#note that s should be <= group.size
s <- all.tls$s[tl]

#Define the number of phenotypic traits characterizing each species
p <- all.tls$p[tl]

#Define lower limit of acceptable true.beta.star values
LL.tbs <- all.tls$LL.tbs[tl]

#Define upper limit of acceptable true.beta.star values,
#choose the value corresponding to the previous choice:
if(LL.tbs == 0.9) UL.tbs <- 1
if(LL.tbs == 0) UL.tbs <- 0.1
if(LL.tbs == 0.1) UL.tbs <- 0.9

#Define the maximum proportion of the sample of each of these species groups
#that will be taken from half of the sampled species
sampling.proportion <- all.tls$sampling.proportion[tl]

#Define the number of specimens that will be sampled from each species group
n <- all.tls$n[tl]

################################################################################
# 2.3. Define parameters for the search of mean phenotypic values of species:
# when adding new species to a group of species, should the search for the mean
# phenotypic values of new species be biased? Bias here means that the search
# focuses within a particular range of values of the probability density of the
# phenotypes of the species already in the species group.

#To simulate species groups with true.beta.star >= 0.9, or 0.1 < true.beta.star 
#< 0.9, searches should be biased towards low values of the probability density 
#of the phenotypes of the species already in the species group:

if(LL.tbs == 0){ bias.mean.to.low.pdf <- FALSE
} else bias.mean.to.low.pdf <- TRUE

#To simulate species groups with true.beta.star <= 0.1, searches should be biased
#towards high values of the probability density of the phenotypes of the species
#already in the species group:

if(LL.tbs == 0){ bias.mean.to.high.pdf <- TRUE
} else bias.mean.to.high.pdf <- FALSE

#Set the upper limit of the value of the probability density beyond which
#no search is done:

if(bias.mean.to.low.pdf==F) mixture.density.meanspace.UL <- c()
if(p==5 & bias.mean.to.low.pdf) mixture.density.meanspace.UL <- 1e-10 
if(p==10 & bias.mean.to.low.pdf) mixture.density.meanspace.UL <- 1e-100
if(p==15 & bias.mean.to.low.pdf) mixture.density.meanspace.UL <- 1e-200 

if(bias.mean.to.high.pdf==F) mixture.density.meanspace.LL <- c()
if(p==5 & bias.mean.to.high.pdf) mixture.density.meanspace.LL <- 0.01 
if(p==10 & bias.mean.to.high.pdf) mixture.density.meanspace.LL <- 0.001 
if(p==15 & bias.mean.to.high.pdf) mixture.density.meanspace.LL <- 0.00001 

################################################################################
# 2.4. Define parameters that determine if and how the phenotypic variance of new
# species to be added to a species group should be shrinked (reduced), as the
# number of unsuccesful searches increases. These parameters determine the number
# of traits whose variance will be shrinked. To shrink the variance of a given trait,
# the original value (sampled from a uniform distriution between 1e-10 and 1) is
# multiplied by a value sampled from the uniform distribution between 0.1 and 0.5.

#To simulate species groups with true.beta.star >= 0.9, the phenotypic variance
#of species to be added to a species group should be shrinked (reduced) as the
#number of unsuccesful searches increases:

if(LL.tbs == 0.9){ reduce.variance.with.tries <- TRUE
} else reduce.variance.with.tries <- FALSE

#Set the parameters that determine how the number of traits with shrinked variance
#increases as the number of unsuccessful searches increases:

if(reduce.variance.with.tries==F)
{
	b1 <- c()
	b2 <- c()
}

if(reduce.variance.with.tries==T)
{
	x <- 1:counter.max
	b1 <- 200 
	b2 <- 2
	if(produce.figures) plot(x, round((p*x^b2)/(b1+x^b2),0), xlab="Unsuccessful serches (Counter)", ylab="Number of traits with shrinked variance", type="l", lwd=2)
	if(produce.figures) abline(v = seq(0, max(x), 10), lty=3, col="gray80")
}

################################################################################
# 2.5. Check the choices of simulation parameters

if(print.simulation.parameters)
{
	#Comment or uncomment the lines below as desired

	#print(paste("LL.tbs =", LL.tbs))
	#print(paste("UL.tbs =", UL.tbs))
	print(paste("group.size =", group.size))
	#print(paste("p =", p))
	#print(paste("s =", s))
	#print(paste("sampling.proportion =", sampling.proportion))
	print(paste("produce.figures =", produce.figures))
	#print(paste("T1 =", T1))
	#print(paste("T2 =", T2))
	print(paste("counter.max =", counter.max))
	print(paste("number.of.species.groups =", number.of.species.groups))
	print(paste("save.running.times =", save.running.times))
	#print(paste("bias.mean.to.low.pdf =", bias.mean.to.low.pdf))
	#print(paste("mixture.density.meanspace.UL =", mixture.density.meanspace.UL))
	#print(paste("bias.mean.to.high.pdf =", bias.mean.to.high.pdf))
	#print(paste("mixture.density.meanspace.LL =", mixture.density.meanspace.LL))
	#print(paste("reduce.variance.with.tries =", reduce.variance.with.tries))
	#print(paste("b1 =", b1))
	#print(paste("b2 =", b2))
	print(paste("root directory: ", root.dir))
}


################################################################################
# 3. Simulate species groups
################################################################################

start.time <- Sys.time() #Store simulation starting time

#Store all combinations of species pairs in a vector
index.species.pairs <- combn(1:group.size, 2)

N <- 1
while(N <= number.of.species.groups)
{
	#Generate vector with names for the mean vectors and variance-covariance matrices for each species in the group
	mean.vector.names <- paste("M.", 1:group.size, sep="")
	covariance.matrix.names <- paste("VCV.", 1:group.size, sep="")
	#Randomly generate the the mean vector and variance-covariance matrix for the first species in the group 
	assign(mean.vector.names[1], runif(p, -1, 1))
	assign(covariance.matrix.names[1], genPositiveDefMat(covMethod="unifcorrmat", rangeVar=c(1e-10,1), dim=p)$Sigma)
	counter <- 0

	#Declare vector which will store true beta star values between all species pairs
	all.pairs.true.beta.star <- rep(NA, times = ncol(index.species.pairs))

	for(i in 2:group.size)
	{
		if(counter > counter.max) break
		counter <- 0
		min.true.beta.star <- LL.tbs - 1 
		max.true.beta.star <- UL.tbs + 1
		while(min.true.beta.star < LL.tbs | max.true.beta.star > UL.tbs)
		{
			counter <- counter + 1
			if(counter > counter.max) break

			if(bias.mean.to.low.pdf==F & bias.mean.to.high.pdf==F) assign(mean.vector.names[i], array(runif(p, -1, 1), dim=c(p,1)))

			if(bias.mean.to.low.pdf)
			{
				mean.space.ran <- array(runif(805260, -1, 1), dim=c(805260/p,p))
				density.meanspace <- matrix(NA, nrow=nrow(mean.space.ran), ncol=(i-1))
				for(g in 1:(i-1))
				{
					density.meanspace[,g] <- dmvnorm(mean.space.ran, mean = eval(as.name(mean.vector.names[g])), sigma = eval(as.name(covariance.matrix.names[g])))
				}
				mixture.density.meanspace <- rowMeans(density.meanspace)
				probability.sample.meanspace <- rep(0, times=length(mixture.density.meanspace))
				if(min(mixture.density.meanspace) < mixture.density.meanspace.UL){
					probability.sample.meanspace[mixture.density.meanspace < mixture.density.meanspace.UL] <- 1
				} else {
					probability.sample.meanspace[order(mixture.density.meanspace)[1:5]] <- 1
				}
				M <-	as.vector(mean.space.ran[sample(1:nrow(mean.space.ran), size=1, prob=probability.sample.meanspace),])
				assign(mean.vector.names[i],M)
			}

			if(bias.mean.to.high.pdf)
			{
				mean.space.ran <- array(runif(805260, -1, 1), dim=c(805260/p,p))
				density.meanspace <- matrix(NA, nrow=nrow(mean.space.ran), ncol=(i-1))
				for(g in 1:(i-1))
				{
					density.meanspace[,g] <- dmvnorm(mean.space.ran, mean = eval(as.name(mean.vector.names[g])), sigma = eval(as.name(covariance.matrix.names[g])))
				}
				mixture.density.meanspace <- rowMeans(density.meanspace)
				probability.sample.meanspace <- rep(0, times=length(mixture.density.meanspace))
				if(max(mixture.density.meanspace) > mixture.density.meanspace.LL){
					probability.sample.meanspace[mixture.density.meanspace > mixture.density.meanspace.LL] <- 1
				} else {
					probability.sample.meanspace[order(mixture.density.meanspace)[(length(mixture.density.meanspace)-5):length(mixture.density.meanspace)]] <- 1
				}
				M <-	as.vector(mean.space.ran[sample(1:nrow(mean.space.ran), size=1, prob=probability.sample.meanspace),])
				assign(mean.vector.names[i],M)
			}

			if(reduce.variance.with.tries == FALSE) assign(covariance.matrix.names[i], genPositiveDefMat(covMethod="unifcorrmat", rangeVar=c(1e-10,1), dim=p)$Sigma)

			if(reduce.variance.with.tries == TRUE)
			{	
				V <- genPositiveDefMat(covMethod="unifcorrmat", rangeVar=c(1e-10,1), dim=p)$Sigma
				number.traits.to.shrink <- round((p*counter^b2)/(b1+counter^b2),0)
				ad <- rep(1, p)
				ad[sample(1:p, number.traits.to.shrink)] <- runif(number.traits.to.shrink, 0.1, 0.5)
				A <- diag(ad, p)
				V <- A%*%V%*%t(A)
				assign(covariance.matrix.names[i], V)
			}

			true.beta.star <- rep(NA, times=(i-1))
			for(j in 1:(i-1))
			{
				true.beta.star[j] <- tbs.ridgeline(eval(as.name(mean.vector.names[i])), eval(as.name(mean.vector.names[j])), eval(as.name(covariance.matrix.names[i])), eval(as.name(covariance.matrix.names[j])))$true.beta.star
			}
			min.true.beta.star <- min(true.beta.star, na.rm=T)
			max.true.beta.star <- max(true.beta.star, na.rm=T)
			
			if(produce.figures)
			{
				if(length(dev.list())>1) dev.set(2)
				plot(ellipse(eval(as.name(covariance.matrix.names[1]))[c(T1, T2), c(T1, T2)], centre=eval(as.name(mean.vector.names[1]))[c(T1, T2)]),
					xlim=c(-3,3), ylim=c(-3,3), xlab=paste("Phenotypic axis ", T1), ylab=paste("Phenotypic axis ", T2),
					cex.lab=1, cex.axis=1, type="l", col=species.color[1], main=paste("p = ", p, ", species group = ", N, ", species = ", i, ", counter = ", counter))
				points(eval(as.name(mean.vector.names[1]))[T1], eval(as.name(mean.vector.names[1]))[T2], col=species.color[1], pch=19, cex=1.5)
				for(h in 1:i)
				{
					points(ellipse(eval(as.name(covariance.matrix.names[h]))[c(T1, T2), c(T1, T2)], centre=eval(as.name(mean.vector.names[h]))[c(T1, T2)]), type="l", col=species.color[h])
					points(eval(as.name(mean.vector.names[h]))[T1], eval(as.name(mean.vector.names[h]))[T2], col=species.color[h], pch=19, cex=1.5)
				}
				#Sys.sleep(5) #use this line to suspend execution for a few seconds and thus have more time to examine figures
				polygon(c(-1,-1,1,1), c(-1,1,1,-1), border="black")
				legend("topright", paste("B.star range :", round(min.true.beta.star,2), "-", round(max.true.beta.star,2)))
				if(length(dev.list())<2) dev.new() else dev.set(3)
				if(counter <2) plot(c(1,counter.max), c(0,1), type="n", bty="n", xlab="Counter", ylab="Beta-star")
				if(counter <2) abline(h=c(0.5, LL.tbs, UL.tbs), lty=3)
				points(counter, min.true.beta.star, col="blue", pch=21)
				points(counter, max.true.beta.star, col="red", pch=21)
				if(LL.tbs > 0 & min.true.beta.star > LL.tbs) points(counter, min.true.beta.star, col="blue", pch=19)
				if(UL.tbs < 1 & max.true.beta.star < UL.tbs) points(counter, max.true.beta.star, col="red", pch=19)				
			}
		}
		#If the true.beta.star between this species and all others in the species group so far
		#falls within the required true.beta.star limits, then add these true.beta.star values
		#to the vector storing true.beta.star values between all species pairs
		add.species.to.group <- "no"
		if(LL.tbs == 0.1 & UL.tbs == 0.9 & min.true.beta.star > LL.tbs & max.true.beta.star < UL.tbs) add.species.to.group <- "yes"
		if(LL.tbs == 0.9 & UL.tbs == 1 & min.true.beta.star >= LL.tbs) add.species.to.group <- "yes"
		if(LL.tbs == 0 & UL.tbs == 0.1 & max.true.beta.star <= UL.tbs) add.species.to.group <- "yes"
		if(add.species.to.group == "yes")
		{
			cols.in.index.species.pairs <- which(index.species.pairs[2,] == i)
			all.pairs.true.beta.star[cols.in.index.species.pairs] <- true.beta.star
		}
	}
	#If s species have been selected, save this species group
	if(length(which(is.na(all.pairs.true.beta.star))) == 0)
	{
		#Save parameters for the species group, create species group index using treatment level indicators, current time, date and a random number
		species.group.tag <- paste("species_", tl.id , "_", format(Sys.time(), "%d%b%Y_%H%M%S"), "_ran_", round(runif(1, min=1, max=1000000)), sep="") 
		mean.vectors <- tapply(mean.vector.names, 1:length(mean.vector.names), FUN = function(x) eval(as.name(x)))
		covariance.matrices <- tapply(covariance.matrix.names, 1:length(covariance.matrix.names), FUN = function(x) eval(as.name(x)))
		#Save, in an R object of class list, the parameters describing each simulated species group
		assign(paste("parameters.", species.group.tag, sep=""),
			list(species.group.tag=species.group.tag,
			mean.vector=mean.vectors,
			covariance.matrix=covariance.matrices,
			B.star=all.pairs.true.beta.star))
		#Check the R object of class list that has the parameters describing the simulated species pair
		#eval(as.name(paste("parameters.", species.group.tag, sep="")))
		save(list=paste("parameters.", species.group.tag, sep=""), file=paste("parameters.", species.group.tag, ".R", sep=""))
		#Remove the R object of class list that has the parameters describing the simulated species pair
		rm(list=c(paste("parameters.", species.group.tag, sep="")))
		N <- N+1
	}
}


################################################################################
# 4. Calculate the time taken to complete the simulation 
################################################################################

all.tls$time.in.minutes[tl] <- difftime(Sys.time(), start.time, units="mins")

if(save.running.times){
	setwd(root.dir)
	try(
		write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE),
		silent=TRUE
	)
}

} #End of for loop which iterated through each treatment level

} #Close brackets around code body
