################################################################################
################################################################################
#
# PURPOSE OF THIS SCRIPT
# This script simulates the process of sampling specimens from species groups
# with two, four, or eight species and then carries out an analysis of the
# phenotypes of the specimens sampled to determine species limits using normal
# mixture models.
#
# DATA FILES REQUIRED TO RUN THIS SCRIPT
# Files describing the phenotypes of simulated species groups must be available
# in a directory having no other files; the number of these files (species groups)
# should not be less than the number of replicates per treatment level.
#
# CONTENTS OF THIS SCRIPT
# 1. Preliminaries.
# 2. Define specific treatment level and directories.
# 3. Sample species groups.
# 4. Fit normal mixture models.
# 5. Create object (list) with the fitted normal mixture models and
#    information on treatment levels.
# 6. Save results.
# 7. Calculate the time taken to sample and fit normal mixture models to all
#    species groups. 
#
################################################################################
################################################################################

{ #Use brackets around code body so that stop() stops script execution immediately

################################################################################
# 1. Preliminaries.
################################################################################

################################################################################
# 1.1. Load packages and setup options for Mclust.

#Load packages
library(mvtnorm) #package to sample multivariate normal distributions
library(mclust)  #package to fit normal mixture models

#Check and set Mclust options as needed (see help for function "Mclust")
mclust.options() #check Mclust options

#Save default options (to restore defaults later, if desired)
#OptMc <- mclust.options()

#Choosing value "VARS" for argument "hcUse" uses the original variables in the
#hierarchical agglomerative clustering that is used in the initialization of
#the expectation maximization algorithm implemented in function "Mclust"
mclust.options(hcUse="VARS")

################################################################################
# 1.2. Define working directories.

#Directory which contains all simulated species groups
root.dir.to.read.species.groups <- "replace with path to directory of all species objects"

#Directory to which all fitted NMMs will be written
root.dir.to.write.NMMs <- "replace with path to directory to write all NMM objects"

################################################################################
# 1.3. Define experimental treatment levels and number of replicates for each.
# Each treatment level is defined by five variables:
# 	s - number of species that will be sampled from the species group
#	p - number of phenotypic traits characterizing each species
#	LL.tbs - lower limit of acceptable true.beta.star values
#	sampling.proportion - value used to determine the proportion of specimens
#		that will be sampled from each species in the species group
#	n - total number of specimens that will be sampled from the species group

#For each experimental variable, define the levels for which the code
#is calibrated; these values describe the full range of experimental treatment
#levels considered in the simulation eperiment and should be the same for every
#run of this script.
s.options <- c(2, 4, 8)
p.options <- c(5, 10, 15)
LL.tbs.options <- c(0, 0.1, 0.9)
sampling.proportion.options <- c(0.5, 0.7, 0.9)
n.options <- seq(50, 1450, 100)

#Check that the full range of values for s (number of species sampled per
#group), p (number of phenotypic traits), LL.tbs (lower limit of
#acceptable true.beta.star values), sampling proportion (proportion indicating
#how specimens are to be distributed when sampling from a species group), and
#n (number of specimens to be sampled from a species group) are those for
#which the code is calibrated.
if(!identical(s.options, c(2, 4, 8)))
	stop("The code is not calibrated to work with s values other than 2, 4 and 8.")
if(!identical(p.options, c(5, 10, 15)))
	stop("The code is not calibrated to work with p values other than 5, 10 and 15.")
if(!identical(LL.tbs.options, c(0, 0.1, 0.9)))
	stop("The code is not calibrated to work with LL.tbs values other than 0, 0.1 and 0.9.")
if(!identical(sampling.proportion.options, c(0.5, 0.7, 0.9)))
	stop("The code is not calibrated to work with sampling.proportion values other than 0.5, 0.7, and 0.9.")
if(any(n.options%%50 != 0 | n.options < 0))
	stop("The code is not calibrated to work with n values other than positive multiples of 50.")

#Enter the specific treatment levels that will be used in sections 2 - 7
#of this script, each vector defined below should be a subset of the full
#range of experimental treatment levels considered in the simulation
#experiment, as defined above.
all.s <- c(2, 8)
all.p <- c(5, 10, 15)
all.LL.tbs <- c(0, 0.9)
all.sampling.proportion <- c(0.5, 0.7, 0.9)
all.n <- c(50, 1450)

#Create dataframe to hold all treatment levels
num.tls <- prod(lengths(list(all.s, all.p, all.LL.tbs, all.sampling.proportion, all.n)))
all.tls <- data.frame(s=numeric(num.tls), p=numeric(num.tls), LL.tbs=numeric(num.tls), sampling.proportion=numeric(num.tls), n=numeric(num.tls), tl.id=character(num.tls), time.in.minutes=numeric(num.tls), stringsAsFactors=FALSE)
if(num.tls > 0)
	all.tls[,1:5] <- expand.grid(all.s, all.p, all.LL.tbs, all.sampling.proportion, all.n)

#Edit selected treatment levels/add additional treatment levels manually, as needed
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
number.replicates <- 181

################################################################################
# 1.4. Create file to update running times for simulations at each treatment
# level

save.running.times <- TRUE

if(save.running.times)
{
	setwd(root.dir.to.write.NMMs)
	running.times.file.name <- paste0("sampling.running.times_", 
		number.replicates, "replicates.per.level_", 
		format(Sys.time(), "%d%b%Y_%H%M%S"), "_ran_", 
		round(runif(1, min=1, max=1000000)), ".csv")
	write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE)
}

################################################################################
# 1.5. Define empirical sampling proportion distributions: the distributions
# s=8 and s=4 hard-coded below are computed from the script 
# "GenerateEmpiricalDistributionForSampling.R", where sampling.proportion is the
# value of half.of.species.proportion in the aforementioned script. (The
# distributions for s=2 are simple to compute and therefore not addressed by the
# aforementioned script.)

proportion.distribution <- function(s, sampling.proportion){
	if(s == 8){
		if(sampling.proportion == 0.9) 
			return(c(0.54, 0.20, 0.10, 0.06, 0.04, 0.02, 0.02, 0.02))
		if(sampling.proportion == 0.7) 
			return(c(0.36, 0.16, 0.10, 0.08, 0.08, 0.08, 0.08, 0.06))
		if(sampling.proportion == 0.5){
			return(c(0.14, 0.14, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12))
		}
	}
	if(s == 4){
		if(sampling.proportion == 0.9) 
			return(c(0.72, 0.18, 0.06, 0.04))
		if(sampling.proportion == 0.7) 
			return(c(0.50, 0.20, 0.16, 0.14))
		if(sampling.proportion == 0.5)
			return(c(0.26, 0.26, 0.24, 0.24))
	}
	if(s == 2){
		if(sampling.proportion == 0.9) 
			return(c(0.9, 0.1))
		if(sampling.proportion == 0.7) 
			return(c(0.7, 0.3))
		if(sampling.proportion == 0.5)
			return(c(0.5, 0.5))
	}
	return(NA)
}


################################################################################
# 2. Define specific treatment levels and directories.
################################################################################

#Toop allows the following sections of this script to execute once for each
#treatment level defined in Section 1
for(tl in 1:num.tls)
{

################################################################################
# 2.1. Obtain values for each experimental variable at this treatment level.

#Define the true number of species per species group
s <- all.tls$s[tl]

#Define the number of phenotypic dimensions for the species in a species group
p <- all.tls$p[tl]

#Define the lower limit of true beta star values among the species in a
#species group
B.star.LL <- all.tls$LL.tbs[tl]

#Define the B.star.range, to be saved along with the fitted NMM
if(B.star.LL == 0) B.star.range = "[0.0, 0.1]"
if(B.star.LL == 0.1) B.star.range = "[0.1, 0.9]"
if(B.star.LL == 0.9) B.star.range = "[0.9, 1.0]"

#Define the sampling.proportion level, which determines the proportion of specimens
#to be sampled from different species in the species group
sampling.proportion <- all.tls$sampling.proportion[tl]

#Define the total number of specimens which should be sampled
n <- all.tls$n[tl]

#Determine the number of specimens per species for this sample
sample.sizes <- round(proportion.distribution(s, sampling.proportion) * n)

#Obtain treatment level ID for this treatment level 
tl.id <- all.tls$tl.id[tl]

################################################################################
# 2.2. Create subdirectories to store fitted NMMs for each treatment level, and
# set working directory to folder containing simulated species groups for the
# corresponding treatment level.

#Set working directory to the folder where subfolders of fitted NMMS will
#be created
setwd(root.dir.to.write.NMMs)

#Create and set working directory to the folder which will store all fitted
#NMMs with the current treatment level of s
s.dir <- paste0("NMMs_s", s)
if(!dir.exists(s.dir))
	dir.create(s.dir)
setwd(s.dir)

#Create the folder which will store fitted NMMs for this treatment level
tl.dir <- paste("NMMs", tl.id, sep="_")
if(!dir.exists(tl.dir))
	dir.create(tl.dir)
NMM.tl.dir <- paste(getwd(), tl.dir, sep="/")

#Change working directory to the directory to read species groups
setwd(root.dir.to.read.species.groups)
if(!dir.exists(paste0("species_s", s)))
	stop("There seems to be no directory with simulated species groups for the current treatment level of s.")
setwd(paste0("species_s", s))
if(!dir.exists(paste("species", tl.id, sep="_")))
	stop("There seems to be no directory with simulated species groups for the current treatment level.")
setwd(paste("species", tl.id, sep="_"))
species.tl.dir <- getwd()

#Ensure that all species groups for this treatment level have been created
if(length(dir()) < number.replicates)
	stop(paste("There seem to be less simulated species groups than required by the number of replicates:
		available number of simulated species groups = ", length(dir()),
		", required number of replicates = ", number.replicates, ", treatment ID = ", tl.id, ".", sep=""))
if(length(dir()) > number.replicates)
	warning(paste("There seem to be more simulated species groups than required by the number of replicates:
		available number of simulated species groups = ", length(dir()),
		", required number of replicates = ", number.replicates, ", treatment ID = ", tl.id, ".", sep=""))


################################################################################
# 3. Sample species groups.
################################################################################

start.time <- Sys.time() #Store starting time for the loop below

for(i in 1:number.replicates){

	#Set working directory to read simulated species groups
	setwd(species.tl.dir)

	#Snsure that the species group parameters match the desired treatment levels
	species.group <- load(dir()[i], verbose=T)
	species.group.tag <- eval(as.name(species.group))$species.group.tag
	if(substr(species.group.tag, 9, 27)!=tl.id)
		stop("The treatment level ID of this species group does not match the directory it is in.")	
	if((B.star.LL == 0 & max(eval(as.name(species.group))$B.star)>0.1)
	| (B.star.LL == 0.1 & (max(eval(as.name(species.group))$B.star)>=0.9 | min(eval(as.name(species.group))$B.star)<=0.1))
	| (B.star.LL == 0.9 & min(eval(as.name(species.group))$B.star)<0.9))
		stop(paste("The beta star range for the ", i,
		"th species group in this directory does not agree with its file name.",
		sep=""))
	if(p != length(eval(as.name(species.group))$mean.vector[[1]]))
		stop(paste("The number of phenotypic dimensions for the ", i,
		"th species group in this directory does not agree with its file name.",
		sep=""))
	
	#Initialize empty data frame to store sample data
	sample.data <- matrix(NA, nrow=sum(sample.sizes), ncol=p)
	colnames(sample.data) <- paste("Trait", 1:p, sep=".")

	#Generate s random values between 1 and 8, so the order of species sampling is random
	sample.species.numbers <- sample(1:8, s)

	#Simulate random sample for each species number in sample.species.numbers
	row.index <- 1
	for(j in 1:s){
		mean <- eval(as.name(species.group))$mean.vector[[sample.species.numbers[j]]]
		var <- eval(as.name(species.group))$covariance.matrix[[sample.species.numbers[j]]]
		sample.data[row.index:(row.index+sample.sizes[j]-1),] <- 
			rmvnorm(sample.sizes[j], mean = mean, sigma = var)
		row.index <- row.index + sample.sizes[j]
	}
	remove(species.group)
	
	
################################################################################
# 4. Fit normal mixture models.
################################################################################

	NMM.fit <- Mclust(sample.data, G = 1:12)
	rm(sample.data)
	#summary(NMM.fit)
	
	
################################################################################
# 5. Create object (of class list) with the fitted normal mixture models and
#    information on treatment levels.
################################################################################

	NMM <- list(species.group.tag = species.group.tag,
	s = s, sample.species.numbers = sample.species.numbers,
	sample.sizes = sample.sizes, Beta.star.range = B.star.range,
	sampling.proportion = sampling.proportion,
	NMM.fit = NMM.fit)


################################################################################
# 6. Save results.
################################################################################

	#Set working directory to write files with the normal mixture models
	setwd(NMM.tl.dir)
	save(NMM, file=paste("NMM.", species.group.tag, ".R", sep=""))
	remove(NMM.fit, NMM)
}


################################################################################
# 7. Calculate the time taken to sample and fit normal mixture models to all
#    species groups.
################################################################################

all.tls$time.in.minutes[tl] <- difftime(Sys.time(), start.time, units="mins")

if(save.running.times){
	setwd(root.dir.to.write.NMMs)
	try(
		write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE),
		silent=TRUE
	)
}

} #End of for loop which iterated through each treatment level

} #Close brackets around code body
