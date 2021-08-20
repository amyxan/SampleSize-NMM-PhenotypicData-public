################################################################################
################################################################################
#
# PURPOSE OF THIS SCRIPT
# This script retrieves or calculates the responses of interest from normal 
# mixture models (NMMs) fitted to (previously acquired) random samples from the
# true phenotypic distribution of simulated species groups. Specifically, for
# each simulated species group it i) retrieves the estimated number of
# phenotypic groups (species), ii) calculates the Goodman and Kruskal tau
# statistic of association between the true and estimated classification of
# specimens into phenotypic groups (species), and iii) calculates all
# Kullback-Leibler divergences between true and estimated phenotypic
# groups. For each replicate, this script saves these responses in a list object.
#
# DATA FILES REQUIRED TO RUN THIS SCRIPT
# Files describing the NMMs fitted to the samples from the simulated species
# groups must be available in a directory having no other files. The files
# containing the corresponding simulated species groups must also be available
# in a separate directory. The fitted NMM files and simulated species groups
# files should reside in subfolders within their respective directories. These
# subfolders should be structured as follows: there should be one subfolder for
# each different number of species simulated per species group. Within each
# subfolder corresponding to a particular number of species per species group,
# there should be one subfolder for each distinct level of the experiment's
# factorial design.
#
# CONTENTS OF THIS SCRIPT
# 1. Preliminaries.
# 2. Obtain/calculate responses of interest from each sample and save responses.
# 3. Calculate the time taken to obtain/calculate reponses.
#
################################################################################
################################################################################

{ #Use brackets around code body so that stop() stops script execution immediately

################################################################################
# 1. Preliminaries.
################################################################################

################################################################################
# 1.1. Load packages and define working directories

#load packages
library(GoodmanKruskal) #package to calculate Goodman and Kruskal tau

#define directories
dir.root.species <- "replace with path to directory of all species objects"
dir.root.NMMs <- "replace with path to directory of all NMM objects"
dir.root.responses <- "replace with path to directory to write all response objects"

################################################################################
# 1.2. Define experimental treatment levels and number of replicates for each.
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

#Check that the possible values for s (number of species sampled per group),
#p (number of phenotypic traits) LL.tbs (lower limit of acceptable true.beta.star
#values), sampling proportion (proportion indicating how specimens are to be
#distributed when sampling from a species group), and n (number of specimens to
#be sampled from a species group) are those for which the code has been
#calibrated.
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

#Enter the specific treatment levels that will be used in sections 2 - 3
#of this script, each vector defined below should be a subset of the full
#range of experimental treatment levels considered in the simulation
#experiment, as defined above.
all.s <- c(8)
all.p <- c(5, 10, 15)
all.LL.tbs <- c(0, 0.1, 0.9)
all.sampling.proportion <- c(0.5, 0.7, 0.9)
all.n <- c(1450)

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
number.replicates <- 181

################################################################################
# 1.3. Create file to update running times for calculations at each treatment
# level

save.running.times <- TRUE

if(save.running.times)
{
	setwd(dir.root.responses)
	running.times.file.name <- paste0("responses.running.times_", 
		number.replicates, "replicates.per.level_", 
		format(Sys.time(), "%d%b%Y_%H%M%S"), "_ran_", 
		round(runif(1, min=1, max=1000000)), ".csv")
	write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE)
}

################################################################################
# 1.4. Function for calculating Kullback-Leibler divergence between two normal 
# distributions.
# Adapted from package "rags2ridges"; original code can be found at
# https://rdrr.io/cran/rags2ridges/src/R/rags2ridges.R#sym-KLdiv

KLdiv <- function(Mtrue, Strue, Mest, Sest) {

	#Determine whether the estimated VCV matrix has any Inf or -Inf values.
	#If it does, return Inf as the KL divergence
	if(any(abs(Sest) == Inf)) {
		return(c(NA, "Inf(s) in fitted VCV"))
	}

	#Solve the estimated VCV matrix, compute the KL divergence, and return
	#the result. Return NA if the matrix is not solvable.
	tryCatch(
		{
			inv.Sest <- solve(Sest)
			KLd <- (sum(diag(inv.Sest%*%Strue))
			+ t(Mest - Mtrue)%*%inv.Sest%*%(Mest - Mtrue)
			- nrow(Strue) - log(det(Strue)) + log(det(Sest)))/2
			return(c(as.numeric(KLd), NA))
		},
		error = function(cond) {
			message <- cond$message
			call <- toString(cond$call)
			return(c(NA, paste0("call: ", call, " message: ", message)))
		}
	)
}


################################################################################
# 2. Obtain/calculate responses of interest from each sample and save responses.
################################################################################

#Begin loop which iterates through each treatment level
for(tl in 1:num.tls)
{

#Store start time for obtaining responses from this treatment level
start.time <- Sys.time()

################################################################################
# 2.1. Define and create directories for a given treatment level.

s <- all.tls$s[tl]

#Define directories for this treatment level
dir.tl.NMMs <- paste0(dir.root.NMMs, "/", "NMMs_s", all.tls$s[tl], "/", "NMMs_", all.tls$tl.id[tl])
dir.tl.species <- paste0(dir.root.species, "/", "species_s", all.tls$s[tl], "/", "species_", all.tls$tl.id[tl])
dir.s.responses <- paste0(dir.root.responses, "/", "responses_s", all.tls$s[tl])
dir.tl.responses <- paste0(dir.s.responses, "/", "responses_", all.tls$tl.id[tl])

#Check that number of species groups and fitted NMM files for this treatment level
#is at least equal to the desired number of replicates
setwd(dir.tl.NMMs)
if(length(dir()) < number.replicates)
	stop(paste0("The number of fitted NMMs for ", all.tls$tl.id[tl], " is less than the specified number of replicates."))
setwd(dir.tl.species)
if(length(dir()) < number.replicates)
	stop(paste0("The number of species groups for ", all.tls$tl.id[tl], " is less than the specified number of replicates."))

#As needed, create folders to store all responses for this treatment level
if(!dir.exists(dir.s.responses))
	dir.create(dir.s.responses)
if(!dir.exists(dir.tl.responses))
	dir.create(dir.tl.responses)

################################################################################
# 2.2. Obtain or calculate all responses for each replicate at the given
# treatment level.

#Begin loop which iterates through all replicates for this treatment level
for(i in 1:number.replicates)
{
	#Set working directory to read fitted NMMs
	setwd(dir.tl.NMMs)

	#Load fitted NMM
	NMM <- eval(as.name(load(dir()[i])))

	#Store estimated number of species
	fitted.num.species <- NMM$NMM.fit$G

	#Calculate and store Goodman-Kruskal tau statistic for classification of specimens
	true.classification <- rep(1:s, NMM$sample.sizes)
	fitted.classification <- NMM$NMM.fit$classification
	GKtau.species <- GKtau(fitted.classification, true.classification)

	#Obtain NMM's means and VCV matrices
	all.fitted.means <- NMM$NMM.fit$parameters$mean
	all.fitted.VCVs <- NMM$NMM.fit$parameters$variance$sigma

	#Set working directory to read species groups
	setwd(dir.tl.species)

	#Load species group
	species <- eval(as.name(load(dir()[i])))

	#Ensure that this species group corresponds to the loaded NMM
	if(!identical(species$species.group.tag, NMM$species.group.tag))
		stop("Loaded species group does not match loaded NMM.")

	#Create matrix which will store KL divergence values between all true and estimated distributions
	KLdiv.matrix <- matrix(ncol=s, nrow=fitted.num.species)
	colnames(KLdiv.matrix) <- paste("true", NMM$sample.species.numbers, sep=".")
	rownames(KLdiv.matrix) <- paste("fit", 1:fitted.num.species, sep=".")

	#Create matrix which will store notes related to possible errors in KL divergence computation
	KLdiv.messages <- KLdiv.matrix

	#Obtain species group's means and VCV matrices
	all.true.means <- species$mean.vector
	all.true.VCVs <- species$covariance.matrix

	#Compute and store KL divergence values between all true and estimated distributions
	for(fit.num in 1:fitted.num.species)
	{
		row.name <- paste("fit", fit.num, sep=".")
		for(true.num in NMM$sample.species.numbers)
		{
			col.name <- paste("true", true.num, sep=".")
			KLdiv.result <- KLdiv(all.true.means[[true.num]], all.true.VCVs[[true.num]], all.fitted.means[,fit.num], all.fitted.VCVs[,,fit.num])
			KLdiv.matrix[row.name,col.name] <- KLdiv.result[1]
			KLdiv.messages[row.name,col.name] <- KLdiv.result[2]
		}
	}

	#Store fitted NMM's variance model name
	fitted.model.name <- NMM$NMM.fit$parameters$variance$modelName

	#Store number of individuals sampled per species
	sample.sizes <- matrix(NMM$sample.sizes, nrow=1)
	colnames(sample.sizes) <- colnames(KLdiv.matrix)

################################################################################
# 2.3. Save response.

	#Set working directory to folder to save responses for this treatment level
	setwd(dir.tl.responses)

	#Save responses as a list
	assign(paste0("responses.", NMM$species.group.tag),
		list(species.group.tag=NMM$species.group.tag,
		sample.sizes=sample.sizes,
		fitted.num.species=fitted.num.species, GKtau=GKtau.species,
		KLdiv.matrix=KLdiv.matrix, KLdiv.messages=KLdiv.messages,
		fitted.model.name=fitted.model.name))
	save(list=paste0("responses.", NMM$species.group.tag), 
		file=paste0("responses.", NMM$species.group.tag, ".R"))

} #End loop which iterates through all replicates for a treatment level


################################################################################
# 3. Calculate the time taken to obtain/calculate reponses.
################################################################################

all.tls$time.in.minutes[tl] <- difftime(Sys.time(), start.time, units="mins")

if(save.running.times){
	setwd(dir.root.responses)
	try(
		write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE),
		silent=TRUE
	)
}

} #End loop which iterates through all treatment levels
} #Close brackets around code body
