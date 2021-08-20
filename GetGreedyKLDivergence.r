################################################################################
################################################################################
# PURPOSE OF THIS SCRIPT
# This script finds and summarizes the KL divergence measure for each true 
# species at each of the specified treatment levels. It does so by examining all
# KL divergences between all possible pairs of true and fitted species'
# phenotypic distributions, and it chooses the pairing that minimizes the KL
# divergence for each species while maintaining no overlap between the species in
# each pair. The results are stored in a table that summarizes the distribution
# of KL divergences assigned to each true species for each treatment level.
# Specifically, this table has one row for each treatment level (defined by the
# number of true species, dimensionality, Beta Star, sampling proportion, and 
# total number of specimens in the sample) and six columns for each true
# species: minimum KL divergence, 1st quartile, median, 3rd quartile, maximum, 
# and number of NAs seen across all replicates. The script should be run with
# one level of s (number of true species) at a time.
#
# DATA FILES REQUIRED TO RUN THIS SCRIPT
# Response objects including all possible Kullback-Leibler divergence measures 
# for a given species group must be available for each replicate in separate 
# directories for each treatment level.
#
# CONTENTS OF THIS SCRIPT
# 1. Preliminaries
# 2. Count number of conforming responses.
# 3. Save table of total counts and running times.
################################################################################
################################################################################

{ #Use brackets around code body so that stop() stops script execution immediately

################################################################################
# 1. Preliminaries
################################################################################

################################################################################
# 1.1. Define working directories

dir.root.responses <- "replace with path to directory of all response objects"
dir.to.write.table <- "replace with path to directory to write KL divergence summary table"

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
#is callibrated; these values describe the full range of experimental treatment
#levels considered in the simulation experiment and should be the same for every
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
#callibrated.
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

#Enter the specific treatment levels that will be used in section 2 of this script,
#each vector defined below should be a subset of the full range of experimental
#treatment levels considered in the simulation experiment, as defined above.
all.s <- c(2)
all.p <- c(5, 10, 15)
all.LL.tbs <- c(0, 0.1, 0.9)
all.sampling.proportion <- c(0.5, 0.7, 0.9)
all.n <- seq(50, 1450, 100)

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
# 1.4. Create file to update running times for simulations at each treatment
# level

save.running.times <- FALSE

if(save.running.times)
{
	setwd(dir.to.write.table)
	running.times.file.name <- paste0("KLdiv.selection.running.times_", 
		number.replicates, "replicates.per.level_", 
		format(Sys.time(), "%d%b%Y_%H%M%S"), "_ran_", 
		round(runif(1, min=1, max=1000000)), ".csv")
	write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE)
}


################################################################################
# 2. Map KL divergences to each true species and summarize results in a table.
################################################################################

#Create data frame to hold KL divergence mesaures for each true species.
KLdiv.summary.table <- all.tls

################################################################################
# 2.1. Set directory to read responses from a particular treatment level, and
# define variables that will aid in summarizing each response.

#Begin loop which iterates through each treatment level
for(tl in 1:num.tls)
{

#Store start time for counting accurate responses from this treatment level
start.time <- Sys.time()

#Define directory containing response files for this treatment level
dir.tl.responses <- paste0(dir.root.responses, "/", "responses_s", all.tls$s[tl], "/", "responses_", all.tls$tl.id[tl])

#Create matrix to hold all selected KL divergences for each true species in
#each replicate
KLdiv.all.reps <- matrix(NA, nrow=number.replicates, ncol=all.tls$s[tl])
colnames(KLdiv.all.reps) <- paste("sampled", 1:all.tls$s[tl], sep=".")

################################################################################
# 2.2. Map KL divergence to each true species and store result in KLdiv.all.reps.

for(i in 1:number.replicates)
{
	#Set working directory to read response file
	setwd(dir.tl.responses)

	#Load response object
	response.obj.name <- load(dir()[i])
	response <- eval(as.name(response.obj.name))
	remove(list=response.obj.name)
	remove(response.obj.name)

	#Store fitted number of species
	G <- response$fitted.num.species

	#Store matrix of all KL divergences
	KLdiv.matrix.trimmed <- response$KLdiv.matrix
	colnames(KLdiv.matrix.trimmed) <- paste("sampled", 1:all.tls$s[tl], sep=".")

	#Obtain minimum and maximum KL divergence via "greedy" selection
	total.distinct.pairs <- min(nrow(KLdiv.matrix.trimmed), ncol(KLdiv.matrix.trimmed))
	num.pairs.selected <- 0
	while(num.pairs.selected<total.distinct.pairs && !all(is.na(KLdiv.matrix.trimmed))){
		#Find index of smallest KL divergence in the trimmed KL divergence matrix
		index.of.min <- which.min(KLdiv.matrix.trimmed)
		arr.ind.of.min <- arrayInd(index.of.min, dim(KLdiv.matrix.trimmed))
		#Save selected KL divergence under the appropriate true species
		KLdiv.all.reps[i, colnames(KLdiv.matrix.trimmed)[arr.ind.of.min[2]]] <- as.numeric(KLdiv.matrix.trimmed[index.of.min])
		#Remove the selected fitted and true species from the KL divergence 
		#matrix by replacing with NA
		KLdiv.matrix.trimmed[arr.ind.of.min[1], ] <- NA
		KLdiv.matrix.trimmed[, arr.ind.of.min[2]] <- NA
		num.pairs.selected <- num.pairs.selected + 1
	}

	#Remove response object
	remove(response)

} #End loop which iterates through all replicates for a treatment level

################################################################################
# 2.3. Store summary of KL divergence results for each true species.

for(sampled.species.num in 1:all.tls$s[tl]) {
    KLdiv.summary.table[tl, paste("KLdiv", sampled.species.num, "min", sep=".")] <- min(KLdiv.all.reps[, sampled.species.num], na.rm=TRUE)
    KLdiv.summary.table[tl, paste("KLdiv", sampled.species.num, "1st.quartile", sep=".")] <- quantile(KLdiv.all.reps[, sampled.species.num], probs=0.25, na.rm=TRUE)
    KLdiv.summary.table[tl, paste("KLdiv", sampled.species.num, "median", sep=".")] <- median(KLdiv.all.reps[, sampled.species.num], na.rm=TRUE)
    KLdiv.summary.table[tl, paste("KLdiv", sampled.species.num, "3rd.quartile", sep=".")] <- quantile(KLdiv.all.reps[, sampled.species.num], probs=0.75, na.rm=TRUE)
    KLdiv.summary.table[tl, paste("KLdiv", sampled.species.num, "max", sep=".")] <- max(KLdiv.all.reps[, sampled.species.num], na.rm=TRUE)
    KLdiv.summary.table[tl, paste("KLdiv", sampled.species.num, "NA.count", sep=".")] <- sum(is.na(KLdiv.all.reps[, sampled.species.num]))
}

#Save running time for this treatment level
all.tls$time.in.minutes[tl] <- difftime(Sys.time(), start.time, units="mins")

} #End loop which iterates through all treatment levels


################################################################################
# 3. Save summary table and running times.
################################################################################

#Set directory to folder to write resulting table
setwd(dir.to.write.table)

#Save summary table
write.table(KLdiv.summary.table, file=paste("KLdivSummaryTable_",
	number.replicates, "reps_",
	format(Sys.time(), "%d%b%Y_%H%M%S"), ".txt", sep=""), sep=",")

#If desired, save running times
if(save.running.times) write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE)

} #Close brackets around code body
