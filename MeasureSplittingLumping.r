################################################################################
################################################################################
# PURPOSE OF THIS SCRIPT
# This script computes the splitting and lumping metrics for each species group
# and each species within the group. Specifically, for each treatment level,
# it summarizes the distribution across all replicates of 
# E[V(estimated | true)]', E[V(true | estimated)]_i', E[V(true | estimated)]',
# and E[V(true | estimated)]_i'. The results for each treatment level are stored
# in one row of the resulting table.
#
# DATA FILES REQUIRED TO RUN THIS SCRIPT
# Files describing the NMM fitted to the samples from the simulated species
# groups must be available in a directory having no other files.
#
# CONTENTS OF THIS SCRIPT
# 1. Preliminaries
# 2. Compute and save splitting and lumping metrics.
# 3. Save summary table and running times.
################################################################################
################################################################################

{ #Use brackets around code body so that stop() stops script execution immediately

################################################################################
# 1. Preliminaries
################################################################################

################################################################################
# 1.1. Load packages and define working directories

#Set directory to read NMMs
dir.root.NMMs <- "replace with path to directory of all NMM objects"

#Set directories to write results
dir.to.write.table <- "replace with path to directory to write splitting and lumping summary tables"

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
# 1.3. Create file to update running times for calculations at each treatment
# level

save.running.times <- TRUE

if(save.running.times)
{
	setwd(dir.to.write.table)
	running.times.file.name <- paste0("SplittingLumpingMetricRunningTimes_", 
		number.replicates, "replicates.per.level_", 
		format(Sys.time(), "%d%b%Y_%H%M%S"), "_ran_", 
		round(runif(1, min=1, max=1000000)), ".csv")
	write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE)
}

################################################################################
# 2. Compute and save splitting and lumping metrics.
################################################################################

#Create data frames to store results
E.V.true.est.summary <- all.tls
E.V.est.true.summary <- all.tls

#Begin loop which iterates through each treatment level
for(tl in 1:num.tls)
{
    #Store start time for obtaining responses from this treatment level
    start.time <- Sys.time()

    #Define directory from which NMM objects should be read
    dir.tl.NMMs <- paste0(dir.root.NMMs, "/", "NMMs_s", all.tls$s[tl], "/", "NMMs_", all.tls$tl.id[tl])

    s <- all.tls$s[tl]
    n <- all.tls$n[tl]

    ################################################################################
    # 2.2. Obtain or calculate all responses for each replicate at the given
    # treatment level.

    #Initialize matrices to hold lumping and splitting measures, respectively, for
    #each replicate
    E.V.true.est.data <- matrix(data = NA, nrow = number.replicates, ncol = s + 1)
    E.V.est.true.data <- matrix(data = NA, nrow = number.replicates, ncol = s + 1)

    #Begin loop which iterates through all replicates for this treatment level
    for(i in 1:number.replicates)
    {
        #Set working directory to read fitted NMMs
        setwd(dir.tl.NMMs)

        #Load fitted NMM
        NMM <- eval(as.name(load(dir()[i])))

        #Obtain number of species per classification
        true.class <- rep(1:s, NMM$sample.sizes)
        est.class <- NMM$NMM.fit$classification
        n_ij <- table(true.class, est.class)
        n_i = rowSums(n_ij)
        n_j = colSums(n_ij)
        n_ij_sq = n_ij^2

        #Compute lumping metrics for each species and the group overall
        E.V.true.est_i <- n_i - rowSums(sweep(n_ij_sq, 2, n_j, FUN = '/'))
        diff.min.max.E.V.true.est_i <- n_i - (n_i^2 / n)
        E.V.true.est_i.prime <- E.V.true.est_i / diff.min.max.E.V.true.est_i
        E.V.true.est.prime <- sum(E.V.true.est_i) / sum(diff.min.max.E.V.true.est_i)

        E.V.true.est.data[i, 1:s] <- E.V.true.est_i.prime
        E.V.true.est.data[i, s + 1] <- E.V.true.est.prime

        #Compute splitting metrics for each species and the group overall
        E.V.est.true_i <- n_i - rowSums(n_ij_sq) / n_i
        diff.min.max.E.V.est.true_i <- n_i - 1
        E.V.est.true.prime <- sum(E.V.est.true_i) / sum(diff.min.max.E.V.est.true_i)
        diff.min.max.E.V.est.true_i[diff.min.max.E.V.est.true_i < 1] <- NA
        E.V.est.true_i.prime <- E.V.est.true_i / diff.min.max.E.V.est.true_i

        E.V.est.true.data[i, 1:s] <- E.V.est.true_i.prime
        E.V.est.true.data[i, s + 1] <- E.V.est.true.prime

        rm(list = "NMM")
    }

    ################################################################################
    # 2.3. Store summary of results distributions in table.

    #Store summary of E[V(estimated | true)]
    E.V.est.true.summary[tl, paste("group.min")] <- 
        min(E.V.est.true.data[, s + 1], na.rm = TRUE)
    E.V.est.true.summary[tl, paste("group.1st.quartile")] <- 
        quantile(E.V.est.true.data[, s + 1], probs = 0.25, na.rm = TRUE)
    E.V.est.true.summary[tl, paste("group.median")] <-
        median(E.V.est.true.data[, s + 1], na.rm = TRUE)
    E.V.est.true.summary[tl, paste("group.3rd.quartile")] <- 
        quantile(E.V.est.true.data[, s + 1], probs = 0.75, na.rm = TRUE)
    E.V.est.true.summary[tl, paste("group.max")] <- 
        max(E.V.est.true.data[, s + 1], na.rm = TRUE)
    E.V.est.true.summary[tl, paste("group.NA.count")] <- 
        sum(is.na(E.V.est.true.data[, s + 1]))

    #Store summary of E[V(true | estimated)]
    E.V.true.est.summary[tl, paste("group.min")] <- 
        min(E.V.true.est.data[, s + 1], na.rm = TRUE)
    E.V.true.est.summary[tl, paste("group.1st.quartile")] <- 
        quantile(E.V.true.est.data[, s + 1], probs = 0.25, na.rm = TRUE)
    E.V.true.est.summary[tl, paste("group.median")] <-
        median(E.V.true.est.data[, s + 1], na.rm = TRUE)
    E.V.true.est.summary[tl, paste("group.3rd.quartile")] <- 
        quantile(E.V.true.est.data[, s + 1], probs = 0.75, na.rm = TRUE)
    E.V.true.est.summary[tl, paste("group.max")] <- 
        max(E.V.true.est.data[, s + 1], na.rm = TRUE)
    E.V.true.est.summary[tl, paste("group.NA.count")] <- 
        sum(is.na(E.V.true.est.data[, s + 1]))

    for(i in 1:s) {
        #Store summary of results for splitting contribution of each species
        E.V.est.true.summary[tl, paste("sp", i, "min", sep = ".")] <- 
            min(E.V.est.true.data[, i], na.rm = TRUE)
        E.V.est.true.summary[tl, paste("sp", i, "1st.quartile", sep = ".")] <- 
            quantile(E.V.est.true.data[, i], probs = 0.25, na.rm = TRUE)
        E.V.est.true.summary[tl, paste("sp", i, "median", sep = ".")] <-
            median(E.V.est.true.data[, i], na.rm = TRUE)
        E.V.est.true.summary[tl, paste("sp", i, "3rd.quartile", sep = ".")] <- 
            quantile(E.V.est.true.data[, i], probs = 0.75, na.rm = TRUE)
        E.V.est.true.summary[tl, paste("sp", i, "max", sep = ".")] <- 
            max(E.V.est.true.data[, i], na.rm = TRUE)
        E.V.est.true.summary[tl, paste("sp" ,i, "NA.count", sep = ".")] <- 
            sum(is.na(E.V.est.true.data[, i]))

        #Store summary of results for lumping contribution of each species
        E.V.true.est.summary[tl, paste("sp", i, "min", sep = ".")] <- 
            min(E.V.true.est.data[, i], na.rm = TRUE)
        E.V.true.est.summary[tl, paste("sp", i, "1st.quartile", sep = ".")] <- 
            quantile(E.V.true.est.data[, i], probs = 0.25, na.rm = TRUE)
        E.V.true.est.summary[tl, paste("sp", i, "median", sep = ".")] <-
            median(E.V.true.est.data[, i], na.rm = TRUE)
        E.V.true.est.summary[tl, paste("sp", i, "3rd.quartile", sep = ".")] <- 
            quantile(E.V.true.est.data[, i], probs = 0.75, na.rm = TRUE)
        E.V.true.est.summary[tl, paste("sp", i, "max", sep = ".")] <- 
            max(E.V.true.est.data[, i], na.rm = TRUE)
        E.V.true.est.summary[tl, paste("sp" ,i, "NA.count", sep = ".")] <- 
            sum(is.na(E.V.true.est.data[, i]))
    }
    
} #End of for loop which iterated through each treatment level

################################################################################
# 3. Save summary table and running times.
################################################################################

setwd(paste0(dir.to.write.table, "/splitting")) 

write.table(E.V.est.true.summary, file=paste("VarEstTrueSummary_",
    number.replicates, "reps_",
    format(Sys.time(), "%d%b%Y_%H%M%S"), ".txt", sep=""), sep=",")

setwd(paste0(dir.to.write.table, "/lumping"))

write.table(E.V.true.est.summary, file=paste("VarTrueEstSummary_",
    number.replicates, "reps_",
    format(Sys.time(), "%d%b%Y_%H%M%S"), ".txt", sep=""), sep=",")

if(save.running.times) {
    setwd(dir.to.write.table)
    write.table(all.tls, file=running.times.file.name, sep=",", row.names=FALSE)
}

} #Close brackets around code body
