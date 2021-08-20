################################################################################
################################################################################
#
# PURPOSE OF THIS SCRIPT
# This script demonstrates how the empirical distributions of specimens among
# species were selected. Specifically, it may be run for species groups of 4 or
# 8 species, where the desired proportion of specimens sampled from the more
# abundant half of the species is 0.5, 0.7, or 0.9. The empircal distributions
# when the mixing proportion is slightly uneven or uneven (i.e. 0.7 or 0.9 of
# the specimens are from the more abundant half of the species) are based 
# roughly on the log-series distribution. When the mixing proportion is even,
# the empirical distribution attempts to sample specimens from species as
# equally as possible.
#
# DATA FILES REQUIRED TO RUN THIS SCRIPT
# This script may run on its own, without other data files.
#
# CONTENTS OF THIS SCRIPT
# 1. Preliminaries: Define the variables that determine the distribution.
# 2. Calculate the proportion of specimens per species for the empirical
#    distribution.
# 3. Plot the empirical distribution and ensure that it is valid for the given
#    constraints.
#
################################################################################
################################################################################

################################################################################
# 1. Preliminaries: Define the variables that determine the distribution.
################################################################################

# Define the sample size which should be used to generate the empirical
# distribution. This should be the smallest sample size that will be used in the
# simulations.
n <- 50

# Define the number of species in each species group.
# Select one of the following:
s <- 4
#s <- 8

# Define the proportion of specimens from the more abundant half of the 
# species. (i.e. if s = 8, then the 4 most abundant species will have a total
# of half.of.species.proportion*sample.size specimens.)
# Select one of the following:
#half.of.species.proportion <- 0.5
half.of.species.proportion <- 0.7
#half.of.species.proportion <- 0.9

################################################################################
# 2. Calculate the proportion of specimens per species for the empirical
#    distribution.
################################################################################

# Arbitrary species IDs
X <- 1:s

if(half.of.species.proportion == 0.5){
	# Allocate specimens to species equally
	specimens.per.species <- rep(n%/%s, times=s)
	remainder <- n%%s
	if(remainder > 0)
		specimens.per.species[1:remainder] <- specimens.per.species[1:remainder] + 1
	proportion.of.specimens.per.species <- specimens.per.species/n

} else { 
	# Allocate specimens to species roughly based on the log-series distribution
	
	# Set omega based on the proportion set above and the number of species per
	# species group. (We found these omega values via guess-and-check, so they are
	# hard-coded below.)
	if(s == 8){
		if(half.of.species.proportion == 0.9) omega <- 0.7539
		if(half.of.species.proportion == 0.7) omega <- 0.9226
	}
	if(s == 4){
		if(half.of.species.proportion == 0.9) omega <- 0.5038
		if(half.of.species.proportion == 0.7) omega <- 0.7960
	}

	# Approximate fractions for each species using log-series distribution
	fractions.per.species <- (-1/log(1-omega))*(omega^X)/X

	# Plot fractions.per.species vs. arbitrary species IDs
	plot(X, fractions.per.species, ylim = c(0,1), xaxt="n",
		main = paste("Log Series Distribution for omega = ", omega, 
		"\n corresponding to the proportion of half of the species = ", 
		half.of.species.proportion), 
		ylab = "Fraction of specimens", xlab = "Species ID",
		type = "o", pch = 19)
	axis(1, at=1:s)

	# Calculate number of specimens per species (prior to allocating the remainder)
	specimens.per.species.initial <- round(n*fractions.per.species)
	# Ensure that the number of specimens belonging to the more common the half of
	# the species is exactly half.of.species.proportion*n
	if(sum(specimens.per.species.initial[1:(s/2)]) != half.of.species.proportion*n){
		stop("The number of specimens allocated to the more common half of the
			species in not exactly half.of.species.proportion*n.
			Suggested fix: deallocate specimens from specimens.per.species.initial[s/2].")
	}
	specimens.per.species.initial

	# Determine how many specimens remain to be allocated
	remaining.specimens <- n - sum(specimens.per.species.initial)
	remaining.specimens

	# Allocate the remaining specimens equally among the less common half of the
	# species.
	equal.allocation.remaining.specimens = floor(remaining.specimens/(s/2))
	add.specimens = rep(c(0, equal.allocation.remaining.specimens), each=s/2)
	remaining.specimens <- remaining.specimens - sum(add.specimens)
	# If equal allocation was not possible, allocate one specimen to each of the
	# last remaining.specimens species. 
	# (Note, remaining.specimens < s/2 at this point.)
	if(remaining.specimens > 0)
		add.specimens[(s - remaining.specimens + 1):s] <- add.specimens[(s - remaining.specimens + 1):s] + 1

	# Calculate the final number of specimens allocated to each species
	specimens.per.species <- specimens.per.species.initial + add.specimens
	specimens.per.species
	# Calculate the proportion of specimens per species
	# This series of proportions defines this empirical distribution
	proportion.of.specimens.per.species <- specimens.per.species/n
	proportion.of.specimens.per.species
}


################################################################################
# 3. Plot the empirical distribution and ensure that it is valid for the given
#    constraints.
################################################################################

plot(X, proportion.of.specimens.per.species, ylim = c(0,1), xaxt="n",
	main = paste("Distribution of specimens \n for the proportion of half of the species =", 
	half.of.species.proportion, "and n =", n),
	ylab = "Proportion of specimens", xlab = "Species ID", 
	type = "o", pch = 19)
axis(1, at=1:s)

# Print warnings if the generated distribution does not use the variables
# values defined at the beginning of the script
if(sum(specimens.per.species) != n){
	warning(paste(sum(specimens.per.species), 
		"specimens were used to create the empirical distribution, but", 
		n, "specimens was the expected number."))
}
if(sum(specimens.per.species[1:(s/2)])/sum(specimens.per.species) != half.of.species.proportion){
	warning(paste("This empirical distribution allocates", 
		sum(specimens.per.species[1:(s/2)])/sum(specimens.per.species), 
		"specimens to the more common half of the species, but the desired proportion was", 
		half.of.species.proportion, "."))
}
