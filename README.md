# Species Delimitation Using Phenotypic Data and Normal Mixture Models: How Many Specimens Should Be Measured?

This repository holds the final implementation of the simulation study in Species Delimitation Using Phenotypic Data and Normal Mixture Models: How Many Specimens Should Be Measured?

## Usage Notes for Replicating the Simulation

This simulation study was carried out using R-3.5.2.

At the top of each script is a description of the script's purpose, any data files needed in order to run the script, and an outline of the script's contents.

Each script defines its parameters in a section labeled "Preliminaries." Prior to running a script, these parameters should be checked and adjusted as necessary to match the parameters of the desired simulation.

### Sequence of Execution
The following orders should be preserved, to ensure that the data files needed in order to run a script are generated before attempting to run that script.

For example, SampleSpeciesGroups.R should run only after SimulateSpecies_p_Dimensions.R has run, as SampleSpeciesGroups.R relies on data files produced by SimulateSpeciesGroups_p_Dimensions.R; similarly, GetResponsesFromNMM.R should run only after SampleSpeciesGroups.R has run; there is no order between SampleSpeciesGroups.R and KLDivergence_UpperBounds.R; etc.

<ol>
<li>SimulateSpeciesGroups_p_Dimensions.R</li>
<li>SampleSpeciesGroups.R</li>
<li>GetResponsesFromNMM.R</li>
<li>GetGreedyKLDivergence.R</li>
</ol>

<ol>
<li>KLDivergence_UpperBounds.R</li>
<li>GetGreedyKLDivergence.R</li>
</ol>

<ol>
<li>SampleSpeciesGroups.R</li>
<li>MeasureSplittingLumping.R</li>
</ol>

## Additional Material
#### GenerateEmpiricalDistributionForSampling.R
GenerateEmpiricalDistributionForSampling.R was used to determine what sampling distributions should be used in the simulation study. The results from this script impacted the implementation coded in SampleSpeciesGroups.R. GenerateEmpircalDistributionForSampling.R is included in this repository to demonstrate how the sampling distributions were chosen; it is not necessary to run GenerateEmpiricalDistributionForSampling.R in order to replicate the simulation study, unless the replicator intends to change the sampling distributions -- if this is the case, looking at GenerateEmpiricalDistributionForSampling.R may help in choosing a new sampling distribution, but the replicator should (also) modify SampleSpeciesGroups.R accordingly.