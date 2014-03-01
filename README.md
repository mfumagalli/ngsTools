#ngsTools

NGS (Next-Generation Sequencing) technologies have revolutionized population genetic research by enabling unparalleled data collection from the genomes or subsets of genomes from many individuals. Current technologies produce short fragments of sequenced DNA called _reads_ that are either de novo assembled or mapped to a pre-existing reference genome. This leads to chromosomal positions being sequenced a variable number of times across the genome, usually referred to as the sequencing depth. Individual genotypes are then inferred from the proportion of nucleotide bases covering each site after the reads have been aligned.

Low sequencing depth, along with high error rates stemming from base calling and mapping errors, cause SNP (Single Nucleotide Polymorphism) and genotype calling from NGS data to be associated with considerable statistical uncertainty. Recently, probabilistic models, which take these errors into account, have been proposed to accurately assign genotypes and estimate allele frequencies (e.g. [Nielsen et al., 2012](http://www.ncbi.nlm.nih.gov/pubmed/22911679); for a review [Nielsen et al., 2011](http://www.ncbi.nlm.nih.gov/pubmed/21587300)).

__ngsTools__ is a collection of programs for population genetics analyses from NGS data, taking into account its statistical uncertainty. The methods implemented in these programs do not rely on SNP or genotype calling, and are particularly suitable for low sequencing depth data. An application note illustrating its application has published [Fumagalli et al., 2014](http://www.ncbi.nlm.nih.gov/pubmed/24458950).

NOTE - this repository is intended for general use as it groups together the latest stable version of each tool. Developers may want to check each tool's main repository.

##Packages

*  [ngsSim](https://github.com/mfumagalli/ngsSim) - Simple sequencing read simulator that can generate data for multiple populations with variable levels of depth, error rates, genetic variability, and individual inbreeding ([Kim et al., 2011](http://www.ncbi.nlm.nih.gov/pubmed/21663684)). It generates mapped reads and the corresponding genotype likelihoods, avoiding mapping uncertainty and being extremely fast.

* [ngsF](https://github.com/fgvieira/ngsF) - This program provides a method to estimate individual inbreeding coefficients using an EM algorithm ([Vieira et al., 2013](http://www.ncbi.nlm.nih.gov/pubmed/23950147)). These can provide insight into a population's mating system or demographic history and, more importantly, they can be used as a prior in ANGSD.

* [ngsPopGen](https://github.com/mfumagalli/ngsPopGen) - Several tools to perform population genetic analyses from NGS data ([Fumagalli et al., 2013](http://www.ncbi.nlm.nih.gov/pubmed/23979584), [Fumagalli., 2013](http://www.ncbi.nlm.nih.gov/pubmed/24260275)).
 * __ngsFst__ - Quantifying population genetic differentiation
 * __ngsCovar__ - Population structure via PCA (principal components analysis)
 * __ngs2dSFS__ - Estimate 2D-SFS from posterior probabilities of sample allele frequencies
 * __ngsStat__ - Estimate number of segregating sites, expected average heterozygosity, and number of fixed differences (if data from 2 populations are provided).
* [ngsUtils](https://github.com/mfumagalli/ngsUtils) - General tools to manipulate data.
 * __GetMergedGeno__ - Merge genotype posterior probabilities files
 * __GetSubGeno__ - Select a subset of genotype posterior probabilities files
 * __GetSubSfs__ - Select a subset of sample allele frequency posterior probabilities files
 * __GetSubSim__ - Select a subset of simulated data files
 * __GetSwitchedGeno__ - Switch major/minor in genotype posterior probabilities files
 * __GetSwitchedSfs__ - Switch major/minor in sample allele frequency posterior probabilities files

## Installation

To download ngsTools and its submodules use:

    % git clone --recursive https://github.com/mfumagalli/ngsTools.git

if you prefer, although it is not recomended, you can download it from the Nielsen lab [webpage](http://cteg.berkeley.edu/~nielsen/resources/software/) and run:

    % tar -xvf  ngsTools_20140114.tar.gz

To install these tools just run:

    % cd ngsTools
    % make
    % make test

NOTE: Test scripts do not work for Mac, yet.

Executables are built into each tool directory in the repository. If you wish to clean all binaries and intermediate files:

    % make clean

To get the latest version of the package:

    % git pull
    % git submodule update
    
NOTE for developers: if you wish to make changes and update the whole package:

    % # in the modified repo
    % git commit -a -m 'Local changes...'
    % git push
    % # in ngsTools main repo
    % git commit -a -m 'Submodules updated'
    % git push

## Input Files

All programs receive as input files produced by ANGSD. In general, these files can contain genotype likelihoods, genotype posterior probabilities, sample allele frequency posterior probabilities or an estimate of the SFS (Site Frequency Spectrum). Please refer to each tool's repository for more explanations and examples on how these tools work.
