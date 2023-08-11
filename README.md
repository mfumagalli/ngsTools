ngsTools
------------------------------
------------------------------

NGS (Next-Generation Sequencing) technologies have revolutionised population genetic research by enabling unparalleled data collection from the genomes or subsets of genomes from many individuals.
Current technologies produce short fragments of sequenced DNA called _reads_ that are either _de novo_ assembled or mapped to a pre-existing reference genome. 
This leads to chromosomal positions being sequenced a variable number of times across the genome.
This parameter is usually referred to as the sequencing depth. 
Individual genotypes are then inferred from the proportion of nucleotide bases covering each site after the reads have been aligned.

Low sequencing depth and high error rates stemming from base calling and mapping errors can cause SNP (Single Nucleotide Polymorphism) and genotype calling from NGS data to be associated with considerable statistical uncertainty. 
Probabilistic models, which take these errors into account, have been proposed to accurately assign genotypes and estimate allele frequencies (e.g. [Nielsen et al. 2012](http://www.ncbi.nlm.nih.gov/pubmed/22911679); for a review [Nielsen et al. 2011](http://www.ncbi.nlm.nih.gov/pubmed/21587300)).

__ngsTools__ is a collection of programs for population genetics analyses from NGS data, taking into account data statistical uncertainty.
The methods implemented in these programs do not rely on SNP or genotype calling, and are particularly suitable for low sequencing depth data. 
An application note illustrating its application has published ([Fumagalli et al. 2014](http://www.ncbi.nlm.nih.gov/pubmed/24458950)).

__NOTE__: this repository is intended for general use as it groups together the latest stable version of each tool. 
Developers (and only them) may want to check each tool's main repository.


Packages
-----------------------------

* [ANGSD](https://github.com/ANGSD/angsd) - Software for analyzing next generation sequencing data taking genotype uncertainty into account by working with genotype probabilities (instead of called genotypes). This is especially useful for low and medium depth data ([Korneliussen et al., 2014](https://www.ncbi.nlm.nih.gov/pubmed/25420514)). **NOTE**: this program is **NOT** developed by `ngsTools` so, if you have any questions about it or encounter any errors/bugs, please check its [wiki](http://www.popgen.dk/wiki/index.php/ANGSD) or contact its authors.

* [ngsSim](https://github.com/mfumagalli/ngsSim) - Simple sequencing read simulator that can generate data for multiple populations with variable levels of depth, error rates, genetic variability, and individual inbreeding ([Kim et al., 2011](http://www.ncbi.nlm.nih.gov/pubmed/21663684)). It generates mapped reads and the corresponding genotype likelihoods, avoiding mapping uncertainty and being extremely fast.

* [ngsF](https://github.com/fgvieira/ngsF) - This program provides a method to estimate individual inbreeding coefficients using an EM algorithm ([Vieira et al., 2013](http://www.ncbi.nlm.nih.gov/pubmed/23950147)). These can provide insight into a population's mating system or demographic history and, more importantly, they can be used as a prior in ANGSD.

* [ngsPopGen](https://github.com/mfumagalli/ngsPopGen) - Several tools to perform population genetic analyses from NGS data ([Fumagalli et al., 2013](http://www.ncbi.nlm.nih.gov/pubmed/23979584), [Fumagalli, 2013](http://www.ncbi.nlm.nih.gov/pubmed/24260275)).
	* __ngsFst__ - Quantifying population genetic differentiation
	* __ngsCovar__ - Population structure via PCA (principal components analysis)
	* __ngs2dSFS__ - Estimate 2D-SFS from posterior probabilities of sample allele frequencies
	* __ngsStat__ - Estimate number of segregating sites, expected average heterozygosity and other nucleotide diversity indexes

* [ngsUtils](https://github.com/mfumagalli/ngsUtils) - General tools to manipulate data produced by ngsTools.
	* __GetMergedGeno__ - Merge genotype posterior probabilities files
	* __GetSubGeno__ - Select a subset of genotype posterior probabilities files
	* __GetSubSim__ - Select a subset of simulated data files
	* __GetSwitchedGeno__ - Switch major/minor in genotype posterior probabilities files

* [ngsDist](https://github.com/fgvieira/ngsDist) is a program that estimates genetic distances from genotype posterior probabilities ([Vieira et al., 2016](http://onlinelibrary.wiley.com/doi/10.1111/bij.12511/abstract)).

* [ngsF-HMM](https://github.com/fgvieira/ngsF-HMM) is a program developed and written by F.G. Vieira to estimate per-individual inbreeding tracts using a two-state Hidden Markov Model ([Vieira et al. 2016](http://bioinformatics.oxfordjournals.org/content/32/14/2096)). It uses a probabilistic framework that takes the uncertainty of genotype's assignation into account; making it specially suited for low-quality or low-coverage datasets. It is not officially part of `ngsTools` so it must be installed separately.

* [ngsLD](https://github.com/fgvieira/ngsLD) is a program that calculates pairwise Linkage Disequilibrium (LD) under a probabilistic framework ([Fox et al., 2019](https://doi.org/10.1093/bioinformatics/btz200)).

* [HMMploidy](https://github.com/SamueleSoraggi/HMMploidy) is a program written by S. Soraggi to infer ploidy levels from low-coverage sequencing data. It is not officially part of `ngsTools` so it must be installed separately.



Installation
--------------------------------------------------

`ngsTools` can be easily installed but some packages have some external dependencies:

* Mandatory:
  * `gcc`: >= 4.9.2 tested on Debian 7.8 (wheezy)
  * `zlib`: v1.2.7 tested on Debian 7.8 (wheezy)
  * `gsl` : v1.15 tested on Debian 7.8 (wheezy)
  * `libbz2.so`, required by **htslib**
  * `liblzma.so`, required by **htslib**
  * `libcurl.so`, required by **angsd**
* Optional (only needed for testing or auxilliary scripts):
  * `md5sum`
  * `samtools`
  * `Perl` packages: `Getopt::Long`, `Graph::Easy`, `Math::BigFloat`, and `IO::Zlib`
  * `R` packages: `optparse`, `tools`, `ggplot2`, `reshape2`, `plyr`, `gtools`, `LDheatmap`, `ape`, `grid`, `methods`, `phangorn`, and `plot3D`

If you have issues with `gsl` package, then on linux make sure you have these packages installed: `gsl-bin` `libgsl-dbg` `libgsl-dev` `libgslcblas0`.

To download ngsTools and its submodules use:

    % git clone --recursive https://github.com/mfumagalli/ngsTools.git

If you prefer, although it is not recomended, you can download a zipped folder on the right side of this page ("Download ZIP"). 

To install these tools just run:

    % cd ngsTools
    % make

To run the tests (currently deprecated):

    % make test

Executables are built into each tool directory in the repository. If you wish to clean all binaries and intermediate files:

    % make clean

To get the latest version of ngsTools package:

    % git pull
    % git submodule update

__NOTE__: for developers only: if you wish to make changes and update the whole package:

    # in the modified repo
    # be sure to be on the master branch: git checkout master
    % git commit -a -m 'Local changes...'
    % git push
    # in ngsTools main repo
    % git commit -a -m 'Submodules updated'
    % git push
    # check that everything went well: 
    % git status

To build the tool set into an [Apptainer](https://apptainer.org/) or
[SingularityCE](https://sylabs.io/singularity/) container you can take the
following definition file as a starting point:

<details>
<summary>
Apptainer/SingularityCE container definition file
</summary>

```
Bootstrap: docker
From: ubuntu:20.04

%post
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -y

    apt install -y git build-essential pkg-config
    apt install -y libz-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libgsl-dev

    git clone --recursive https://github.com/mfumagalli/ngsTools.git
    cd ngsTools

    # this will build a specific commit (in this case commit 6505f80) #
    # alternatively you can build the latest commit by commenting out the next two lines
    git checkout 6505f80
    git submodule update --init --recursive

    make

%environment
    export LC_ALL=C

%runscript
    export PATH=/ngsTools/angsd:$PATH
    export PATH=/ngsTools/ngsDist:$PATH
    export PATH=/ngsTools/ngsF:$PATH
    export PATH=/ngsTools/ngsF-HMM:$PATH
    export PATH=/ngsTools/ngsLD:$PATH
    export PATH=/ngsTools/ngsPopGen:$PATH
    export PATH=/ngsTools/ngsSim:$PATH
    export PATH=/ngsTools/ngsUtils:$PATH

    $@

%labels
    Author Radovan Bast

%help
    You can build this container with:
    $ sudo singularity build container.sif container.def
    This is how I use this container image:
    $ ./container.sif ngsDist [other arguments]
    $ ./container.sif [some other NGS tool]
```
</details>


Input Files
-----------------------------------

All programs receive as input files produced by ANGSD. 
In general, these files can contain genotype likelihoods, genotype posterior probabilities, sample allele frequency posterior probabilities or an estimate of the SFS (Site Frequency Spectrum). 
Please refer to each tool's repository or the comprehensive [Tutorial](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md) for more explanations and examples on how these tools work.


Tutorial
-------------------------------------------

A tutorial on some analyses using ngsTools/ANGSD from BAM files can be found [here](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).
In this Tutorial, you will find how to filter your data, assess population structure and estimate summary statistics using these tools for low-depth data.
For most cases, you will find all the information you need here.


Additional information
--------------------------------------------

__Authors__

Matteo Fumagalli & Filipe G. Vieira.
Other programmers and developers: Tyler Lynderoth, Rasmus Nielsen.
Some lines of code have been and adapted from: Thorfinn Korneliussen, Anders Albrechtsen, Jacob Crawford, Dean Ousby, Martin Sykora, Leo Diaz.

__Updates__

If you want to be updated about new releases and fixed bugs please follow the github repository.
For specific queries on the code, please use github features to raise issues.
For informal questions (better on ngsTools and not ANGSD) feel free to contact [Matteo Fumagalli](https://www.qmul.ac.uk/sbbs/staff/matteo-fumagalli.html).

__Citation__

ngsTools package can be cited as:

	ngsTools: methods for population genetics analyses from next-generation sequencing data.
	Fumagalli M, Vieira FG, Linderoth T, Nielsen R.
	Bioinformatics. 2014 May 15;30(10):1486-7

ANGSD can be cited as:

	ANGSD: Analysis of Next Generation Sequencing Data.
	Korneliussen T, Albrechtsen A, Nielsen R
	BMC Bioinformatics. 2014 Nov 25;15(1):356

	SNP calling, genotype calling, and sample allele frequency estimation from New-Generation Sequencing data.
	Nielsen R, Korneliussen T, Albrechtsen A, Li Y, Wang J
	PLoS One. 2012;7(7):e37558

FST and PCA methods can be cited as:

	Quantifying Population Genetic Differentiation from Next-Generation Sequencing Data.
	Fumagalli M, Vieira FG, Korneliussen TS, Linderoth T, Huerta-Sánchez E, Albrechtsen A, Nielsen R
	Genetics. 2013 Nov;195(3):979-92

Inbreeding estimation can be cited as:

	Estimating inbreeding coefficients from NGS data: impact on genotype calling and allele frequency estimation.
	Vieira FG, Fumagalli M, Albrechtsen A, Nielsen R
	Genome Res. 2013 Nov;23(11):1852-61

	Estimating IBD tracts from low coverage NGS data
	Vieira FG, Albrechtsen A and Nielsen R
	Bioinformatics. 2016; 32:2096-2102

Nucleotide diversity estimates from NGS data implemented here have been proposed in:

	Sequencing of 50 human exomes reveals adaptation to high altitude.
	Yi X, Liang Y, Huerta-Sanchez E, Jin X, Cuo ZX, Pool JE, Xu X, Jiang H, Vinckenbosch N, Korneliussen TS, Zheng H, Liu T, He W, Li K, Luo R, Nie X, Wu H, Zhao M, Cao H, Zou J, Shan Y, Li S, Yang Q, Asan, Ni P, Tian G, Xu J, Liu X, Jiang T, Wu R, Zhou G, Tang M, Qin J, Wang T, Feng S, Li G, Huasang, Luosang J, Wang W, Chen F, Wang Y, Zheng X, Li Z, Bianba Z, Yang G, Wang X, Tang S, Gao G, Chen Y, Luo Z, Gusang L, Cao Z, Zhang Q, Ouyang W, Ren X, Liang H, Zheng H, Huang Y, Li J, Bolund L, Kristiansen K, Li Y, Zhang Y, Zhang X, Li R, Li S, Yang H, Nielsen R, Wang J, Wang J
	Science. 2010 Jul 2;329(5987):75-8

	Calculation of Tajima's D and other neutrality test statistics from low depth next-generation sequencing data.
	Korneliussen TS, Moltke I, Albrechtsen A, Nielsen R
	BMC Bioinformatics. 2013 Oct 2;14(1):289

	Assessing the effect of sequencing depth and sample size in population genetics inferences.
	Fumagalli M
	PLoS One. 2013 Nov 18;8(11):e79667

Estimation of genetic distances have been described here:

	Improving the estimation of genetic distances from Next-Generation Sequencing data.
	Vieira FG, Lassalle F, Korneliussen TS, Fumagalli M
	Biological Journal of the Linnean Society. Special Issue: Collections-Based Research in the Genomic Era. 117(1):139–149

LD estimation has been published here:

	Fox EA, Wright AE, Fumagalli M, Vieira FG. ngsLD: evaluating linkage disequilibrium using genotype likelihoods. Bioinformatics. 2019 Oct 1;35(19):3855-3856. doi: 10.1093/bioinformatics/btz200. PMID: 30903149.
	
HMMploidy can be cited as:

	Soraggi, Samuele; Rhodes, Johanna; Altinkaya, Isin; Tarrant, Oliver; Balloux, Francois; Fisher, Matthew C; Fumagalli, Matteo. HMMploidy: inference of ploidy levels from short-read sequencing data. Peer Community Journal, Volume 2 (2022), article no. e60. doi : 10.24072/pcjournal.178.



