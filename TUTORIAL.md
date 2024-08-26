
A tutorial for some basic analyses using ngsTools/ANGSD
===============

Installation
--------------------

We recommend to install ANGSD separately following the instructions [here](http://popgen.dk/angsd/index.php/Download_and_installation):

    git clone https://github.com/samtools/htslib.git
    git clone https://github.com/ANGSD/angsd.git 
    cd htslib;make;cd ../angsd ;make HTSSRC=../htslib;cd ..

Please be sure you are using the most updated version of ngsTools. In doubt please run: 

    git checkout master
    git pull
    git submodule update

and in case you need to update something, type:

    make clean
    make
	
inside the ngsTools directory.

Settings
----------

In this tutorial we will be using several programs including ngsTools, ANGSD and NGSadmix to perform population genetics analyses from low-depth sequencing data.
Please note that [ANGSD](http://popgen.dk/angsd/index.php/Main_Page#Overview) and [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) have not been developed by us and therefore questions on these tools should be addressed to their Authors.
However, given the utility of such tools, we felt the need to include them to present a more comprehensive view on the application of this probabilistic approach to process NGS data in population genetics.
Finally, we are using [SAMtools](http://samtools.sourceforge.net/) for indexing files, [FastMe](http://www.atgc-montpellier.fr/fastme/) for plotting trees and [R](https://www.r-project.org/) for manipulating and plotting results. 
This tutorial has been tested on a Linux Ubuntu machine with ANGSD version 0.917-126-gb1d9615 (htslib: 1.4-30-g6a50863), SAMtools version 1.4.1, FastME version 2.1.5, R version 3.4.0 (with packages: methods, optparse, ggplot2, ape, phangorn, plot3D, plot3Drgl, rgl).

Please note that R scripts provided for plot here are for illustrative purposes with the example data sets only.
They should be modified if you want to use them with your data set.

First, set directories to all required programs depending on where you installed them, for instance these are my paths:

    NGSTOOLS=~/Software/ngsTools
    ANGSD=~/Software/angsd
    NGSADMIX=~/Software/NGSadmix/NGSadmix

    SAMTOOLS=~/Software/samtools-1.4.1/samtools
    HTSLIB=~/Software/htslib-1.4.1
    FASTME=~/Software/fastme-2.1.5/src/fastme


Second, create all directories where you will be working:

    mkdir Tutorial
    cd Tutorial
    mkdir Data
    mkdir Results

Our goal in this tutorial is to show how to go from BAM files to summary statistics using ngsTools/ANGSD.
Namely we will cover the following topics:
* basic data filtering
* population structure
* inbreeding
* site frequency spectrum
* population genetic differentiation
* nucleotide diversity

Data
----------

As an illustration, we will use 30 BAM files of human samples (of African, European, and Native American descent), a reference genome, and a putative ancestral sequence.
BAM files have a mean depth of around 8X.
The human data represents a small genomic region (1Mbp on chromosome 11) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).
All data is publicly available.

A pipeline to retrieve such data is provided [here](https://github.com/mfumagalli/ngsTools/blob/master/Scripts/data.sh).
You need to have 'gunzip' and 'wget' installed in your /usr/bin to run this.
You also need to specify to paths to 'samtools' (tested with version 1.4.1) and 'bgzip' (tested with htslib version 1.4.1).
```
    cp $NGSTOOLS/Files/*.txt .
    bash $NGSTOOLS/Scripts/data.sh $SAMTOOLS $HTSLIB/bgzip 
```

Now we have 30 BAM files at low/medium depth, a reference and an ancestral sequence in FASTA format.
```
	cat Data/download.log
```
and the list with BAM files has been written to 'ALL.bamlist'.
```
	cat ALL.bamlist
```
Check that all your BAM files in 'Data/download.log' have been downloaded correctly (file size different than 0).

As a note for the general use, in case an ancestral sequence is not available, analyses on FST, PCA, nucleotide diversity (but not the number of fixed differences) can be carried out using the reference sequence to polarise your data. 
Please be aware that, under this scenario, some quantities (e.g. the unfolded joint site frequency spectrum) will be nonsense.
We also discourage to use the folded site frequency spectrum option to estimate FST/PBS.


Basic filtering using ANGSD
----------------------

Here we will use ANGSD to analyse our data, for filtering and for generating files as input for ngsTools.
To see a full list of options in ANGSD type:
```
$ANGSD/angsd
```
and you should see something like:
```
	-> angsd version: 0.917-126-gb1d9615 (htslib: 1.4-30-g6a50863) build(Jun  7 2017 14:32:36)

	-> angsd version: 0.917-126-gb1d9615 (htslib: 1.4-30-g6a50863) build(Jun  7 2017 14:32:36)
	-> Please use the website "http://www.popgen.dk/angsd" as reference
	-> Use -nThreads or -P for number of threads allocated to the program
Overview of methods:
	-GL		Estimate genotype likelihoods
	-doCounts	Calculate various counts statistics
	-doAsso		Perform association study
	-doMaf		Estimate allele frequencies
	-doError	Estimate the type specific error rates
	-doAncError	Estimate the errorrate based on perfect fastas
	-HWE_pval		Est inbreedning per site or use as filter
	-doGeno		Call genotypes
	-doFasta	Generate a fasta for a BAM file
	-doAbbababa	Perform an ABBA-BABA test
	-sites		Analyse specific sites (can force major/minor)
	-doSaf		Estimate the SFS and/or neutrality tests genotype calling
	-doHetPlas	Estimate hetplasmy by calculating a pooled haploid frequency

	Below are options that can be usefull
	-bam		Options relating to bam reading
	-doMajorMinor	Infer the major/minor using different approaches
	-ref/-anc	Read reference or ancestral genome
	-doSNPstat	Calculate various SNPstat
	-cigstat	Printout CIGAR stat across readlength
	many others

For information of specific options type: 
	./angsd METHODNAME eg 
		./angsd -GL
		./angsd -doMaf
		./angsd -doAsso etc
		./angsd sites for information about indexing -sites files
Examples:
	Estimate MAF for bam files in 'list'
		'./angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1'

```

ANGSD can accept several input files, as described [here](http://popgen.dk/angsd/index.php/Input):

* BAM/CRAM
* Pileup
* Genotype likelihood/probability files
* VCF

We can do some basic filtering of our data directly with ANGSD.
These filters are based on:

* quality and depth, see [here](http://www.popgen.dk/angsd/index.php/Filters)
* SNP quality, see [here](http://popgen.dk/angsd/index.php/SnpFilters)
* sites, see [here](http://popgen.dk/angsd/index.php/Sites)

If the input file is in BAM format, the possible options are:
```
$ANGSD/angsd -bam
...
---------------
parseArgs_bambi.cpp: bam reader:
	-bam/-b		(null)	(list of BAM/CRAM files)
	-i		(null)	(Single BAM/CRAM file)
	-r		(null)	Supply a single region in commandline (see examples below)
	-rf		(null)	Supply multiple regions in a file (see examples below)
	-remove_bads	1	Discard 'bad' reads, (flag >=256) 
	-uniqueOnly	0	Discards reads that doesn't map uniquely
	-show		0	Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
	-minMapQ	0	Discard reads with mapping quality below
	-minQ		13	Discard bases with base quality below
	-trim		0	Number of based to discard at both ends of the reads
	-trim		0	Number of based to discard at 5' ends of the reads
	-trim		0	Number of based to discard at 3' ends of the reads
	-only_proper_pairs 1	Only use reads where the mate could be mapped
	-C		0	adjust mapQ for excessive mismatches (as SAMtools), supply -ref
	-baq		0	adjust qscores around indels (as SAMtools), supply -ref
	-checkBamHeaders 1	Exit if difference in BAM headers
	-doCheck	1	Keep going even if datafile is not suffixed with .bam/.cram
	-downSample	0.000000	Downsample to the fraction of original data
	-nReads		50	Number of reads to pop from each BAM/CRAMs
	-minChunkSize	250	Minimum size of chunk sent to analyses

Examples for region specification:
		chr:		Use entire chromosome: chr
		chr:start-	Use region from start to end of chr
		chr:-stop	Use region from beginning of chromosome: chr to stop
		chr:start-stop	Use region from start to stop from chromosome: chr
		chr:site	Use single site on chromosome: chr
```

Some basic filtering consists in removing, for instance, reads with low quality and/or with multiple hits, and this can be achieved using the parameters ```-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1```.

Also, you may want to remove reads with low mapping quality and sites with low quality or covered by few reads (very low depth).
However, it is necessary to know the overall distribution of per-site depth, in order to avoid filtering too many sites.
We first derive the distribution of quality scores and depth on our data set using ```-doQsDist 1 -doDepth 1 -doCounts 1```.

Before proceeding, we need to specify the reference and ancestral sequences.
```
	REF=Data/ref.fa.gz
	ANC=Data/anc.fa.gz
```

We are now printing the distribution of quality scores and per-site depths (global and per-sample).
```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL.qc -r 11\
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
	-minMapQ 20 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500 &> /dev/null
```
As input we give the list of BAM files with option `-b` and then specify the references sequence with `-ref` and the prefix for output files with `-out`.
Please note that we are analysing only part of chromosome 11, indicated by `-r 11`.
Additionally, ```-baq 1``` computes base alignment quality as explained here ([BAQ](http://samtools.sourceforge.net/mpileup.shtml)) to rule out false SNPs close to INDELS, and ```-trim 0``` means that we are not trimming the ends of reads.
With ```-minMapQ 20``` we filter out reads with low mapping quality.
Finally, ```-maxDepth 500``` means that all sites with depth equal or greater than this value will be binned together, and ```-P 4``` means that I am using 4 threads.

You can have a look at the files generated:
```
ls Results/*
...
Results/ALL.qc.arg  Results/ALL.qc.depthGlobal  Results/ALL.qc.depthSample  Results/ALL.qc.qs
...
```
and open them:
```
# counts of quality scores
less -S Results/ALL.qc.qs
# counts of per-sample depth
less -S Results/ALL.qc.depthSample 
wc -l Results/ALL.qc.depthSample # 30 Results/ALL.qc.depthSample
# counts of global depth
less -S Results/ALL.qc.depthGlobal 
```

It is convenient to compute the percentiles of these distributions (and visualize them) in order to make an informative decision on the threshold values we will use for our filtering.
We provide you with a quick R script to pursue this.
```
Rscript $NGSTOOLS/Scripts/plotQC.R Results/ALL.qc 2> /dev/null
```
Have a look at the output files:
```
less -S Results/ALL.qc.info
evince Results/ALL.qc.pdf
```

We may also want to remove sites where half of the individual have no data. This is achieved by the `-minInd` option.
After inspecting these output files, a possible choice of parameters may be:

Parameter | Meaning |
--- | --- |
-minMapQ 20 | minimum mapping quality of 20 |
-minQ 20 | minimum base quality of 20 |
-minInd 15 | use only sites with data from at least 15 individuals |
-setMinDepth 60 | minimum total depth |
-setMaxDepth 400 | maximum total depth |

Please note that ANGSD can also compute more sophisticated metrics to filter out SNPs, as described [here](http://popgen.dk/angsd/index.php/SnpFilters), mostly based on:

* strand bias
* deviation from HWE
* quality score bias

Moreover, additional filtering should be considered.
For instance transitions (A<->G, C<->T) are more likely than transversions, so we expect the ts/tv ratio to be greater than 0.5.
However, we are not discussing this additional filtering options in this tutorial.


Population structure
---------------------------------------

Suppore we want to investigate the population structure our samples: PEL (Peruvians), TSI (Europeans), LWK (Africans).
One solution would be to perform a Principal Component Analysis (PCA) or a Multidimensional Scaling (MDS) or some clustering based on genetic distances among samples.
We are here showing how to perform these analyses using ngsTools/ANGSD in case of low-depth data.

To do this, we first need to assign genotype probabilities at each site for each individual.
The specific option in ANGSD is `-doGeno`.

```
$ANGSD/angsd -doGeno
...
-doGeno	0
	1: write major and minor
	2: write the called genotype encoded as -1,0,1,2, -1=not called
	4: write the called genotype directly: eg AA,AC etc 
	8: write the posterior probability of all possible genotypes
	16: write the posterior probability of called genotype
	32: write the posterior probabilities of the 3 gentypes as binary
	-> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
	-postCutoff=0.333333 (Only genotype to missing if below this threshold)
	-geno_minDepth=-1	(-1 indicates no cutof)
	-geno_maxDepth=-1	(-1 indicates no cutof)
	-geno_minMM=-1.000000	(minimum fraction af major-minor bases)
	-minInd=0	(only keep sites if you call genotypes from this number of individuals)

	NB When writing the posterior the -postCutoff is not used
	NB geno_minDepth requires -doCounts
	NB geno_maxDepth requires -doCounts
```

Therefore, if we set `-doGeno 2`, genotypes are coded as 0,1,2, as the number of alternate alleles.
If we want to print the major and minor alleles as well then we set `-doGeno 3`.

To calculate the posterior probability of genotypes we need to define a model.
```
$ANGSD/angsd -doPost

...
-doPost	0	(Calculate posterior prob 3xgprob)
	1: Using frequency as prior
	2: Using uniform prior
	3: Using SFS as prior (still in development)
	4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an
...
```
`-doPost 1` uses the estimate per-site allele frequency as a prior for genotype proportions, assuming Hardy Weinberg Equilibrium.

When the assumption of HWE is not valid, you can use an estimate of the inbreeding coefficient, for instance calculated using [ngsF](https://github.com/fgvieira/ngsF) and using the option:
```
...
	-indFname	(null) (file containing individual inbreedcoeficients)
...
```
We will discuss later how to deal with inbred species.

For most cases, we want to restrict this analysis on a set of putative polymorphic sites (SNPs), as non-variable sites (across all samples) will not carry information regarding population structure or differentiation.
The rationale for assigning probabilities of being variable at each site is based on the estimation of the allele frequencies.

ANGSD has an option to estimate allele frequencies called `-doMaf`:

```
$ANGSD/angsd -doMaf
...
-doMaf	0 (Calculate persite frequencies '.mafs.gz')
	1: Frequency (fixed major and minor)
	2: Frequency (fixed major unknown minor)
	4: Frequency from genotype probabilities
	8: AlleleCounts based method (known major minor)
	NB. Filedumping is supressed if value is negative
...
Filters:
	-minMaf  	-1.000000	(Remove sites with MAF below)
	-SNP_pval	1.000000	(Remove sites with a pvalue larger)
	-rmTriallelic	0.000000	(Remove sites with a pvalue lower)
Extras:
	-ref	(null)	(Filename for fasta reference)
	-anc	(null)	(Filename for fasta ancestral)
	-eps	0.001000 [Only used for -doMaf &8]
	-beagleProb	0 (Dump beagle style postprobs)
	-indFname	(null) (file containing individual inbreedcoeficients)
	-underFlowProtect	0 (file containing individual inbreedcoeficients)
NB These frequency estimators requires major/minor -doMajorMinor

```

Therefore, the estimation of allele frequencies requires the specification of how to assign the major and minor alleles (if biallelic).
```
$ANGSD/angsd -doMajorMinor
...
        -doMajorMinor   0
        1: Infer major and minor from GL
        2: Infer major and minor from allele counts
        3: use major and minor from a file (requires -sites file.txt)
        4: Use reference allele as major (requires -ref)
        5: Use ancestral allele as major (requires -anc)
        -rmTrans: remove transitions 0
        -skipTriallelic 0
```

Finally, you need to specify which genotype likelihood model to use.
```
$ANGSD/angsd -GL
...
	-GL=0: 
	1: SAMtools
	2: GATK
	3: SOAPsnp
	4: SYK
	5: phys
	6: Super simple sample an allele type GL. (1.0,0.5,0.0)
	-trim		0		(zero means no trimming)
	-tmpdir		angsd_tmpdir/	(used by SOAPsnp)
	-errors		(null)		(used by SYK)
	-minInd		0		(0 indicates no filtering)

Filedumping:
	-doGlf	0
	1: binary glf (10 log likes)	.glf.gz
	2: beagle likelihood file	.beagle.gz
	3: binary 3 times likelihood	.glf.gz
	4: text version (10 log likes)	.glf.gz
```

If we are interested in looking at allele frequencies only for sites that are actually variable in our sample, we need to perform a SNP calling first.
There are two main ways to call SNPs using ANGSD with these options:
```
        -minMaf         0.000000        (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
```
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

## Principal Component Analysis

Back to our example, we are now performing a principal component analyses (PCA) without relying on called genotypes, but rather by taking their uncertainty into account.
More specifically, the next program we are going to use (ngsTools) takes as input genotype probabilities in binary format, so we need to specify `-doGeno 32`.
Also, we are using a HWE-based prior with `-doPost 1`.
Recalling also our choice for data filtering, our command line is:
```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL -r 11 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 60 -setMaxDepth 400 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
	-doGeno 32 -doPost 1 &> /dev/null
```
Unzip the results (but you cannot open it since it is in binary format)
```
gunzip Results/ALL.geno.gz
```

We are going to use `ngsCovar`, which estimates the covariance matrix between individuals based on genotype probabilities.
Then this matrix will be decomposed into principal components which will be investigated for population structure analyses.
Note that although `ngsCovar` can account for SNPs uncertanity, we find that it is faster to perform a light SNP filtering first (as we did using ```-SNP_pval 1e-3```) than using all sites.

If you type:
```
$NGSTOOLS/ngsPopGen/ngsCovar
```
you will see a list of possible options.
For instance, we need to define how many sites we have.
To retrieve such value, we can inspect the file with allele frequencies:
```
less -S Results/ALL.mafs.gz
NSITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $NSITES
```
In my case I have 4133 sites.

Now we can perform a PCA by estimating the covariance matrix first:
```
$NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/ALL.geno -outfile Results/ALL.covar -nind 30 -nsites $NSITES -call 0 -norm 0 &> /dev/null
```
with the options `-call 0` meaning that we do not perform genotype calling and `-norm 0` that we are not normalising by allele frequency.
The latter may give more weight to low frequency variants which are harder to estimate.

Look at the output file:
```
less -S Results/ALL.covar
```
which represents a symmetric matrix of NxN with N individuals.

Finally, we perform an eigenvector decomposition and plot the resulting map.
First we need to create a plink cluster-like file defining the labelling (population) for each sample.
For instance, in my case I can generate it with:
```
Rscript -e 'write.table(cbind(seq(1,30),rep(1,30),c(rep("LWK",10),rep("TSI",10),rep("PEL",10))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="Results/ALL.clst", quote=F)'
```
Run and plot:
```
Rscript $NGSTOOLS/Scripts/plotPCA.R -i Results/ALL.covar -c 1-2 -a Results/ALL.clst -o Results/ALL.pca.pdf
evince Results/ALL.pca.pdf
```
where the parameter `-c 1-2` specifies that we are plotting only the first and second component.
On the screen, you will see a series of numbers.
These are the percentage of explained variance for each component.

## Genetic distances

An alternative approach would be to compute genetic distances first, and then perform a Multi Dimensional Scaling (MDS) on those.
We are using [ngsDist](https://github.com/fgvieira/ngsDist) to estimate pairwise genetic distances from genotype probabilities.

Again, we run ANGSD to compute genotype psoterior probabilities assuming HWE and specifying the output format with `-doGeno`:
```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL -r 11 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 60 -setMaxDepth 400 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null
```

Record how many sites we retrieve (although this should be equal what found earlier):
```
NSITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $NSITES
```

For plotting purposes, we now create a file with labels indicating the population of interest for each sample.
```
Rscript -e 'cat(paste(rep(c("LWK","TSI","PEL"),each=10), rep(1:10, 3), sep="_"), sep="\n", file="Data/pops.label")'
cat Data/pops.label
```

With [ngsDist](https://github.com/fgvieira/ngsDist) we can compute pairwise genetic distances without relying on individual genotype calls.
```
$NGSTOOLS/ngsDist/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 30 -n_sites $NSITES -labels Data/pops.label -o Results/ALL.dist -n_threads 4 &> /dev/null
less -S Results/ALL.dist
```

We can visualise the pairwise genetic distances in form of a tree (in Newick format) using FastME:
```
$FASTME -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null
cat Results/ALL.tree
```
We can use some R packages to plot the resulting tree.
```
Rscript $NGSTOOLS/Scripts/plotTree.R Results/ALL.tree
evince Results/ALL.tree.pdf
```

One can even generate bootstrapped replicates of genetic distances in order to get a measure of confidence intervals.
This can be achieved by using the following options:
```
$NGSTOOLS/ngsDist/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 30 -n_sites $NSITES -labels Data/pops.label -o Results/ALL.boot.dist -n_threads 4 -n_boot_rep 20 -boot_block_size 20 &> /dev/null
```
which will generate 20 replicates by randomly sampling with replacemente blocks of 20 SNPs (since we called SNPs earlier, otherwise will indicate the genomic length).
Again, we can plot these results on a form of a tree:
```
$FASTME -D 21 -i Results/ALL.boot.dist -o Results/ALL.boot.tree -m b -n b &> /dev/null
Rscript $NGSTOOLS/Scripts/plotTreeBoots.R Results/ALL.boot.tree
evince Results/ALL.boot.tree.pdf
```
with boostrapped values shown on branches.
Please note that we specifiy 21 tree in FASTME as the output file consists of the original tree plus the 20 boostrapped replicates.

From these calculated genetics distances (not the boostrapped ones), we can also perform a MDS analysis and investigate the population genetic structure of our samples.
```
NSAMPLES=30
tail -n +3 Results/ALL.dist | head -n $NSAMPLES | Rscript --vanilla --slave $NGSTOOLS/Scripts/getMDS.R --no_header --data_symm -n 4 -m "mds" -o Results/ALL.mds &> /dev/null
less -S Results/ALL.mds
```
The first line gives the proportion of explained variance by each component, and the other lines specify the coordinates of each sample on each component.
You can plot the resulting components:
```
Rscript $NGSTOOLS/Scripts/plotMDS.R -i Results/ALL.mds -c 1-2 -a Results/ALL.clst -o Results/ALL.mds.pdf
evince Results/ALL.mds.pdf
```

## Admixture proportions

Admixture proportions can be estimated from genotype likelihoods using [NGSadmix](http://popgen.dk/angsd/index.php/NGSadmix), which requires input files in BEAGLE format.
This can be accomplished in ANGSD with the following command:
```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL -r 11 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 60 -setMaxDepth 400 -doCounts 1 \
        -GL 1 -doMajorMinor 4 -doMaf 1 -skipTriallelic 1 \
	-doGlf 2 -SNP_pval 1e-6 &> /dev/null
```
Note that, as an illustration, here we use a more stringent cutoff for SNP calling and we impose of the the alleles to be the reference one (-doMajorMinor 4).

Assuming we want to test for 3 ancestral components, admixture proportions can be obtained with:
```
K=3
$NGSADMIX -likes Results/ALL.beagle.gz -K $K -outfiles Results/ALL -P 4 -minMaf 0 &> /dev/null
```
and results are stored in the following files 'Results/ALL.qopt' and 'Results/ALL.fopt.gz', the former containing the inferred proportions for each individual.
You can plot the admixture proportions:
```
Rscript $NGSTOOLS/Scripts/plotAdmix.R -i Results/ALL.qopt -o Results/ALL.admix.pdf &> /dev/null
evince Results/ALL.admix.pdf
```

Inbreeding
------------------------------------------

For some of the previous analyses, we assumed that our populations were in HWE.
If this is not the case, you may want to adjust the prior by including an estimate of the inbreeding coefficient.
This can be achieved from genotype likelihoods using [ngsF](https://github.com/fgvieira/ngsF), part of ngsTools.

First, again we need to calculate genotype likelihoods in a format accepted by ngsF.
You also need to perform a SNP calling and ideally perform this analysis only on unlinked sites.
This can be achieved (for instance) by randomly sampling sites at a given distance (not shown here).
Please also note that here we are performing this analysis on single populations and not for the pooled sample.
```
for POP in LWK TSI PEL;
do
	echo $POP
	$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -out Results/$POP -r 11 \
        	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
        	-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
        	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        	-doGlf 3 -SNP_pval 1e-3 &> /dev/null
done
```

Inbreeding coefficients are calculated by first estimating reliable starting values and then performing a deep search.
We can do this using ngsF:
```
NSAMS=10
for POP in LWK TSI PEL;
do
	NSITES=`zcat Results/${POP}.mafs.gz | tail -n+2 | wc -l`
	echo $POP $NSAMS $NSITES

	# preliminary search
	zcat Results/$POP.glf.gz | $NGSTOOLS/ngsF/ngsF --n_ind $NSAMS --n_sites $NSITES --glf - --out Results/$POP.approx_indF --approx_EM --init_values u --n_threads 4 &> /dev/null

	zcat Results/$POP.glf.gz | $NGSTOOLS/ngsF/ngsF --n_ind $NSAMS --n_sites $NSITES --glf - --out Results/$POP.indF --init_values Results/$POP.approx_indF.pars --n_threads 4 &> /dev/null

	cat Results/$POP.indF

done
```

We can even use a routine in ngsF that computes 20 initial searches to find the best starting point (this will create a temporary directory in /home/scratch):
```
NSAMS=10
for POP in LWK TSI PEL;
do
	NSITES=`zcat Results/${POP}.mafs.gz | tail -n+2 | wc -l`
        echo $POP $NSAMS $NSITES

	zcat Results/$POP.glf.gz > Results/$POP.glf
	$NGSTOOLS/ngsF/ngsF.sh --n_ind $NSAMS --n_sites $NSITES --glf Results/$POP.glf --out Results/$POP.indF &> /dev/null

	cat Results/$POP.indF

done
```

If the estimate inbreeding coefficients are much larger than 0, you may want to correct your prior on genotypes and allele frequencies accordingly.
For instance, if we want to calculate genotype posterior probabilities taking into account the individual estimated inbreeding coefficients, we can use options `-doSaf 2 -indF ...indF`.
Please note that we also need to specify an ancestral sequence with `-anc`.
```
for POP in LWK TSI PEL;
do
	echo $POP
	$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/$POP.inbred -r 11 \
        	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
        	-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
        	-GL 1 -doMajorMinor 1 -doMaf -1 -skipTriallelic 1 \
        	-SNP_pval 1e-3 \
        	-doGeno 20 -doPost 1 \
		-doSaf 2 -indF Results/$POP.indF &> /dev/null
done
```
Note that, as an illustration here, the option `-doMaf 1` will not print out a file with the allele frequencies.
Other analyses involving allele frequencies, as we will see later, can be done incorporating individual inbreeding coefficients.

Additionally, one can estimate per-individual inbreeding tracts via a two-state Hidden Markov Model (HMM) using [ngsF-HMM](https://github.com/fgvieira/ngsF-HMM).
You again need genotype likelihoods generated with `-doGlf 3` option.
```
NSAMS=10
for POP in LWK TSI PEL;
do
        echo $POP

	# uncomment if you don't have these files already
	$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -out Results/$POP -r 11 \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
                -doGlf 3 -SNP_pval 1e-3 &> /dev/null

	zcat Results/$POP.glf.gz > Results/$POP.glf

	NSITES=$((`zcat Results/$POP.mafs.gz | wc -l`-1))

	echo $POP $NSITES

	$NGSTOOLS/ngsF-HMM --geno Results/$POP.glf --loglkl --n_ind $NSAMS --n_sites $NSITES --freq r --freq_est 2 --indF 0.1,0.1 --out Results/$POP.F-HMM --n_threads 4 &> /dev/null
done
```


Summary statistics using ANGSD
-----------------------------------

Nucleotide diversity indexes and measures of population differentiation can be estimated taking data uncertainty into account both with ANGSD and ngsTools.
First, we show how to estimate such summary statistics using ANGSD.

## The Site Frequency Spectrum (SFS)

One of the most important aspect of data analysis for population genetics is the estimate of the Site Frequency Spectrum (SFS). 
SFS records the proportions of sites at different allele frequencies. It can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state. 
SFS is informative on the demography of the population or on selective events (when estimated at a local scale).

We use ANGSD to estimate SFS using on example dataset, using the methods described [here](http://www.ncbi.nlm.nih.gov/pubmed/22911679).
Details on the implementation can be found [here](http://popgen.dk/angsd/index.php/SFS_Estimation).
Briefly, from sequencing data one computes genotype likelihoods (as previously described). 
From these quantities ANGSD computes posterior probabilities of Sample Allele Frequency (SAF), for each site. 
Finally, an estimate of the SFS is computed.

These steps can be accomplished in ANGSD using `-doSaf 1/2` options and the program `realSFS`.

```
$ANGSD/angsd -doSaf
...
-doSaf		0
	1: perform multisample GL estimation
	2: use an inbreeding version
	3: calculate genotype probabilities
	4: Assume genotype posteriors as input (still beta) 
	-doThetas		0 (calculate thetas)
	-underFlowProtect	0
	-fold			0 (deprecated)
	-anc			(null) (ancestral fasta)
	-noTrans		0 (remove transitions)
	-pest			(null) (prior SFS)
	-isHap			0 (is haploid beta!)
NB:
	  If -pest is supplied in addition to -doSaf then the output will then be posterior probability of the sample allelefrequency for each site
```

The SFS is typically computed for each population separately.
We need to slightly modify the filtering options as now each population has 20 samples. 
So now we set `-minInd 20 -setMinDepth 20 -setMaxDepth 200`.
Also, we want to estimate the unfolded SFS and we use a putative ancestral sequence to polarise our alleles (to ancestral and derived states).

We cycle across all populations:
```
for POP in LWK TSI PEL
do
	echo $POP
	$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/$POP -r 11 \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
		-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
		-GL 1 -doSaf 1 &> /dev/null
done
```

Have a look at the output file.
```
$ANGSD/misc/realSFS print Results/PEL.saf.idx | less -S
```
These values represent the sample allele frequency likelihoods at each site.
So the first value (after the chromosome and position columns) is the likelihood of having 0 copies of the derived allele, the second indicates the probability of having 1 copy and so on.
Note that these values are in log format and scaled so that the maximum is 0.

The next step would be to use these likelihoods and estimate the overall SFS.
This is achieved by the program `realSFS`.
```
$ANGSD/misc/realSFS
-> ---./realSFS------
	-> EXAMPLES FOR ESTIMATING THE (MULTI) SFS:

	-> Estimate the SFS for entire genome??
	-> ./realSFS afile.saf.idx 

	-> 1) Estimate the SFS for entire chromosome 22 ??
	-> ./realSFS afile.saf.idx -r chr22 

	-> 2) Estimate the 2d-SFS for entire chromosome 22 ??
	-> ./realSFS afile1.saf.idx  afile2.saf.idx -r chr22 

	-> 3) Estimate the SFS for the first 500megabases (this will span multiple chromosomes) ??
	-> ./realSFS afile.saf.idx -nSites 500000000 

	-> 4) Estimate the SFS around a gene ??
	-> ./realSFS afile.saf.idx -r chr2:135000000-140000000 

	-> Other options [-P nthreads -tole tolerence_for_breaking_EM -maxIter max_nr_iterations -bootstrap number_of_replications]

	-> See realSFS print for possible print options
	-> Use realSFS print_header for printing the header
	-> Use realSFS cat for concatenating saf files

	->------------------
	-> NB: Output is now counts of sites instead of log probs!!
	-> NB: You can print data with ./realSFS print afile.saf.idx !!
	-> NB: Higher order SFS's can be estimated by simply supplying multiple .saf.idx files!!
	-> NB: Program uses accelerated EM, to use standard EM supply -m 0 
	-> Other subfunctions saf2theta, cat, check

```

This command will estimate the SFS for each population:
```
for POP in LWK TSI PEL
do
        echo $POP
        $ANGSD/misc/realSFS Results/$POP.saf.idx -P 4 2> /dev/null > Results/$POP.sfs
done
```
The output will be saved in Results/POP.sfs files.

You can now have a look at the output file, for instance for the African (LWK) samples:
```
cat Results/LWK.sfs
```
The first value represent the expected number of sites with derived allele frequency equal to 0, the second column the expected number of sites with frequency equal to 1 and so on.

How many values do you expect?
```
awk -F' ' '{print NF; exit}' Results/LWK.sfs 
```
Indeed this represents the unfolded spectrum, so it has 2N+1 values with N diploid individuals.
Please note that this maximum likelihood estimation of the SFS should be performed at the whole-genome level to have enough information for the algorithm to converge.
However, for practical reasons, here we could not use large genomic regions.

You can plot the SFS for each pop using this simple R script.
```
Rscript $NGSTOOLS/Scripts/plotSFS.R Results/LWK.sfs-Results/TSI.sfs-Results/PEL.sfs LWK-TSI-PEL 0 Results/ALL.sfs.pdf
evince Results/ALL.sfs.pdf
```
where the second parameter is the labels for each population and the third parameter 0 is a boolean whether the SFS should be folded (1) or not (0).
Since we provide an ancestral sequence here we can plot the unfolded spectrum.

If, on the other hand, we don't provide an ancestral sequence, we can use the reference sequence to polarise our results and the fold the spectrum afterwards.
This can be achieved, for instance, with the following lines:
```
for POP in LWK TSI PEL
do
        echo $POP
        $ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $REF -out Results/${POP}.ref -r 11 \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doSaf 1 &> /dev/null
	$ANGSD/misc/realSFS Results/$POP.ref.saf.idx -P 4 2> /dev/null > Results/$POP.ref.sfs
	
done
Rscript $NGSTOOLS/Scripts/plotSFS.R Results/LWK.ref.sfs-Results/TSI.ref.sfs-Results/PEL.ref.sfs LWK-TSI-PEL 1 Results/ALL.ref.sfs.pdf
evince Results/ALL.ref.sfs.pdf
```

It is sometimes convenient to generate bootstrapped replicates of the SFS, by sampling with replacements genomic segments.
This could be used for instance to get confidence intervals when using the SFS for demographic inferences.
This can be achieved in ANGSD using:
```
$ANGSD/misc/realSFS Results/PEL.saf.idx -bootstrap 10  2> /dev/null > Results/PEL.boots.sfs
cat Results/PEL.boots.sfs
```
This command may take some time.
The output file has one line for each boostrapped replicate.

It is very useful to estimate a multi-dimensional SFS, for instance the joint SFS between 2 populations (2D).
This can be used for making inferences on their divergence process (time, migration rate and so on) or even used as prior information for estimating population genetic differentiation.

An important issue when doing this is to be sure that we are comparing the exactly same corresponding sites between populations.
ANGSD does that automatically and considers only a set of overlapping sites.
The 2D-SFS between all populations and PEL, for instance, are computed with:
```
for POP in LWK TSI
do
	$ANGSD/misc/realSFS -P 4 Results/$POP.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/$POP.PEL.sfs
done
```
while beteen LWK and TSI is:
```
$ANGSD/misc/realSFS -P 4 Results/LWK.saf.idx Results/TSI.saf.idx 2> /dev/null > Results/LWK.TSI.sfs
```

The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] [0,2] ... [1,0] [1,1] ... and so on until [20, 20] (assuming 10 diploids per population).
```
less -S Results/LWK.PEL.sfs
```
You can plot it, but you need to define how many individuals you have per population (please note that this script assumes that your samples are diploid):
```
Rscript $NGSTOOLS/Scripts/plot2DSFS.R Results/LWK.PEL.sfs LWK-PEL 10-10
evince Results/LWK.PEL.sfs.pdf
```
This script masks the non-variant cells and it is based on unfolded data.
At this stage, it should not be used for folded data.

You can even estimate SFS with higher order of magnitude (more than 2 populations).
This command may take some time (and should be run on a much larger genomic region).
```
$ANGSD/misc/realSFS -P 4 Results/LWK.saf.idx Results/TSI.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/LWK.TSI.PEL.sfs
```

There is an (unsupported) routine to generate a plot for a 3D-SFS.
You need install the following program:
```
cd ~/Software # or go to wherever you installed your programs
git clone https://github.com/lpmdiaz/pop3D
cd pop3D
mkdir bin
make
```
Now go back to where you are running this tutorial and add:
```
POP3D=~/Software/pop3D
```
You can now plot the 3D-SFS with the lines (where 10 is the number of individuals per population and 1 is the window to consider, as we are not doing a sliding-windows scan):
```
$POP3D/bin/parse3Dsfs Results/LWK.TSI.PEL.sfs Results/ALL.parsed3Dsfs 10 10 10 1
Rscript $POP3D/scripts/plot3Dsfs.R Results/ALL.parsed3Dsfs Results/ALL.3Dsfs
evince Results/ALL.3Dsfs.pdf
```
The interactive session is not fully supported.
Please note that frequencies are in log scale.


## Population genetic differentiation (FST/PBS)

Here we are going to calculate allele frequency differentiation using the PBS (population branch statistic) and FST metrics.
Again, we can achieve this by avoid genotype calling using ANGSD directly from the sample allele frequencies likelihoods.
Note that we are using the 2D-SFS calculated above as prior information for the joint allele frequency probabilities at each site.
As an illustration, PEL is our target population, while LWK and TSI are reference populations.

Specifically, we are computing a slinding windows scan, with windows of 50kbp and a step of 10kbp.
This can be achieved using the following commands.
This first command computes per-site FST indexes:
```
$ANGSD/misc/realSFS fst index Results/LWK.saf.idx Results/TSI.saf.idx Results/PEL.saf.idx -sfs Results/LWK.TSI.sfs -sfs Results/LWK.PEL.sfs -sfs Results/TSI.PEL.sfs -fstout Results/PEL.pbs -whichFST 1 &> /dev/null
```
and you can have a look at their values:
```
$ANGSD/misc/realSFS fst print Results/PEL.pbs.fst.idx | less -S
```
where columns are: chromosome, position, (a), (a+b) values for the three FST comparisons, where FST is defined as a/(a+b).
Note that FST on multiple SNPs is calculated as sum(a)/sum(a+b).

Then, the next command performs a sliding-window analysis
```
$ANGSD/misc/realSFS fst stats2 Results/PEL.pbs.fst.idx -win 50000 -step 10000 -whichFST 1 > Results/PEL.pbs.txt
```

You can have a look at the output file:
```
less -S Results/PEL.pbs.txt 
```
The header is:
```
region	chr	midPos	Nsites	Fst01	Fst02	Fst12	PBS0	PBS1	PBS2
```
Where are interested in the column `PB2` which gives the PBS values assuming PEL (coded here as 2) being the target population.
Note that negative PBS and FST values are equivalent to 0.
We are also provided with the individual FST values.
You can see that high values of PBS2 are indeed associated with high values of both Fst02 and Fst12 but not Fst01.

Please note that if you give only 2 populations in input, only the FST will be calculated.
Therefore, FST between a pair of populations can calculated by repeating the steps above by using only 2 .saf.idx input files.
For instance, FST between LWK and TSI is achieved with:
```
$ANGSD/misc/realSFS fst index Results/LWK.saf.idx Results/TSI.saf.idx -sfs Results/LWK.TSI.sfs -fstout Results/LWK.TSI -whichFST 1
$ANGSD/misc/realSFS fst stats2 Results/LWK.TSI.fst.idx -win 50000 -step 10000 -whichFST 1 > Results/LWK.TSI.fst.txt
less -S Results/LWK.TSI.fst.txt
```


## Nucleotide diversity

You may be interested in assessing levels of nucleotide diversity within a particular population.
Again, we can achieve this using ANGSD by estimating levels of diversity without relying on called genotypes.

The procedure is similar to what done for FST/PBS, and the SFS is again used as a prior to compute allele frequencies probabilities. 
From these quantities, expectations of various diversity indexes are compute.
This can be achieved using the following pipeline, assuming we are computing such indexes for all populations (but separately).

First we compute the allele frequency posterior probabilities and associated statistics (-doThetas) using the SFS as prior information (-pest)
```
for POP in LWK TSI PEL
do
	echo $POP
	$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/$POP -r 11 \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doSaf 1 \
		-doThetas 1 -pest Results/$POP.sfs &> /dev/null
done
```

Then we need to index these file and perform a sliding windows analysis using a window length of 50kbp and a step size of 10kbp, as an example (note that the option -cChr is twice the number of diploid individuals).
```
for POP in LWK TSI PEL
do
	echo $POP
	$ANGSD/misc/thetaStat do_stat Results/$POP.thetas.idx &> /dev/null
	# perform a sliding-window analysis
	$ANGSD/misc/thetaStat do_stat Results/$POP.thetas.idx -win 50000 -step 10000 -outnames Results/$POP.thetas &> /dev/null
done
```
Values in this output file are the sum of the per-site estimates for the whole window.
For instance:
```
less -S Results/PEL.thetas.pestPG
```

Finally, you may also be interested in estimating allele frequencies for single SNPs of interest.
In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option.
Assume that these are the SNPs we are interested in (chromosome and genomic position 1-based):
- 11 61627960 <br>
- 11 61631510 <br>
- 11 61632310 <br>
- 11 61641717 <br>
- 11 61624414 <br>
- 11 61597212 <br>

The file with these positions need to be formatted as (chromosome positions).
```
> Data/snps.txt
echo 11 61627960 >> Data/snps.txt
echo 11 61631510 >> Data/snps.txt
echo 11 61632310 >> Data/snps.txt
echo 11 61641717 >> Data/snps.txt
echo 11 61624414 >> Data/snps.txt
echo 11 61597212 >> Data/snps.txt
```
We need to index this file in order for ANGSD to process it.
```
$ANGSD/angsd sites index Data/snps.txt
```

We are interested in calculating the derived allele frequencies, so are using the ancestral sequence to polarise the alleles with '-doMajorMinor 5'.
Note that here we change the filtering (more relaxed) since we are interested in outputting all sites.
```
for POP in LWK TSI PEL
do
        echo $POP
        $ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/$POP -r 11 \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -sites Data/snps.txt &> /dev/null
done
zcat Results/LWK.mafs.gz Results/TSI.mafs.gz Results/PEL.mafs.gz
```

Summary statistics using ngsTools
-----------------------------------

VERY IMPORTANT NOTE: we recommend the use of ANGSD to calculate summary statistics. 
While the rationale behind these estimations are very similar with ANGSD and ngsTools, ANGSD provides a much better implementation of them.

Most of the summary statistics can be now estimated using ANGSD using the commands described above.
For the sake of completeness, here we show how we can perform similar analyses using ngsTools, although ANGSD may be faster and require less memory than ngsTools.
Note that the methods behind how ANGSD and ngsTools estimate such quantities are very similar (but results may differ slightly).

## Population genetic differentiation (FST)

We first compute the sample allele frequency likelihoods using ANGSD.
```
for POP in LWK TSI PEL
do
        echo $POP
        $ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/$POP -r 11 \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doSaf 1 &> /dev/null
done
```

Assume we are interested in calculating FST between TSI and PEL.
Unlike in ANGSD, here we first need to get the subset of overlapping sites (unfiltered for both populations), and we will store this information in the file `Data/intersect.txt`.
```
	$ANGSD/misc/realSFS print Results/TSI.saf.idx Results/PEL.saf.idx | cut -f 1-2 > Data/intersect.txt
	$ANGSD/angsd sites index Data/intersect.txt
	NSITES=`wc -l Data/intersect.txt | cut -f 1 -d " "`
	echo $NSITES
```

We then compute the sample allele frequency likelihoods only for the overlapping (valid) sites.
```
for POP in TSI PEL
do
        echo $POP
        $ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/$POP -r 11 \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doSaf 1 \
		-sites Data/intersect.txt &> /dev/null
done
```

We then estimate the 2D-SFS to be used as prior using ngsTools.
```
NSITES=`wc -l Data/intersect.txt | cut -f 1 -d " "` # if not already done
zcat Results/TSI.saf.gz > Results/TSI.saf
zcat Results/PEL.saf.gz > Results/PEL.saf
$NGSTOOLS/ngsPopGen/ngs2dSFS -postfiles Results/TSI.saf Results/PEL.saf -outfile Results/TSI.PEL.2dsfs -nind 10 10 -nsites $NSITES -maxlike 1 -relative 1
```
The output is a matrix giving the proportion of sites with a given joint allele frequency.

ngsTools implements a very simple estimator of the SFS, where the most likely joint allele frequency is recorded at each site.
This tends to overestimate singletons, for instance.
If you have enough sites, then it is recommended to use the 2D-SFS estimated in ANGSD instead.
Recalling what previously shown, the command is (note that -sites option is kept for consistency but it is not really needed here):
```
$ANGSD/misc/realSFS -P 4 Results/TSI.saf.idx Results/PEL.saf.idx -sites Data/intersect.txt 2> /dev/null > Results/TSI.PEL.sfs
```
If you want to use this 2D-SFS calculated by ANGSD, then you need to convert it into the input file for ngsTools (note that you need to specify the number of individuals for both populations):
```
Rscript $NGSTOOLS/Scripts/SFSangsd2tools.R Results/TSI.PEL.sfs 10 10 > Results/TSI.PEL.angsd.2dsfs
```

We can now calculate per-site FST values (using for instance the 2D-SFS calculated by ANGSD).
```
$NGSTOOLS/ngsPopGen/ngsFST -postfiles Results/TSI.saf Results/PEL.saf -priorfile Results/TSI.PEL.angsd.2dsfs -nind 10 10 -nsites $NSITES -outfile Results/TSI.PEL.fst
```
The output has the following header: a, (a+b), correcting factor (you can ignore this), FST, probability of being variable.
Note that FST is equal (ignoring the correcting factor) to a/(a+b), and multiple-sites estimates are sum(a)/sum(a+b).
Also note that negative FST values mean FST equal to 0. 
Very large negative or positive per-site values are simply due to numerical corrections and they tend to appear when a site is highly unlikely to be variable.
These will not affect multiple-sites estimates and can be safely ignored.

These example scripts will produce a plot and text file with sliding windows values (size of 50kbp, step of 10kbp).
Note that, as an illustration, we also filter out sites with a probability of being variable less than 0.90, although other options are valid too.
```
Rscript $NGSTOOLS/Scripts/plotFST.R -i Results/TSI.PEL.fst -o Results/TSI.PEL.scan.fst -p Data/intersect.txt -w 50000 -s 10000 -t 0.90
less -S Results/TSI.PEL.scan.fst.txt
evince Results/TSI.PEL.scan.fst.pdf
```

## Nucleotide diversity

We use the sample allele frequency probabilities (.saf files) calculated on previous steps and estimate the marginal SFS to be used as priors.
Again, you might have already calculated these SFS and the -sites option is not really necessary but it is kept for consistency.
```
for POP in LWK TSI PEL
do
	echo $POP
	$ANGSD/misc/realSFS -P 4 Results/$POP.saf.idx -sites Data/intersect.txt > Results/$POP.sfs &> /dev/null
done
```

Assuming we are interested in TSI and PEL, we can now calculate some summary statistics, namely number of segregating sites, expected heterozygosity, number of fixed differences and dxy.
Please note that the latter 2 statistics have not been properly tested and they are shown to be 
For instance, dxy been shown to be over-estimated and should be used only for inspecting the distribution and not to make inferences based on its absolute values.
In case you want to estimate dxy, you can find in `ngsTools/ngsPopGen/scripts` folder a perl script written by [Nagarjun Vijay](https://lsa.umich.edu/eeb/people/postdoctoral-fellows/nagarju.html) and a R script written by [Joshua Penalba][https://joshuapenalba.com/] which calculate Dxy from ANGSD allele frequency files. Please see that script for additional help.

We can estimate summary statistics for two populations with:
```
zcat Results/TSI.saf.gz > Results/TSI.saf
zcat Results/PEL.saf.gz > Results/PEL.saf
NSITES=`wc -l Data/intersect.txt | cut -f 1 -d " "` # if not already done
$NGSTOOLS/ngsPopGen/ngsStat -npop 2 -postfiles Results/TSI.saf Results/PEL.saf -nsites $NSITES -nind 10 10 -outfile Results/TSI.PEL.stats.txt
```
while if we are interested in only one population the command would be:
```
$NGSTOOLS/ngsPopGen/ngsStat -npop 1 -postfiles Results/PEL.saf -nsites $NSITES -nind 10 -outfile Results/PEL.stats.txt
```

Have a look at the output files:
```
less -S Results/TSI.PEL.stats.txt
less -S Results/PEL.stats.txt
```
The output file (for a 2-populations analysis) has the following header: position start, position end, segregating sites (pop 1), heterozygosity (pop 1), segregating sites (pop 2), heterozygosity (pop 2), fixed differences, dxy.
If only 1 population is analysed, this is the header: position start, position end, segregating sites, heterozygosity.

These example scripts will produce a plot and text file with sliding windows values, for 2 populations and only 1 population.
```
Rscript $NGSTOOLS/Scripts/plotSS.R -i Results/TSI.PEL.stats.txt -p Data/intersect.txt -o Results/TSI.PEL.scan.stats -n TSI-PEL -w 50000 -s 10000
Rscript $NGSTOOLS/Scripts/plotSS.R -i Results/PEL.stats.txt -p Data/intersect.txt -o Results/PEL.scan.stats -n PEL -w 50000 -s 10000
less -S Results/TSI.PEL.scan.stats.txt
less -S Results/PEL.scan.stats.txt
evince Results/TSI.PEL.scan.stats.pdf
evince Results/PEL.scan.stats.pdf
```

Additional help
-----------------------------

Further information and more example can be found at the individual web page for [ANGSD](http://popgen.dk/wiki/index.php/ANGSD).
Also, [ANGSD-wrapper](https://github.com/mojaveazure/angsd-wrapper) is a utility to run ANGSD and ngsTools developed by the [Ross-Ibarra Lab](http://www.rilab.org/) at UC Davis.
Please contact [me](https://www.imperial.ac.uk/people/m.fumagalli) for questions and feedback on this tutorial.


