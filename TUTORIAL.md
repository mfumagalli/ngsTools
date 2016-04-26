
A short tutorial for some basic analyses using ngsTools (plus ANGSD and NGSadmix)
===============

Please be sure you are using the most updated version of ngsTools. In doubt please run: 

    git pull
    git submodule update
    make clean
    make
	
inside the ngsTools directory.

Compatibility issues with ANGSD
---------------

We are currently working on updating ngsTools with the latest version of ANGSD! 
Check back soon or subscribe to our google mailing list.

Settings
----------

In this tutorial we will be using several programs including ngsTools, ANGSD and NGSadmix to perform population genetics analyses from low-depth sequencing data.
Please note that [ANGSD](http://popgen.dk/angsd/index.php/Main_Page#Overview) and [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) have not been developed by us and therefore questions on these tools should be addressed to their Authors.
However, given the utility of such tools, we felt the need to include them to present a more comprehensive view on the application of this probabilistic approach to process NGS data in population genetics.
Finally, we are using [SAMtools](http://samtools.sourceforge.net/) for indexing files, [FastMe](http://www.atgc-montpellier.fr/fastme/) for plotting trees and [R](https://www.r-project.org/) for manipulating and plotting results. This tutorial has been tested with SAMtools version 1.3.1, FastME version 2.1.4, R version 3.2.5.

First, set directories to all required programs depending on where you installed them:

    NGSTOOLS=/data/Software/ngsTools
    ANGSD=$NGSTOOLS/angsd
    NGSADMIX=/data/data/Software/NGSadmix/NGSadmix

    SAMTOOLS=/data/data/Software/samtools-1.3/samtools
    FASTME=/data/data/Software/fastme-2.1.4/binaries/fastme-2.1.4-linux64

Second, create all directories where you will be working:

    mkdir Tutorial
    cd Tutorial
    mkdir Data
    mkdir Results


Data
----------

As an illustration, we will use 60 BAM files of human samples (of African, European, and Native American descent), a reference genome, and a putative ancestral sequence.
BAM files have been downsampled to a mean depth of around 4X.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).
All data is publicly available.

A pipeline to retrieve such data is provided [here](https://github.com/mfumagalli/Weggis/blob/master/Files/data.sh).
You need to have 'samtools' (tested with version 1.3.1), 'bgzip' (tested with version 1.2.1-69-gb79f40a), 'gunzip' and 'wget' installed in your /usr/bin to run this.
Otherwise edit the first line of 'Scripts/data.sh' file to set the appropriate paths.
```
    cp $NGSTOOLS/Files/*.txt .
    bash $NGSTOOLS/Scripts/data.sh
```

Now we have 60 BAM files at low depth, a reference and an ancestral sequence in FASTA format.
```
	cat Data/download.log
```
and the list with BAM files has been written to 'ALL.bamlist'.
```
	cat ALL.bamlist
```

As a note for the general use, in case an ancestral sequence is not available, analyses on FST, PCA, nucleotide diversity (but not the number of fixed differences) can be carried out using the reference sequence to polarise your data. 
Please be aware that, under this scenario, some quantities (e.g. the unfolded joint site frequency spectrum) will be nonsense.

Please note that, since we are randomly subsampling reads here, your results in this tutorial may (slightly) differ from what written here. 

Basic filtering using ANGSD
----------------------

Here we will use ANGSD to analyse our data, for filtering and for generating files as input for ngsTools.
To see a full list of options in ANGSD type:
```
$ANGSD/angsd
```
and you should see something like:
```
...
Overview of methods:
	-GL		Estimate genotype likelihoods
	-doCounts	Calculate various counts statistics
	-doAsso		Perform association study
	-doMaf		Estimate allele frequencies
	-doError	Estimate the type specific error rates
	-doAncError	Estimate the errorrate based on perfect fastas
	-doHWE		Est inbreedning per site
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
	-uniqueOnly	0	Discards reads that doesnt map uniquely
	-show		0	Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
	-minMapQ	0	Discard reads with mapping quality below
	-minQ		13	Discard bases with base quality below
	-trim		0	Number of based to discard at both ends of the reads
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
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 1000 &> /dev/null
```
As input we give the list of BAM files with option `-b` and then specify the references sequence with `-ref` and the prefix for output files with `-out`.
Additionally, ```-C 50``` reduces the effect of reads with excessive mismatches, while ```-baq 1``` computes base alignment quality as explained here ([BAQ](http://samtools.sourceforge.net/mpileup.shtml)) to rule out false SNPs close to INDELS, and ```-trim 0``` means that we are not trimming the ends of reads.
With ```-minMapQ 20``` we filter out reads with low mapping quality.
Finally, ```-maxDepth 1000``` means that all sites with depth equal or greater than 1000 will be binned together, and ```-P 4``` means that I am using 4 threads.

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
wc -l Results/ALL.qc.depthSample # 60 Results/ALL.qc.depthSample
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
-minMap 20 | minimum mapping quality of 20 |
-minQ 10 | minimum base quality of 10 |
-minInd 40 | use only sites with data from at least 40 individuals |
-setMinDepth 40 | minimum total depth |
-setMaxDepth 240 | maximum total depth |

Please note that ANGSD can also compute more sophisticated metrics to filter out SNPs, as described [here](http://popgen.dk/angsd/index.php/SnpFilters), mostly based on:

* strand bias
* deviation from HWE
* quality score bias

Moreover, additional filtering should be considered.
For instance transitions (A<->G, C<->T) are more likely than transversions, so we expect the ts/tv ratio to be greater than 0.5.
However, we are not discussing this additional filtering options in this tutorial.


Population structure
---------------------------------------

Suppore we want to investigate the population structure our samples: PEL (Peruvians), TSI (Europeans), LWK (Africans), and CHB (East Asians).
One solution would be to perform a Principal Component Analysis (PCA) or a Multidimensional Scaling (MDS) or some clustering based on genetic distances among samples.
We are here showing how to perform these analyses using ngsTools/ANGSD in case of low-depth data.

To do this, we first need to assign genotype probabilities at each site for each individual.
The specific option in ANGSD is `-doGeno`.

```
$ANGSD/angsd -doGeno
...
-doGeno 0
        1: write major and minor
        2: write the called genotype encoded as -1,0,1,2, -1=not called
        4: write the called genotype directly: eg AA,AC etc
        8: write the posterior probability of all possible genotypes
        16: write the posterior probability of called gentype
        32: write the posterior probability of called gentype as binary
        -> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
        -postCutoff=0.333333 (Only genotype to missing if below this threshold)
        -geno_minDepth=-1       (-1 indicates no cutof)
        -geno_maxDepth=-1       (-1 indicates no cutof)
        -geno_minMM=-1.000000   (minimum fraction af major-minor bases)
        -minInd=0       (only keep sites if you call genotypes from this number of individuals)

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
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
...
```
`-doPost 1` uses the estimate per-site allele frequency as a prior for genotype proportions, assuming Hardy Weinberg Equilibrium.
When the assumption of HWE is not valid, you can use an estimate of the inbreeding coefficient, for instance calculated using [ngsF](https://github.com/fgvieira/ngsF).

For most cases, we want to restrict this analysis on a set of putative polymorphic sites (SNPs), as non-variable sites (across all samples) will not carry information regarding population structure or differentiation.
The rationale for assigning probabilities of being variable at each site is based on the estimation of the allele frequencies.

ANGSD has an option to estimate allele frequencies called `-doMaf`:

```
$ANGSD/angsd -doMaf
...
-doMaf  0 (Calculate persite frequencies '.mafs.gz')
        1: Frequency (fixed major and minor)
        2: Frequency (fixed major unknown minor)
        4: Frequency from genotype probabilities
        8: AlleleCounts based method (known major minor)
        NB. Filedumping is supressed if value is negative
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
Filters:
        -minMaf         -1.000000       (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
        -rmTriallelic   0.000000        (Remove sites with a pvalue lower)
Extras:
        -ref    (null)  (Filename for fasta reference)
        -anc    (null)  (Filename for fasta ancestral)
        -eps    0.001000 [Only used for -doMaf &8]
        -beagleProb     0 (Dump beagle style postprobs)
        -indFname       (null) (file containing individual inbreedcoeficients)
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
        -trim           0               (zero means no trimming)
        -tmpdir         angsd_tmpdir/   (used by SOAPsnp)
        -errors         (null)          (used by SYK)
        -minInd         0               (0 indicates no filtering)

Filedumping:
        -doGlf  0
        1: binary glf (10 log likes)    .glf.gz
        2: beagle likelihood file       .beagle.gz
        3: binary 3 times likelihood    .glf.gz
        4: text version (10 log likes)  .glf.gz
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
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 10 -minInd 40 -setMinDepth 40 -setMaxDepth 240 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3\
        -doGeno 32 -doPost 1 &> /dev/null
```
Unzip the results (but you cannot open it since it is in binary format)
```
gunzip Results/ALL.geno.gz
```

We are going to use `ngsCovar`, which estimates the covariance matrix between individuals based on genotype probabilities.
Then this matrix will be decomposed into principal componenets which will be investigated for population structure analyses.
Note that although `ngsCovar` can account for SNPs uncertanity, we find that it is faster to perform a light SNP filtering first (as we did using ```-SNP_pval 1e-3```) than using all sites.

If you type
```
$NGSTOOLS/ngsPopGen/ngsCovar
```
you will see a list of possible options.

For instance, we need to define how many sites we have.
To retrieve such value, we can inspect the file with allele frequencies:
```
less -S Results/ALL.mafs.gz
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```

Now we can perform a PCA by estimating the covariance matrix first:
```
$NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/ALL.geno -outfile Results/ALL.covar -nind 80 -nsites $N_SITES -call 0 -norm 0 &> /dev/null
```
with the options `-call 0` meaning that we do not perform genotype calling and `-norm 0` that we are not normalising by allele frequency.
The latter may give more weight to low frequency variants which are harder to estimate.

Look at the output file:
```
less -S Results/ALL.covar
```
which represents a matrix of NxN with N individuals giving the covariance.
Note that this matrix is symmetric.

Finally, we perform an eigenvector decomposition and plot the resulting map:
```
# create a cluster-like file defining the labelling (population) for each sample
Rscript -e 'write.table(cbind(seq(1,80),rep(1,80),c(rep("LWK",20),rep("TSI",20),rep("CHB",20),rep("PEL",20))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="Results/ALL.clst", quote=F)'
# run and plot
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
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 10 -minInd 40 -setMinDepth 40 -setMaxDepth 240 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null
```

Record how many sites we retrieve.
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```

For plotting purposes, we now create a file with labels indicating the population of interest for each sample.
```
Rscript -e 'cat(paste(rep(c("LWK","TSI","CHB","PEL"),each=20), rep(1:20, 4), sep="_"), sep="\n", file="Data/pops.label")'
cat Data/pops.label
```

With [ngsDist](https://github.com/fgvieira/ngsDist) we can compute pairwise genetic distances without relying on individual genotype calls.
```
$NGSTOOLS/ngsDist/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 80 -n_sites $N_SITES -labels Data/pops.label -o Results/ALL.dist -n_threads 4 &> /dev/null
less -S Results/ALL.dist
```

We can visualise the pairwise genetic distances in form of a tree (in Newick format) using FastME:
```
$FASTME -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null
cat Results/ALL.tree
```
We can use some R packages to plot the resulting tree.
```
Rscript -e 'library(ape); library(phangorn); pdf(file="Results/ALL.tree.pdf"); plot(read.tree("Results/ALL.tree"), cex=0.5); dev.off();' &> /dev/null
evince Results/ALL.tree.pdf
```

From these distances, we can also perform a MDS analysis and investigate the population genetic structure of our samples.
```
N_SAMPLES=80
tail -n +3 Results/ALL.dist | head -n $N_SAMPLES | Rscript --vanilla --slave $NGSTOOLS/Scripts/getMDS.R --no_header --data_symm -n 4 -m "mds" -o Results/ALL.mds &> /dev/null
less -S Results/ALL.mds
```
The first line gives the proportion of explained variance by each component, and the other lines specify the coordinates of each sample on each component.

OLD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This is just a possible arbitrary filtering setting. Please refer to ANGSD website for more options. 
In this scenario, we set a minimum mapping quality of 30 and base quality of 20. 
We furthermore filter out sites where we have data for less than 6 individuals and a global depth less than 20 or greater than 100.
We analyse only the first chromosome, using the `-r 1:` options.
Options `-GL 1 -doSaf 1 -doMaf 1 -doMajorMinor 1` specify how we calculate genotype likelihoods, major and minor alleles, and allele frequencies.
Finally, we performed a SNP calling by removing sites where the estimated minor allele frequency is less than 0.05 (equal to the frequency of singletons with 10 diploids).
Despite all methods work even when all sites are considered, excluding sites which are clearly not variable can help reduce computational time, memory space, and prevent some numerical instabilities.
We do this filtering for the whole sample as well as for both populations.

    $ANGSD/angsd -b bam.pop1.filelist -anc chimpHg19.fa -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 3 -out test.pop1.angsd -P 5 -setMinDepth 10 -setMaxDepth 50 -r 1: -GL 1 -doSaf 1

    $ANGSD/angsd -b bam.pop2.filelist -anc chimpHg19.fa -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 3 -out test.pop2.angsd -P 5 -setMinDepth 10 -setMaxDepth 50 -r 1: -GL 1 -doSaf 1

As a side note, to convert files from latest versions of ANGSD (>0.800) to older versions, you can to convert .saf files using `ANGSD realSFS print` tool, as `$ANGSD/misc/realSFS print pop1.saf.idx -oldout 1 > pop1.saf`. Please note that ngsTools has anyway full compatibility only with ANGSD <0.800.

In case you analyse more than one population, you first need to get the subset of overlapping sites, stored in the file `intersect.txt`.

    $ANGSD/misc/realSFS print test.pop1.angsd.saf.idx test.pop2.angsd.saf.idx | cut -f 1-2 > intersect.txt
    $ANGSD/angsd sites index intersect.txt
    N_SITES=`wc -l intersect.txt | cut -f 1 -d " "`


Population genetic differentiation - FST
---------------

We first need to compute sample allele frequencies likelihoods (stored in .saf files) using ANGSD, only for the overlapping filtered sites.

    $ANGSD/angsd -b bam.pop1.filelist -anc chimpHg19.fa -out test.pop1 -P 5 -r 1: -GL 1 -doSaf 1 -sites intersect.txt
    $ANGSD/angsd -b bam.pop2.filelist -anc chimpHg19.fa -out test.pop2 -P 5 -r 1: -GL 1 -doSaf 1 -sites intersect.txt

We then estimate the 2D-SFS to be used as prior.

    N_SITES=`wc -l intersect.txt | cut -f 1 -d " "`
    zcat test.pop1.saf.gz > test.pop1.saf
    zcat test.pop2.saf.gz > test.pop2.saf
    $NGSTOOLS/ngsPopGen/ngs2dSFS -postfiles test.pop1.saf test.pop2.saf -outfile test.pops.2dsfs -nind 5 5 -nsites $N_SITES

Please note that ANGSD can also estimate it (using a Maximum Likelihood approach).

    $ANGSD/misc/realSFS test.pop1.saf.idx test.pop2.saf.idx > test.pops.realSFS.2dsfs
    Rscript $NGSTOOLS/scripts/convertSFS.R test.pops.realSFS.2dsfs > test.pops.realSFS.2dsfs

We can now calculate per-site FST values.

    $NGSTOOLS/ngsPopGen/ngsFST -postfiles test.pop1.saf test.pop2.saf -priorfile test.pops.2dsfs -nind 5 5 -nsites $N_SITES -outfile test.pops.fst

These example scripts will produce a plot and text file with sliding windows values.

    Rscript $NGSTOOLS/scripts/plotFST.R -i test.pops.fst -o test.pops.fst -p intersect.txt -w 1 -s 1
    less -S test.pops.fst.txt
    evince test.pops.fst.pdf


Nucleotide diversity
----------------------------

We use the sample allele frequency probabilities (.saf files) calculated on the previous step and estimate the marginal SFS to be used as priors.

    $ANGSD/misc/realSFS -P 5 test.pop1.saf.idx > test.pop1.sfs
    $ANGSD/misc/realSFS -P 5 test.pop2.saf.idx > test.pop2.sfs

From these priors, we calculate the sample allele frequency posterior probabilities for each population.

    $ANGSD/angsd -bam bam.pop1.filelist -out test.pop1 -doSaf 1 -pest test.pop1.sfs -anc chimpHg19.fa -GL 1 -P 5 -r 1: -sites intersect.txt
    $ANGSD/angsd -bam bam.pop2.filelist -out test.pop2 -doSaf 1 -pest test.pop2.sfs -anc chimpHg19.fa -GL 1 -P 5 -r 1: -sites intersect.txt

We can now calculate some summary statistics, namely number of segregating sites, expected heterozygosity, number of fixed differences and dxy (the latter has been shown to be over-estimated and should be used only for inspecting the distribution and not to make inferences based on its absolute values).

    $NGSTOOLS/ngsPopGen/ngsStat -npop 2 -postfiles test.pop1.saf test.pop2.saf -nsites $N_SITES -nind 5 5 -outfile test.stat

These example scripts will produce a plot and text file with sliding windows values, for 2 populations and only 1 population.

    Rscript $NGSTOOLS/scripts/plotSS.R -i test.stat -p intersect.txt -o test.stat.pdf -n pop1-pop2 -w 5000 -s 1000
    Rscript $NGSTOOLS/scripts/plotSS.R -i test.stat -p intersect.txt -o test.stat.pop1.pdf -n pop1 -w 5000 -s 1000


Population structure - PCA
-----------------------------

We first calculate genotype posterior probabilities, assuming HWE, using ANGSD.

    $ANGSD/angsd -b bam.filelist -nInd 10 -doGeno 32 -doPost 1 -out test.pops -P 5 -r 1: -sites intersect.txt -GL 1 -doMajorMinor 1 -doMaf 1
    zcat test.pops.geno.gz > test.pops.geno

We can now estimate the covariance matrix.

    $NGSTOOLS/ngsPopGen/ngsCovar -probfile test.pops.geno -outfile test.covar -nind 10 -nsites $N_SITES -call 0

For plotting purposes, we create a dummy PLINK cluster file.

    Rscript -e 'write.table(cbind(seq(1,10),rep(1,10),c(rep("A",3),rep("B",3),rep("C",4))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="test.pops.clst", quote=F)'

This script will calculate principal components and plot them.

    Rscript $NGSTOOLS/scripts/plotPCA.R -i test.covar -c 1-2 -a test.pops.clst -o test.pca.pdf
    evince test.pca.pdf

Genetic distances
---------------

We first calculate genotype posterior probabilities, assuming HWE, using ANGSD. Please note that ngsDist is compatible with the most recent versions of ANGSD unlike ngsTools.

    $ANGSD/angsd -b bam.filelist -anc chimpHg19.fa -sites intersect.txt -r 1: -out test.pops -doMajorMinor 1 -doPost 1 -doMaf 1 -doGeno 8 -GL 1 -minMaf 0.05
    N_SITES=$((`zcat test.pops.mafs.gz | wc -l`-1))

We can now estimate pairwise genetic distances.
Here we perform 100 bootstraps, each one with size of 5 sites (these are variable sites since we filtered out non-variable sites).

    cut -d " " -f 1 test.pops.clst > pops.label
    $NGSDIST/ngsDist -verbose 0 -geno test.pops.geno.gz -probs -n_ind 10 -n_sites $N_SITES -labels pops.label -out test.pops.dist -n_threads 10 -n_boot_rep 100 -boot_block_size 5

In case you are interested in producing trees out of these distances, you can used [FastME](http://www.atgc-montpellier.fr/fastme/).

    perl -p -i -e 's/\t/ /g' test.pops.dist
    $FASTME/fastme -D 101 -i test.pops.dist -o test.pops.tree

These scripts will produce a plot of the estimated tree.

    Rscript $NGSTOOLS/scripts/plotTree.R test.pops.tree
    evince test.pops.tree.pdf
