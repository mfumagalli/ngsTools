
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
    FASTME=/data/data/Software/fastme-2.1.4/binaries

Second, create all directories where you will be working:

    mkdir Tutorial
    cd Tutorial
    mkdir Data
    mkdir Results

Data
----------

As an illustration, we will use 80 BAM files of human samples (of African, European, East Asian, and Native American descent), a reference genome, and a putative ancestral sequence.
BAM files have been downsampled to a mean depth of around 2X.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).
All data is publicly available.

A pipeline to retrieve such data is provided [here](https://github.com/mfumagalli/Weggis/blob/master/Files/data.sh).
You need to have 'samtools' (tested with version 1.3.1), 'gunzip' and 'wget' installed in your /usr/bin to run this.
Otherwise edit the first line of 'data.sh' file to set the appropriate paths.

    cp $NGSTOOLS/Files/* .
    bash $NGSTOOLS/Scripts/data.sh phase3_bamlist.txt

As a note for the general use, in case an ancestral sequence is not available, analyses on FST, PCA, nucleotide diversity (but not the number of fixed differences) can be carried out using the reference sequence to polarise your data. 
Please be aware that, under this scenario, some quantities (e.g. the unfolded joint site frequency spectrum) will be nonsense.




Filtering using ANGSD
----------------------

Check the distribution of depth and quality scores.

    $ANGSD/angsd -b bam.filelist -doQsDist 1 -doCounts 1 -maxDepth 100 -doDepth 1 -out bam.qc

    Rscript $NGSTOOLS/scripts/plotQC.R bam.qc

    less bam.qc.info
    evince bam.qc.pdf

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
