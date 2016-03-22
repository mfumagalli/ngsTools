
A short tutorial for some basic analyses using ngsTools
===============

Please be sure you are using the most updated version of ngsTools. In doubt please run: 

	git pull
	make clean
	make
	
inside the ngsTools directory.


Settings
----------

Set directories to installed programs:

        NGSTOOLS=~/ngsTools
	ANGSD=$NGSTOOLS/angsd
	NGSDIST=$NGSTOOLS/ngsDist
	SAMTOOLS=samtools
	FASTME=~/fastme-2.1.4/binaries

Data
----------

Dowload example datasets (taken from ANGSD website) of 10 individuals.
Index BAM and fasta files.
Create a subset of 5 individuals into 2 populations.

	wget http://popgen.dk/software/download/angsd/bams.tar.gz
	tar -xvf bams.tar.gz

	for i in bams/*.bam; do $SAMTOOLS index $i; done

	ls bams/*.bam > bam.filelist

	ls bams/*.bam | head -n 5 > bam.pop1.filelist
	ls bams/*.bam | tail -n 5 > bam.pop2.filelist

	wget wget http://dna.ku.dk/~thorfinn/hg19ancNoChr.fa.gz
	zcat hg19ancNoChr.fa.gz > chimpHg19.fa

	$SAMTOOLS faidx chimpHg19.fa

As a note for the general use, in case an ancestral sequence is not available, analyses on FST, PCA, nucleotide diversity (but not number of fixed differences) can be carried out using the reference sequence to polarise your data. Please be aware that, under this scenario, some quantities (e.g. the unfolded joint site frequency spectrum) will be nonsense.


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

Please note that ANGSD can also estimate it (using a Maximum Likelihood approach). You then need to convert its output to be read by ngsTools.

        $ANGSD/misc/realSFS test.pop1.saf.idx test.pop2.saf.idx > test.pops.realSFS.2dsfs
        #Rscript $NGSTOOLS/scripts/convertSFS.R test.pops.realSFS.2dsfs > test.pops.realSFS.2dsfs

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

We first calculate genotype posterior probabilities, assuming HWE, using ANGSD.

	$ANGSD/angsd -b bam.filelist -anc chimpHg19.fa -sites intersect.txt -r 1: -out test.pops -doMajorMinor 1 -doPost 1 -doMaf 1 -doGeno 8 -GL 1 -minMaf 0.05
	N_SITES=$((`zcat test.pops.mafs.gz | wc -l`-1))

We can now estimate pairwise genetic distances.
Here we perform 100 bootstraps, each one with size of 5 sites (these are variable sites since we filtered out non-variable sites).

        cut -d " " -f 1 test.pops.clst > pops.label
        $NGSDIST/ngsDist -verbose 0 -geno test.pops.geno.gz -probs -n_ind 10 -n_sites $N_SITES -labels pops.label -out test.pops.dist -n_threads 10 -n_boot_rep 100 -boot_block_size 5

In case you are interested in producing trees out of these distances, you can used FastME.

   	perl -p -i -e 's/\t/ /g' test.pops.dist
	$FASTME/fastme -D 101 -i test.pops.dist -o test.pops.tree

These scripts will produce a plot of the estimated tree.

	Rscript $NGSTOOLS/scripts/plotTree.R test.pops.tree
	evince test.pops.tree.pdf
