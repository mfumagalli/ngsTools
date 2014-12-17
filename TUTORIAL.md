
A short tutorial to some basic analyses using ngsTools from BAM files
===============

Settings
----------

Set directories to installed programs:

	ANGSD=~/angsd
	SAMTOOLS=samtools
	NGSTOOLS=~/ngsTools
	NGSDIST=~/ngsDist

Data
----------

Dowload example datasets 

	wget http://popgen.dk/software/download/angsd/bams.tar.gz
	tar xf bams.tar.gz

	for i in bams/*.bam; do $SAMTOOLS index $i; done

	ls bams/*.bam > bam.filelist

	ls bams/*.bam | head -n 5 > bam.pop1.filelist
	ls bams/*.bam | tail -n 5 > bam.pop2.filelist

	wget wget http://dna.ku.dk/~thorfinn/hg19ancNoChr.fa.gz
	zcat hg19ancNoChr.fa.gz > chimpHg19.fa

	$SAMTOOLS faidx chimpHg19.fa



Filtering
----------------------

	$ANGSD/angsd -b bam.filelist -doQsDist 1 -doCounts 1 -maxDepth 100 -doDepth 1 -out bam.qc

	Rscript plotQC.R bam.qc

	less bam.qc.info
	evince bam.qc.pdf

	$ANGSD/angsd -b bam.filelist -anc chimpHg19.fa -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 6 -out test.pops.angsd -P 5 -setMinDepth 20 -setMaxDepth 100 -r 1: -GL 1 -doSaf 1 -doMaf 2 -minMaf 0.05 -doMajorMinor 1

        $ANGSD/angsd -b bam.pop1.filelist -anc chimpHg19.fa -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 3 -out test.pop1.angsd -P 5 -setMinDepth 10 -setMaxDepth 50 -r 1: -GL 1 -doSaf 1

	$ANGSD/angsd -b bam.pop2.filelist -anc chimpHg19.fa -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 3 -out test.pop2.angsd -P 5 -setMinDepth 10 -setMaxDepth 50 -r 1: -GL 1 -doSaf 1

	gunzip -c test.pop1.angsd.saf.pos.gz test.pop2.angsd.saf.pos.gz | sort -S 50% | uniq -d | sort -k1,1 -S 50% | gzip > tmp.gz
	gunzip -c tmp.gz test.pops.angsd.saf.pos.gz | sort -S 50% | uniq -d | sort -k1,1 -S 50% > intersect.txt
        rm tmp* test.pop?.angsd.*

	$ANGSD/angsd sites index intersect.txt
	N_SITES=`wc -l intersect.txt | cut -f 1 -d " "`


Population genetic differentiation - FST
---------------

	$ANGSD/angsd -b bam.pop1.filelist -anc chimpHg19.fa -out test.pop1 -P 5 -r 1: -GL 1 -doSaf 1 -sites intersect.txt
	$ANGSD/angsd -b bam.pop2.filelist -anc chimpHg19.fa -out test.pop2 -P 5 -r 1: -GL 1 -doSaf 1 -sites intersect.txt

	#$ANGSD/misc/realSFS 2dsfs test.pop1.saf test.pop2.saf 10 10 > test.pops.angsd.log.2dsfs
	#Rscript convertSFS.R test.pops.log.2dsfs > test.pops.angsd.2dsfs

	N_SITES=`wc -l intersect.txt | cut -f 1 -d " "`
	$NGSTOOLS/ngsPopGen/ngs2dSFS -postfiles test.pop1.saf test.pop2.saf -outfile test.pops.2dsfs -nind 5 5 -nsites $N_SITES

	$NGSTOOLS/ngsPopGen/ngsFST -postfiles test.pop1.saf test.pop2.saf -priorfile test.pops.2dsfs -nind 5 5 -nsites $N_SITES -outfile test.pops.fst

	Rscript plotFST.R -i test.pops.fst -o test.pops.fst -p intersect.txt -w 1 -s 1
	less -S test.pops.fst.txt
	evince test.pops.fst.pdf


Nucleotide diversity
----------------------------

	$ANGSD/angsd -b bam.pop1.filelist -anc chimpHg19.fa -out test.pop1 -P 5 -r 1: -GL 1 -doSaf 1 -sites intersect.txt
        $ANGSD/angsd -b bam.pop2.filelist -anc chimpHg19.fa -out test.pop2 -P 5 -r 1: -GL 1 -doSaf 1 -sites intersect.txt

	$ANGSD/misc/realSFS test.pop1.saf 10 -P 5 > test.pop1.sfs
	$ANGSD/misc/realSFS test.pop2.saf 10 -P 5 > test.pop2.sfs

	$ANGSD/angsd -bam bam.pop1.filelist -out test.pop1 -doSaf 1 -pest test.pop1.sfs -anc chimpHg19.fa -GL 1 -P 5 -r 1: -sites intersect.txt
	$ANGSD/angsd -bam bam.pop2.filelist -out test.pop2 -doSaf 1 -pest test.pop2.sfs -anc chimpHg19.fa -GL 1 -P 5 -r 1: -sites intersect.txt

	N_SITES=`wc -l intersect.txt | cut -f 1 -d " "`

	$NGSTOOLS/ngsPopGen/ngsStat -npop 2 -postfiles test.pop1.saf test.pop2.saf -nsites $N_SITES -nind 5 5 -outfile test.stat

	Rscript plotSS.R -i test.stat -p intersect.txt -o test.stat.pdf -n pop1-pop2 -w 5000 -s 1000
	Rscript plotSS.R -i test.stat -p intersect.txt -o test.stat.pop1.pdf -n pop1 -w 5000 -s 1000


Population structure - PCA
-----------------------------

	$ANGSD/angsd -b bam.filelist -nInd 10 -doGeno 32 -doPost 1 -out test.pops -P 5 -r 1: -sites intersect.txt -GL 1 -doMajorMinor 1 -doMaf 2
	zcat test.pops.geno.gz > test.pops.geno
	N_SITES=`wc -l intersect.txt | cut -f 1 -d " "`

	$NGSTOOLS/ngsPopGen/ngsCovar -probfile test.pops.geno -outfile test.covar -nind 10 -nsites $N_SITES -call 0

	Rscript -e 'write.table(cbind(seq(1,10),rep(1,10),c(rep("A",3),rep("B",3),rep("C",4))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="test.pops.clst", quote=F)'

	Rscript plotPCA.R -i test.cover -c 1-2 -a test.pops.clst -o test.pca.pdf
	evince test.pca.pdf


Genetic distances
---------------


	$ANGSD/angsd -b bam.filelist -anc chimpHg19.fa -sites intersect.txt -r 1: -out test.pops -doMajorMinor 1 -doPost 1 -doMaf 1 -doGeno 8 -GL 1 -minMaf 0.05

	N_LINES=`zcat test.pops.mafs.gz | wc -l`
	N_SITES=`expr $NS - 1`

        $NGSDIST/ngsDist -verbose 0 -geno test.pops.geno.gz -probs -n_ind 10 -n_sites $N_SITES -labels pops.label -out_prefix test.pops -n_threads 10 -n_boot_rep 100 -boot_block_size 5

	$NGSDIST/FastME/fastme_linux64 -d 101 -i test.pops.dist -o test.pops.tree -m b -n b

	Rscript plotTree.R test.pops.tree
	evince test.pops.tree.pdf




