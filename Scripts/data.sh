
## Pipeline to download and process the data to be used for this tutorial

# file with BAM names $1 is given as input

# set path if not in your /usr/bin
SAMTOOLS=`which samtools`
WGET=`which wget`
GUNZIP=`which gunzip`
BGZIP=`which bgzip`
PARALLEL=`which parallel`
# echo Is this your path to samtools? $SAMTOOLS
# echo Is this your path to wget? $WGET
# echo Is this your path to gunzip? $GUNZIP
# echo Is this your path to bgzip? $BGZIP

if [ SAMTOOLS == "" ];
then
	echo Please install samtools and add it to your PATH environment variable
fi

if [ WGET == "" ];
then
	echo Please install wget and add it to your PATH environment variable.
fi

if [ GUNZIP == "" ];
then
	echo Please install gunzip and add it to your PATH environment variable.
fi

if [ BGZIP == "" ];
then
	echo Please install bgzip and add it to your PATH environment variable.
fi

if [ PARALLEL == "" ];
then
	echo Please install GNU parallel and add it to your PATH environment variable.
fi

CHROM=11

echo Retrieving file names...
NS=20 # first 20 samples
for POP in LWK TSI PEL;
do

	INPUT=$POP.txt
	INDS=`cat $INPUT`
	OUTPUT=$POP.BAMs.txt
	> $OUTPUT
	> tmp

	for i in $INDS;
	do
		grep $i phase3_bamlist.txt >> tmp
	done
	head -n $NS tmp > $OUTPUT
	rm tmp
done

echo Downloading BAM files...
for POP in LWK TSI PEL;
do
	mkdir Data/$POP.BAMs
	echo $POP
	INDLIST=`cat $POP.BAMs.txt`
	for i in $INDLIST;
	do
		NAME=`echo -n $i | tail -c 58`
		echo $NAME
		$SAMTOOLS view -bs 0.5 -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/$i $CHROM:61000000-62000000 > Data/$POP.BAMs$NAME 2> /dev/null
	done
done
rm *.bai
# this has created files and folders in Data/PEL.BAMs/* and TSI and LWK

# create file with list of BAMs
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL.bamlist
ls Data/LWK.BAMs/*.bam > LWK.bamlist
ls Data/TSI.BAMs/*.bam > TSI.bamlist
ls Data/PEL.BAMs/*.bam > PEL.bamlist

# index all .bam files
cat ALL.bamlist | parallel "samtools index {}"

# download ancestral sequence
echo Downloading and processing ancestral sequence...
wget ftp://ftp.ensembl.org/pub/release-65/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e65.tar.bz
tar xjf homo_sapiens_ancestor_GRCh37_e65.tar.bz
cp homo_sapiens_ancestor_GRCh37_e65/homo_sapiens_ancestor_$CHROM.fa Data/tmp.fa
sed "1s/.*/>11/" Data/tmp.fa > Data/anc.fa
rm Data/tmp.fa
$BGZIP Data/anc.fa
$SAMTOOLS faidx Data/anc.fa.gz
rm -rf homo_sapiens_ancestor_GRCh37_e65*

# download reference sequence
echo Downloading and processing reference sequence...
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz &> /dev/null
zcat human_g1k_v37.fasta.gz > Data/tmp.fa 2> /dev/null
$SAMTOOLS faidx Data/tmp.fa
$SAMTOOLS faidx Data/tmp.fa $CHROM > Data/ref.fa
$BGZIP Data/ref.fa
$SAMTOOLS faidx Data/ref.fa.gz
rm Data/tmp.fa

echo Done!
ls -lh Data/* > Data/download.log
echo Open Data/download.log to see which files have been generated.

exit



