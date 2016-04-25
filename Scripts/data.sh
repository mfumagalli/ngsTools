
## Pipeline to download and process the data to be used for this tutorial

# file with BAM names $1 is given as input

# set path if not in your /usr/bin
SAMTOOLS=samtools
WGET=wget
GUNZIP=gunzip
echo Is this your path to samtools? $SAMTOOLS
echo Is this your path to wget? $WGET
echo Is this your path to gunzip? $GUNZIP

echo Retrieving file names...
NS=20 # first 20 samples
for POP in LWK TSI CHB PEL;
do

	INPUT=$POP.txt
	INDS=`cat $INPUT`
	OUTPUT=$POP.BAMs.txt
	> $OUTPUT
	> tmp

	for i in $INDS;
	do
		grep $i $1 >> tmp
	done
	head -n $NS tmp > $OUTPUT
	rm tmp
done

echo Downloading BAM files...
for POP in LWK TSI CHB PEL;
do
	mkdir Data/$POP.BAMs
	echo $POP
	INDLIST=`cat $POP.BAMs.txt`
	for i in $INDLIST;
	do
		NAME=`echo -n $i | tail -c 58`
		echo $NAME
		$SAMTOOLS view -s 0.25 -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/$i 11:61000000-62000000 > Data/$POP.BAMs$NAME 2> /dev/null
	done
done
rm *.bai
# this has created files and folders in Data/PEL.BAMs/* and TSI and LWK and CHB

# create file with list of BAMs
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/CHB.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL.bamlist
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL_noCHB.bamlist
ls Data/LWK.BAMs/*.bam > LWK.bamlist
ls Data/TSI.BAMs/*.bam > TSI.bamlist
ls Data/CHB.BAMs/*.bam > CHB.bamlist
ls Data/PEL.BAMs/*.bam > PEL.bamlist

# download ancestral sequence
echo Downloading and processing ancestral sequence...
$WGET http://dna.ku.dk/~thorfinn/hg19ancNoChr.fa.gz &>/dev/null
$GUNZIP hg19ancNoChr.fa.gz &>/dev/null
$SAMTOOLS faidx hg19ancNoChr.fa
mv hg19ancNoChr.* Data/.

# download reference sequence
echo Downloading and processing reference sequence...
$WGET http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz &>/dev/null
$GUNZIP hs37d5.fa.gz &>/dev/null
$SAMTOOLS faidx hs37d5.fa
mv hs37d5.* Data/.

echo Done!
ls -lh Data/* > Data/download.log
echo Open Data/download.log to see what files have been generated.

exit



