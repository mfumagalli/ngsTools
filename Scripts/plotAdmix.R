
# this is taken and modified from NGSadmix website

# Usage: Rscript -i infile.qopt -a annotation.file -o outfile.pdf 

library(optparse)

# colorblind friendly: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add: scale_fill_manual(values=cbPalette)
# To use for line and point colors, add: scale_colour_manual(values=cbPalette)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file (output from NGSadmix)'),
			#make_option(c('-a','--annot_file'), action='store', type='character', default=NULL, help='Annotation file with individual classification (2 column TSV with ID and ANNOTATION)'),
			make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

# Admixture
admix <- t(as.matrix(read.table(opt$in_file)))
K <- nrow(admix)
if (K > 8) stop("Maximum number of colours for admixture proportions is 8.")

# Annotation
# Read annot file
#annot <- read.table(opt$annot_file, sep="\t", header=T); # note that plink cluster files are usually tab-separated
#pop <- as.character(annot$CLUSTER)

pdf(file=opt$out_file)

## CHANGE HERE!!!

barplot(admix, col=cbPalette[1:K], space=0, border=NA, xlab="Individuals", ylab="Admixture")

#pop<-read.table("pop.info",as.is=T)
#admix<-t(as.matrix(read.table("myoutfiles.qopt")))
#admix<-admix[,order(pop[,1])]
#pop<-pop[order(pop[,1]),]
#h<-barplot(admix,col=1:3,space=0,border=NA,xlab="Individuals",ylab="admixture")
#text(tapply(1:nrow(pop),pop[,1],mean),-0.05,unique(pop[,1]),xpd=T)

invisible(dev.off())












