
# Usage: Rscript plotFST.R -i infile.fst -o outfile -p positions.txt -w win -s step

library(methods)
library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file'),
make_option(c('-p','--pos_file'), action='store', type='character', default=NULL, help='Input position file'),
make_option(c('-w','--window'), action='store', type='character', default=1, help='Window length'),
make_option(c('-s','--step'), action='store', type='character', default=1, help='Step size'),
make_option(c('-t','--th'), action='store', type='character', default=0, help='Minimum probability of being variable'),
make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
)

opt <- parse_args(OptionParser(option_list = option_list))

# Read input file
values <- read.table(opt$in_file, stringsAsFact=F);
ind <- which(values[,5]>=as.numeric(opt$th))
pos <- read.table(opt$pos_file, stringsAsFact=F, head=F)
cat("There are", nrow(values), "values.\n");
if (nrow(pos)!=nrow(values)) stop("Dimensions of fst values and positions must match. Terminate.\n");
values <- values[ind,]
pos <- pos[ind,]
cat("After removing non-variable (based on your -th value) sites, now there are",nrow(values),"values.\n")
cat("Overall FST:",sum(values[,1])/sum(values[,2]),"\n");

# how many chroms?
uc <- sort(unique(pos$V1))
cpos <- cbind(pos, V3=rep(0, nrow(pos)))
rm(pos)

# set cumulative positions for plotting purpose
cpos$V3[which(cpos$V1==uc[1])]=cpos$V2[which(cpos$V1==uc[1])]
if (length(uc)>1) {
	for (i in 2:length(uc)){
		offset <- max(as.numeric(as.character(cpos$V3)), na.rm=T)+1e5
		cpos$V3[which(cpos$V1==uc[i])] <- cpos$V2[which(cpos$V1==uc[i])]+offset
	}
}

# plot
win <- as.numeric(opt$window);
step <- as.numeric(opt$step);

# if no window scan, assign each value
values$V4[values$V4<0] <- 0
values$V4[values$V4>1] <- NA
if (opt$window==1) {

	df <- data.frame(cbind(Chrom=cpos[,1], Pos=cpos[,2], cum_Pos=cpos[,3], Value=values[,4]));
	df$Chrom <- factor(df$Chrom, levels=uc)
	df[,2:4] <- sapply(df[,2:4], as.character)
	df[,2:4] <- sapply(df[,2:4], as.numeric)


} else {
	fst <- wpos <- wwpos <- c()
	# Windows
	for (i in uc) {
		indi <- which(cpos[,1]==i)
		start <- seq(min(cpos[indi,3]), max(cpos[indi,3]), step);
		end <- start+win-1;
		cchrom <- rep(i, length(start))
		wpos <- c(wpos,round(start+(win/2))); # position of the window in the plot (center)
		for (j in 1:length(start)) {
			ipos <- which(cpos[,3]>=start[j] & cpos[,3]<=end[j])
			fst <- c(fst,sum(values[ipos,1])/sum(values[ipos,2]))
		}
		start <- seq(min(cpos[indi,2]), max(cpos[indi,2]), step);
                end <- start+win-1;
                wwpos <- c(wwpos,round(start+(win/2)));
	}

	fst[which(fst<0)] <- 0

	df <- data.frame(cbind(Chrom=cchrom, Pos=wwpos, cum_Pos=wpos, Value=fst));
        df$Chrom <- factor(df$Chrom, levels=uc)
        df[,2:4] <- sapply(df[,2:4], as.character)
        df[,2:4] <- sapply(df[,2:4], as.numeric)

}

# Data
write.table(df, file=paste(opt$out_file,".txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)
# Plot

qplot(data=df, x=cum_Pos, y=Value, col=Chrom)
ggsave(paste(opt$out_file,".pdf",sep="",collapse=""))
unlink("Rplots.pdf", force=TRUE)


