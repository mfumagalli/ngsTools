
# Usage: Rscript -i infile.stat -o outfile -p positions.txt -n name1-name2 -w window_size -s step_size
# Usage: Rscript -i infile.stat -o outfile -p positions.txt -n name1 -w window_size -s step_size # if only 1 pop

# Output:
# outfile.pdf & output.txt: plot and text with sliding windows values

library(methods)
library(grid)
library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file'),
                    make_option(c('-n','--names'), action='store', type='character', default=1-2, help='Name(s) of population(s)'),
                    make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file'),
			make_option(c('-p','--pos_file'), action='store', type='character', default=NULL, help='Input position file'),
			make_option(c('-w','--window'), action='store', type='character', default=1, help='Window length'),
			make_option(c('-s','--step'), action='store', type='character', default=1, help='Step size')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

# from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# How many pops
pops <- as.character(strsplit(opt$names, "-", fixed=TRUE)[[1]]);
npop=length(pops);
cat("Detected", npop, "population(s).\n")

# Read input file
values <- read.table(opt$in_file, stringsAsFact=F, head=F);
pos <- read.table(opt$pos_file, stringsAsFact=F, head=F)

if (nrow(pos)!=nrow(values)) stop("Dimensions of values and positions must match. Terminate.\n");

# how many chroms?
uc=sort(unique(pos$V1))
pos=cbind(pos, V3=rep(0, nrow(pos)))

# set cumulative positions for plotting purpose
pos$V3[which(pos$V1==uc[1])]=pos$V2[which(pos$V1==uc[1])]
if (length(uc)>1) {
        for (i in 2:length(uc)){
                offset=max(as.numeric(as.character(pos$V3)), na.rm=T)+1e5
                pos$V3[which(pos$V1==uc[i])]=pos$V2[which(pos$V1==uc[i])]+offset
        }
}


win=as.numeric(opt$window);
step=as.numeric(opt$step);

avalues=achr=apos=c()
      
for (i in uc) {

	cpos=pos[which(pos$V1==i),]
	# create dataframe of windows values here
	# Windows
	win=as.numeric(opt$window);
	step=as.numeric(opt$step);
	start=seq(min(cpos$V3), max(cpos$V3), step);
	end=start+win-1;
	wpos=round(start+(win/2)); # position of the window in the plot (center)

	# windows values
	sub=matrix(NA, nrow=length(start), ncol=(ncol(values)-2)); # S, pi (S2, pi2, diff, dxy)
	for (j in 1:length(start)) {
		ipos=which(cpos$V3>=start[j] & cpos$V3<=end[j])
		sub[j,]=apply(MAR=2, X=values[ipos,3:(ncol(values))], FUN=sum, na.rm=T)
	}

	avalues=rbind(avalues, sub)
	apos=c(apos, wpos)
	achr=c(achr, rep(i,length(wpos)))
}


# Plot
if (npop==2) {
  title <- "";

  # Data
  df = data.frame(cbind( Pop=c(rep(pops[1],length(apos)),rep(pops[2],length(apos))), Pos=rep(apos,2), Segr.sites=c(avalues[,1], avalues[,3]), Exp.heterozygosity=c(avalues[,2],avalues[,4]), Fixed.differences=(rep(avalues[,5],2)), Dxy=(rep(avalues[,6],2)),  Pops=c(rep(paste(pops[1],"\n",pops[2]),length(apos)*2)) ) );
  df[,2:6] = sapply(df[,2:6], as.character)
  df[,2:6] = sapply(df[,2:6], as.numeric)

  p1 = ggplot(data=df, aes(x=Pos, y=Segr.sites, color=Pop)) + geom_line() + ggtitle(title)
  p2 = ggplot(data=df, aes(x=Pos, y=Exp.heterozygosity, color=Pop)) + geom_line() + ggtitle(title)
  p3 = ggplot(data=df, aes(x=Pos, y=Fixed.differences, color=Pops)) + geom_line() + ggtitle(title)
  p4 = ggplot(data=df, aes(x=Pos, y=Dxy, color=Pops)) + geom_line() + ggtitle(title)

  pdf(paste(opt$out_file,".pdf",sep="",collapse=""));
  multiplot(p1, p2, p3, p4, cols=1)
  null <- dev.off();

  df=df[,1:6]

}

if (npop==1) {
  df = data.frame(cbind( Pop=rep(pops[1],length(apos)), Pos=apos, Segr.sites=avalues[,1], Exp.heterozygosity=avalues[,2]));
  df[,2:4] = sapply(df[,2:4], as.character)
  df[,2:4] = sapply(df[,2:4], as.numeric)

  title <- "";

  p1 = ggplot(data=df, aes(x=Pos, y=Segr.sites, color=Pop)) + geom_line() + ggtitle(title)
  p2 = ggplot(data=df, aes(x=Pos, y=Exp.heterozygosity, color=Pop)) + geom_line() + ggtitle(title)

  pdf(paste(opt$out_file,".pdf",sep="",collapse=""));
  multiplot(p1, p2, cols=1)
  null <- dev.off();

  df=df[,1:4]
}

write.table(df, file=paste(opt$out_file,".txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)

