
## Script provided by Dean Ousby and adapted

library(plot3D)

args <- commandArgs(T)
fin <- args[1]
pops <- unlist(strsplit(args[2],split="-"))
nPops <- as.numeric(unlist(strsplit(args[3], split="-")))*2
rm(args)

ANGSD.2D.SFS <- scan(paste(fin, sep=""), quiet=T)
ANGSD.2D.SFS <- t(matrix(ANGSD.2D.SFS, nrow=nPops[2]+1, ncol=nPops[1]+1))

# mask non-variant sites
ANGSD.2D.SFS[1,1] <- 0
ANGSD.2D.SFS[nrow(ANGSD.2D.SFS),ncol(ANGSD.2D.SFS)] <- 0

fout <- paste(fin,".pdf", sep="", collapse="")

pdf(file=fout)

hist3D(x = seq(0,1,length.out = nrow(ANGSD.2D.SFS)), y = seq(0,1,length.out=ncol(ANGSD.2D.SFS)), ANGSD.2D.SFS, cex.lab=1.2, xlab=pops[1], ylab=pops[2], zlab="Frequency", main=paste("2D-SFS"),pin=c(10,0), cex.main=1.4, zlim=c(0,max(ANGSD.2D.SFS)))

invisible(dev.off())



