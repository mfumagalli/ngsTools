
## Script provided by Dean Ousby and adapted

library(plot3D)

args=commandArgs(T)
fin=args[1]
nPop1=as.numeric(args[2])*2
nPop2=as.numeric(args[3])*2
rm(args)

ANGSD.2D.SFS <- scan(paste(fin, sep=""), quiet=T)
ANGSD.2D.SFS <- t(matrix(ANGSD.2D.SFS, nrow=nPop1+1, ncol=nPop2+1))

# Since we used 1 MB region in our examples, the non-polymorphic sites (0,0) will greatly outweigh the variant sites, therefore if you don't silence these out they will mask the rest of the SFS. Try with and without this step and you'll see what I mean

# mask non-variant sites
ANGSD.2D.SFS[1,1] <- 0
ANGSD.2D.SFS[nrow(ANGSD.2D.SFS),ncol(ANGSD.2D.SFS)] <- 0

fout=paste(fin,".pdf", sep="", collapse="")

pdf(file=fout)

hist3D(x = seq(0,1,length.out = nrow(ANGSD.2D.SFS)), y = seq(0,1,length.out=ncol(ANGSD.2D.SFS)), ANGSD.2D.SFS, cex.lab=1.2, xlab="LWK",ylab="TSI",zlab="SFS", main=paste("2D-SFS"),pin=c(10,0), cex.main=1.4, zlim=c(0,max(ANGSD.2D.SFS)))

dev.off()

cat("Output file:", fout, "\n")


