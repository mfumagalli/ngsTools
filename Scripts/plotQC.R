
# quick script to compute percentiles and plot distributions of quality scores and depths

fin=commandArgs(T)

cat("", file=paste(fin,".info",sep="",collapse=""))

pdf(paste(fin,".pdf",sep="",collapse=""))
par(mfrow=c(3,1))
## barplot q-scores
qs=read.table(paste(fin,".qs",sep="",collapse=""), head=T, stringsAsFactors=F)
barplot(height=qs$counts, names.arg=qs$qscore, xlab="Q-score", ylab="Counts")
qs=cbind(qs, perc=cumsum(qs$counts/1e4)/sum(qs$counts/1e4,na.rm=T))
write.table( qs, row.names=F, col.names=T, quote=F, sep="\t", file=paste(fin,".info",sep="",collapse=""), append=T)

## global depth
dep=as.numeric(scan(paste(fin,".depthGlobal", sep="",collapse=""),what="char", quiet=T))
barplot(height=dep, names.arg=seq(1,length(dep))-1, xlab="Global Depth", ylab="Counts")
cat("\nGlobal_depth\tpercentile\n", file=paste(fin,".info",sep="",collapse=""), append=T)
write.table( cbind(seq(1,length(dep))-1,cumsum(dep)/sum(dep)), row.names=F, col.names=F, quote=F, sep="\t", file=paste(fin,".info",sep="",collapse=""), append=T)

## sample depth
deps=read.table(paste(fin,".depthSample", sep="",collapse=""),head=F, stringsAsFactors=F)
## per sample
depp=matrix(NA, nrow=nrow(deps), ncol=ncol(deps))
for (i in 1:nrow(depp)) {
depp[i,]=apply(X=deps[i,], FUN=sum, MAR=2)
}

# this will plot only the first 'xl' bins of depth
xl=10
barplot(depp[,1:xl], names.arg=c(seq(1,ncol(deps))-1)[1:xl], xlab="Sample Depth", ylab="Counts", beside=T)

dev.off()


