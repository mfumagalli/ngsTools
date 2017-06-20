
# from ANGSD website, adapted

args <- commandArgs(T)
fins <- unlist(strsplit(args[1], split="-"))
pops <- unlist(strsplit(args[2], split="-"))
fold <- as.numeric(args[3])
fout <- args[4]
rm(args)

cat("Populations:", pops, "\n")
cat("Maximum frequency: ")

#function to normalize
norm <- function(x) x/sum(x)

pdf(file=fout)

pvars <- rep(NA, length(fins))

maxAF <- 0
for (i in 1:length(fins)) {
	sfs <- scan(fins[i], quiet=T)
	if (fold==0) {
		maxAF <- max(maxAF, length(sfs))
	} else {
		maxAF <- max(maxAF, ceiling(length(sfs)/2))
	}
}

ncols <- maxAF-2
if (fold) ncols <- maxAF-1
msfs <- matrix(NA, nrow=length(pops), ncol=ncols)

foldSFS <- function(sfs) {
	nInd <- (length(sfs)-1)/2
	folded <- sfs[1:(nInd+1)]
	for (i in 1:(length(folded)-1)) folded[i] <- sfs[i] + sfs[length(sfs)-i+1]
	folded
}

for (i in 1:length(fins)) {

	sfs <- scan(fins[i], quiet=T)
	if (fold) sfs <- foldSFS(sfs)
	cat(length(sfs)-1," ")

	#the variability as percentile
	pvar <- (1-norm(sfs)[1]-norm(sfs)[length(sfs)])*100
	if (fold) pvar <- (1-norm(sfs)[1])*100
	pvars[i] <-  paste(pops[i], " (",length(sfs)-1, "); Var.=",round(pvar,3),"%",sep="",collapse="")

	#the variable categories of the sfs
	if (fold) {
		msfs[i,1:(length(sfs)-1) ] <- norm(sfs[-c(1)])
	} else {
		msfs[i,1:(length(sfs)-2) ] <- norm(sfs[-c(1,length(sfs))]) 
	}
}

xlab <- "Derived allele frequency"
if (fold) xlab <- "Minor allele frequency"
barplot(msfs, beside=T, legend=pvars, xlab=xlab, names=1:ncol(msfs), ylab="Proportions", main="SFS")

invisible(dev.off())

cat("\n")

