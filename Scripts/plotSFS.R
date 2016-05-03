
# from ANGSD website, adapted
# this is written knowning that the path is Results/

fins=commandArgs(T)

pops=unlist(strsplit(unlist(strsplit(fins, split="Results/")), split=".sfs"))
cat("Populations:", pops, "\n")
fout=paste("Results/", paste(pops, sep="",collapse="_"), ".pdf", sep="", collapse="")

#function to normalize
norm <- function(x) x/sum(x)

pdf(file=fout)

msfs=pvars=c()
for (i in 1:length(fins)) {
	#read data
	sfs <- (scan(fins[i], quiet=T))
	cat(length(sfs)," ")
	#the variability as percentile
	pvar <- (1-norm(sfs)[1]-norm(sfs)[length(sfs)])*100
	pvars[i]= paste(pops[i], "; Variability=",round(pvar,3),"%",sep="",collapse="")
	#the variable categories of the sfs
	sfs <- norm(sfs[-c(1,length(sfs))]) 
	msfs <- rbind(msfs, sfs)
}

barplot(msfs, beside=T, legend=pvars, xlab="Chromosomes", names=1:length(sfs), ylab="Proportions", main="SFS")

dev.off()

cat("Output file:",fout,"\n")


