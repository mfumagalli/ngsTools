
options(scipen=20)

args <- commandArgs(T)
fin <- args[1]
n1 <- as.numeric(args[2])
n2 <- as.numeric(args[3])
rm(args)

# Read input SFS
sfs <- read.table(fin, head=F, stringsAsFactors=F)
# Normalize SFS
sfs <- sfs/sum(sfs)
# Format output
sfs <- format(round(sfs, 6), nsmall=6)
# Print formated SFS
write.table(matrix(sfs,nrow=((2*n1)+1),ncol=((2*n1)+1), sep="\t", row.names=F, col.names=F)

