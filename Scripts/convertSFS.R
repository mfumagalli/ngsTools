
options(scipen=20)

fin <- commandArgs(T)
# Read input SFS
sfs <- read.table(fin, head=F, stringsAsFactors=F)
# Number of categories
n_categ <- length(sfs)
# Normalize SFS
sfs <- sfs/sum(sfs)
# Format output
sfs <- format(round(sfs, 6), nsmall=6)
# Print formated SFS
write.table(matrix(sfs,nrow=sqrt(n_categ),ncol=sqrt(n_categ)), sep="\t", row.names=F, col.names=F)

