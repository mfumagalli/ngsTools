
options(scipen=20)

fin=commandArgs(T)

sfs=exp(read.table(fin, head=F, stringsAsFactors=F))

sfs[sfs<1e-10]=0
sfs=sfs/sum(sfs)

write.table(sfs, sep="\t", row.names=F, col.names=F)



