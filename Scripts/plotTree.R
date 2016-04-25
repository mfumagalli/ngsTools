
fin=commandArgs(T)

library(phangorn)

trees=read.tree(fin, skip=2)
tree=read.tree(fin)[[1]]

pdf(file=paste(fin,".pdf",sep="",collaps=""))
plotBS(tree, trees, type="phylo")
dev.off()



