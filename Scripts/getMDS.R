
library(optparse)
library(plyr)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default="stdin", help='Input file'),
                    make_option(c('--no_header'), action='store_true', type='logical', default=FALSE, help='Input file has no header'),
                    make_option(c('--data_symm'), action='store_true', type='logical', default=FALSE, help='Input matrix is symmetric'),
                    make_option(c('-m','--method'), action='store', type='character', default="PCA", help='What method to run: PCA, fastPCA, MDS. [%default]'),
                    make_option(c('-n','--n_comp'), action='store', type='numeric', default=10, help='Number of components to calculate. [%default]'),
                    make_option(c('-r','--s_ref'), action='store', type='numeric', default=NULL, help='Use first "N" samples as ref for PCA calculations (PCA only!).'),
                    make_option(c('-p','--s_proj'), action='store', type='numeric', default=NULL, help='Column labels, from "s_ref", to just project onto the previous calculated PCA components (PCA only)'),
                    make_option(c('-x','--s_excl'), action='store', type='character', default=NULL, help='File or CSV of samples to exclude'),
                    make_option(c('-d','--data'), action='store', type='character', default="general", help='Type of data: geno, general (affects how the variance is calculated). [%default]'),
                    make_option(c('-N','--miss_char'), action='store', type='character', default="NA", help='Missing data character. [%default]'),
                    make_option(c('-o','--out_file'), action='store', type='character', default="", help='Output file')
)
opt <- parse_args(OptionParser(option_list = option_list))

#################################################################################
## Functions written by Martin Sikora
## Get pca from a matrix. Rownames are snpIds, colnames sampleIds
## Optionally build space from subset of samples & project the rest
## If genotypes, the formula of patterson et al is used, which is based
## on the expected variance in allele frequency due to drift (sqrt(p * (1 - p))

## loadings are calculated as in Zou et al. Hum Hered 2010 June; 70(1): 9–22
getPcaGT <- function(gt, idsPca = colnames(gt), pcs = 10, var = "geno"){
    doProjection <- !all(colnames(gt) %in% idsPca)

    ## filter variant & missing snps in samples for building pca space
    nNonMissing <- rowSums(!is.na(gt[, idsPca]))
    idxFixed <- (rowSums(gt[, idsPca], na.rm = TRUE) / (2 * nNonMissing)) %in% c(0,1) ## fixed snps in nonmissing genotypes
    idxMissing <- rowSums(is.na(gt[, idsPca])) == length(idsPca) ## snps with all missing data
    gt <- gt[!idxMissing & !idxFixed,]

    ## center & normalize
    gtMean <- rowMeans(gt[, idsPca], na.rm = TRUE)
    if(var == "general"){
      gtVar <- var(gtMean)
    }else if(var == "geno"){
      gtVar <- (gtMean / 2)*(1 - gtMean / 2)
    }else{
      cat("ERROR: invalid data type!", fill=TRUE)
      quit("no")
    }
    gtStd <- sqrt(gtVar)
    gt1 <- (gt - gtMean) / gtStd
    gt1[is.na(gt1)] <- 0 ## set missing values to have mean gt (0 after normalizing)

    ## make pca
    pca <- svd(gt1[, idsPca])
    rownames(pca$u) <- rownames(gt)
    rownames(pca$v) <- idsPca

    res <- pca$v[, 1:pcs]
    colnames(res) <- paste("PC", 1:pcs, sep = "")

    ## calculate loadings for each SNP
    l <- sapply(1:ncol(pca$u), function(x){
        r <- pca$u[, x] * sqrt(pca$d[x] * gtVar)
    })

    ## project if necessary
    if(doProjection){
        idsProj <- colnames(gt)[!(colnames(gt) %in% idsPca)]
        proj <- t(gt1[, idsProj]) %*% pca$u[, 1:pcs]
        proj <- t(t(proj) / pca$d[1:pcs])
        rownames(proj) <- idsProj
        res <- rbind(res, proj)
    }

    expl <- round(100 * pca$d^2 / sum(pca$d^2), 2)[1:pcs]
    return(list(summary = list(pca = res, explained = expl, loadings = l[, 1:pcs]), raw = pca ))
}


## faster version using rARPACK, but restricted to small number of PCAs
## % variance explained not calculated since only partial svd is calculated
getPcaGTFast <- function(gt, idsPca = colnames(gt), pcs = 10, var = "geno"){
  require(rARPACK)

  doProjection <- !all(colnames(gt) %in% idsPca)

  ## filter variant & missing snps in samples for building pca space
  nNonMissing <- rowSums(!is.na(gt[, idsPca]))
  idxFixed <- (rowSums(gt[, idsPca], na.rm = TRUE) / (2 * nNonMissing)) %in% c(0,1) ## fixed snps in nonmissing genotypes
  idxMissing <- rowSums(is.na(gt[, idsPca])) == length(idsPca) ## snps with all missing data
  gt <- gt[!idxMissing & !idxFixed,]

  ## center & normalize
  gtMean <- rowMeans(gt[, idsPca], na.rm = TRUE)
  if(var == "general"){
    gtVar <- var(gtMean)
  }else if(var == "geno"){
    gtVar <- (gtMean / 2)*(1 - gtMean / 2)
  }else{
    cat("ERROR: invalid data type!", fill=TRUE)
    quit("no")
  }
  gtStd <- sqrt(gtVar)
  gt1 <- (gt - gtMean) / gtStd
  gt1[is.na(gt1)] <- 0 ## set missing values to have mean gt (0 after normalizing)

  ## make pca
  pca <- svds(gt1[, idsPca], k = pcs)
  rownames(pca$u) <- rownames(gt)
  rownames(pca$v) <- idsPca

  res <- pca$v
  colnames(res) <- paste("PC", 1:pcs, sep = "")

  ## project if necessary
  if(doProjection){
    idsProj <- colnames(gt)[!(colnames(gt) %in% idsPca)]
    proj <- t(gt1[, idsProj]) %*% pca$u
    proj <- t(t(proj) / pca$d[1:pcs])
    rownames(proj) <- idsProj
    res <- rbind(res, proj)
  }
  
  return(list(summary = list(pca = res), raw = pca ))
}


#################################################################################
# Print parameters
cat('# Input file:', opt$in_file, fill=TRUE)
cat('# Has input file a header:', !opt$no_header, fill=TRUE)
cat('# Is input matrix symmetric?:', opt$data_symm, fill=TRUE)
cat('# Method:', opt$method, fill=TRUE)
opt$method <- toupper(opt$method)
cat('# Number components to calculate:', opt$n_comp, fill=TRUE)
cat('# Reference samples:', opt$s_proj, fill=TRUE)
cat('# Projection samples:', opt$s_proj, fill=TRUE)
cat('# Excluded samples:', opt$s_excl, fill=TRUE)
cat('# Type of data:', opt$data, fill=TRUE)
opt$data <- tolower(opt$data)
cat('# Missing data char:', opt$miss_char, fill=TRUE)
cat('# Out file:', opt$out_file, fill=TRUE)


### Read data
data <- read.table(opt$in_file, row.names=1, header=!opt$no_header, stringsAsFactors=FALSE, check.names=FALSE, na.strings=c(opt$miss_char))

if(opt$no_header && opt$data_symm)
  colnames(data) = rownames(data)
if(ncol(data) > nrow(data)){
  cat('# Found ',ncol(data),' columns and ',nrow(data),' rows. Assuming SNPs in columns...', fill=TRUE)
  data <- t(data)
}

# Exclude samples
if( !is.null(opt$s_excl) ){
  cat("# \tExcluding samples...", fill=TRUE)
  if(file.exists(opt$s_excl))
    opt$s_excl <- read.table(opt$s_excl, header=FALSE, stringsAsFactors=FALSE, check.names=FALSE)[,1]
  else
    opt$s_excl <- unlist(strsplit(opt$s_excl, ","))
  data <- data[!(rownames(data) %in% opt$s_excl),]
  if(opt$data_symm)
    data <- data[,!(colnames(data) %in% opt$s_excl)]
}

# Set default Reference and Projection samples
if(is.null(opt$s_ref)){
  opt$s_ref <- ncol(data)
}
if(is.null(opt$proj)){
  opt$s_proj <- ncol(data)
}

if(opt$method == "MDS"){
  if(ncol(data) != nrow(data))
    stop("\tERROR: Number of rows and columns do not match!", fill=TRUE)

  colnames(data) <- rownames(data);
  mds <- cmdscale(as.matrix(data), k=opt$n_comp, eig=TRUE, add=TRUE)

  #rank of the matrix...independent eigenvectors
  #r <- rank(mds.x);
  # get the % explained by each dimension: TO CHECK!!!!
  explained <- 100 * mds$eig / sum(mds$eig);

  out_matrix <- mds$points
  colnames(out_matrix) <- paste("D",1:opt$n_comp,sep="")
  colnames(out_matrix) <- paste(colnames(out_matrix), explained[1:opt$n_comp], sep="_")
}else{
  # Select samples
  s_ref <- colnames(data)[1:opt$s_ref]
  s_proj <- c()
  if(opt$s_ref < ncol(data)){
    s_proj <- colnames(data)[opt$s_ref+1:opt$s_proj]
  }
  
  # Calculate PCA
  if(opt$method == "PCA"){
    pca <- getPcaGT(data[, c(s_ref,s_proj)], idsPca=s_ref, pcs=opt$n_comp, var=opt$data)
  }else if(opt$method == "FASTPCA"){
    pca <- getPcaGTFast(data[, c(s_ref,s_proj)], idsPca=s_ref, pcs=opt$n_comp, var=opt$data)
  }
  out_matrix <- pca$summary$pca
  colnames(out_matrix) <- paste(colnames(out_matrix), pca$summary$explained, sep="_")
  rownames(out_matrix) <- rownames(data)
}

# Save results
write.table(out_matrix, opt$out_file, quote=FALSE, sep="\t")
