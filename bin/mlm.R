#!/usr/bin/env Rscript

## Run mlm or manova in simulated datasets

##  0. Load libraries, parse arguments

library(optparse)
library(mlm)
library(data.table)
library(BEDMatrix)

option_list = list(
  make_option(c("-p","--pheno"), type="character", 
              help="phenotype data obtained by simulatePT.R", metavar="character"),
  make_option(c("-g","--geno"), type="character", 
              help="genotype data obtained by simulateGT.nf", metavar="character"),
  make_option(c("-c", "--covariates"), type="character", default = NULL,
              help="Covariates (genotype PCs computed by plink2)", metavar="character"),
  make_option(c("-k", "--number"), type="numeric", default=20,
              help="Number of PCs used to correct", metavar="numeric"),
  make_option(c("--mlm"), type="character", 
              help="Output (MLM p-values) file name)", metavar="character"),
  make_option(c("--manova"), type="character",
              help="Optional output (MANOVA p-values) file name", metavar="character"),
  make_option(c("--scale"), type="character", action = "store_true", default = FALSE,
              help="Scale response variables [default %default]", metavar="logical")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null (opt$geno) || is.null (opt$pheno) || is.null(opt$mlm) ){
    print_help(opt_parser)
    stop("Required I/O files must be supplied\n", call.=FALSE)
}

set.seed(123)

Y <- as.matrix(fread(opt$pheno, data.table = F))

if(sum(Y) == 0){ # Variant has been filtered due to MAF < 0.01 as in GEMMA (default)
  write.table(NULL, file = opt$mlm, col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(NULL, file = opt$manova, col.names = F, row.names = F, quote = F, sep = "\t")
  quit('no', status = 0)
}

X <- as.matrix(BEDMatrix(opt$geno, simple_names = T))
id <- colnames(X)

covariate_file <- opt$covariates

if (!is.null(covariate_file)){
  # Subset covariates
  k <- opt$number
  covariates <- fread(covariate_file, data.table = F)[, -c(1:2)]
  covariates <- covariates[, 1:k]
}

## 1. Run mlm/manova

if(opt$scale){ # WARNING: Asymptotic null may not hold
    Y <- scale(Y)
}

if (is.null(covariate_file)){
    res <- tryCatch( {mlm(Y ~ ., data = data.frame(snp = X[, 1]), type = "I", subset = "snp")$aov.tab[1, 6]}, 
                     error = function(e){return(NA)} )
} else {
    res <- tryCatch( {mlm(Y ~ ., data = data.frame(covariates, snp = X[, 1]), type = "I", subset = "snp")$aov.tab[1, 6]}, 
                     error = function(e){return(NA)} )
}
 
res <- cbind(id, res)
write.table(res, file = opt$mlm, col.names = F, row.names = F, quote = F, sep = "\t")

if(!is.null(opt$manova)){
    if (is.null(covariate_file)){
        res_manova <- tryCatch( {summary(manova(Y ~ ., data.frame(snp = X[, 1])))$stats[1, 6]}, 
                                error = function(e){return(NA)} )
    } else {
        res_manova <- tryCatch( {summary(manova(Y ~ ., data = data.frame(covariates, snp = X[, 1])))$stats["snp", 6]}, 
                                error = function(e){return(NA)} )
    }
    res_manova <- cbind(id, res_manova)
    write.table(res_manova, file = opt$manova, col.names = F, row.names = F, quote = F, sep = "\t")
}

