#!/usr/bin/env Rscript

## Run mlm or manova (with and w/o GAMMA transformation) in simulated datasets

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
  make_option(c("-k","--kinship"), type="character", 
              help="kinship matrix obtained by simulateGT.nf", metavar="character"),
  make_option(c("-v", "--vc"), type="character", 
              help="Variance components computed by GEMMA", metavar="character"),
  make_option(c("-c", "--covariates"), type="character", 
              help="Covariates (genotype PCs computed by plink2)", metavar="character"),
  make_option(c("-n", "--number"), type="numeric", default=20,
	      help="Number of PCs used to correct", metavar="numeric"),
  make_option(c("-t", "--transformation"), type="character", default="none",
              help="Transformation: none or GAMMA [default %default]", metavar="character"),
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
X <- as.matrix(BEDMatrix(opt$geno, simple_names = T))

id <- colnames(X)

# Apply same filter as in GEMMA regarding MAF (default)
maf <- apply(X, 2, function(x){
  tbl <- table(x)
  maf <- sum(c(tbl[3], tbl[2]/2), na.rm = T) / sum(tbl, na.rm = T)
})
X <- X[, maf >= 0.01, drop = F]
id <- id[maf >= 0.01]

transf <- opt$t
covariate_file <- opt$covariates

if(transf == "GAMMA"){
  Rg <- fread(opt$kinship, data.table = F)
  vc <- colMeans(read.table(opt$vc))

  # GAMMA functions
  eigen_solve <- function(K, tol = 1e-12) {
    a <- eigen(K)$vectors
    b <- eigen(K)$values
    b[b < tol] <- tol
    b <- 1/sqrt(b)
    return(a%*%diag(b)%*%t(a))
  }
  rotate <- function(A, Sigma) {
    U <- eigen_solve(Sigma)
    tU <- t(U)
    UA <- tU%*%A
    return(UA)
  }
  
  # Transform phenotypes and genotypes
  Vg <- vc[1]; Ve <- vc[2] 
  Sigma <- Vg*Rg + Ve*diag(nrow(Y))
  Y <- rotate(Y, Sigma)		        
  X <- rotate(X, Sigma)

} else if (transf == "PCA"){
  # Subset covariates
  k <- opt$number
  covariates <- read.table(covariate_file)[, -c(1:2)]
  covariates <- covariates[, 1:k]
}

## 1. Run mlm/manova


 if (transf != "PCA"){
   res <- apply(X, 2, function(x){mlm(Y ~ x)$aov.tab[1,6]})
 } else {
   res <- apply(X, 2, function(x){mlm(Y ~ ., data = data.frame(x, covariates))$aov.tab[1,6]})
 } 
 res <- cbind.data.frame(id, res)
 write.table(res, file = opt$mlm, col.names = F, row.names = F, quote = F, sep = "\t")

 if(!is.null(opt$manova)){
   if (transf != "PCA"){
     res_manova <- apply(X, 2, function(x){tryCatch( {summary(manova(Y ~ x))$stats[1,6]}, 
	     error = function(e){return(NA)} ) })
   } else {
     res_manova <- apply(X, 2, function(x){tryCatch( {summary(manova(Y ~ ., data = data.frame(x, covariates)))$stats[1,6]},
	     error = function(e){return(NA)} )})
   }
   res_manova <- cbind.data.frame(id, res_manova)
   write.table(res_manova, file = opt$manova, col.names = F, row.names = F, quote = F, sep = "\t")
 } 
