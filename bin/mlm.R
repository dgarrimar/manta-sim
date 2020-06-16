#!/usr/bin/env Rscript

## Run mlm (with and w/o GAMMA transformation) in simulated datasets

##  0. Load libraries, parse arguments

library(optparse)
library(mlm)
library(data.table)

option_list = list(
  make_option(c("-p","--pheno"), type="character", default=NULL,
              help="phenotype data obtained by simulatePT.R", metavar="character"),
  make_option(c("-g","--geno"), type="character", default=NULL,
              help="genotype data obtained by simulateGT.nf", metavar="character"),
  make_option(c("-k","--kinship"), type="character", default=NULL,
              help="kinship matrix obtained by simulateGT.nf", metavar="character"),
  make_option(c("-v", "--vc"), type="character", default=NULL,
              help="Variance components computed by GEMMA", metavar="character"),
  make_option(c("-c", "--covariates"), type="character", default=NULL,
              help="Covariates (genotype PCs computed by plink2)", metavar="character"),
  make_option(c("-t", "--transformation"), type="character", default="none",
              help="Transformation: none or GAMMA", metavar="character"),
  make_option(c("-l", "--level"), type="numeric", default=0.05,
              help="Significance level", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output (MLM p-values) file name", metavar="character")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null (opt$geno) || is.null (opt$pheno) || is.null(opt$output) ){
  print_help(opt_parser)
  stop("Required I/O files must be supplied\n", call.=FALSE)
}

set.seed(123)

Y <- as.matrix(fread(opt$pheno, data.table = F))
X <- fread(opt$geno, data.table = F, sep = ",")

rs <- X[, 1]
X[, c(1:3)] <- NULL
X <- t(X)

# Apply same filter as in GEMMA regarding MAF (default)
maf <- apply(X, 2, function(x){
  tbl <- table(x)
  maf <- sum(c(tbl[3], tbl[2]/2), na.rm = T) / sum(tbl, na.rm = T)
})
X <- X[, maf >= 0.01]
rs <- rs[maf >= 0.01]

########################
########################
X <- X[, 1:100]
rs <- rs[1:100]
########################
########################

transf <- opt$t
covariate_file <- opt$covariates

if(transf == "GAMMA"){
  Rg <- opt$kinship
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
  k <- 3
  covariates <- read.table(covariate_file)[,-1]
  covariates <- covariates[, c(1:k)]
}

## 1. Run mlm

 # Run
 res <- c()
 if (is.null(covariate_file)){
   res <- apply(X, 2, function(x){mlm(Y ~ x)$aov.tab[1,6]})
 } else {
   res <- apply(X, 2, function(x){mlm(Y ~ ., data = data.frame(x, covariates))$aov.tab[1,6]})
 }
 
 res <- cbind.data.frame(rs, res)
 
 # Print TIE/POWER
 cat(mean(res[,2] < opt$level )) 
 
 # Save
 write.table(res, file = opt$output, col.names = F, row.names = F, sep = "\t")

