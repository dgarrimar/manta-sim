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
              help="Transformation: none, PCA or GAMMA [default %default]", metavar="character"),
  make_option(c("--mlm"), type="character", 
              help="Output (MLM p-values) file name)", metavar="character"),
  make_option(c("--manova"), type="character",
              help="Optional output (MANOVA p-values) file name", metavar="character"),
  make_option(c("--scale"), type="character", action = "store_true", default = FALSE,
              help="Scale response variables [default %default]", metavar="logical"),
  make_option(c("--runtime"), type="character", action = "store_true", default = FALSE,
              help="Report running time for each analysis step [default %default]", metavar="logical")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null (opt$geno) || is.null (opt$pheno) || is.null(opt$mlm) ){
    print_help(opt_parser)
    stop("Required I/O files must be supplied\n", call.=FALSE)
}

set.seed(123)

t0_readXY <- Sys.time()
Y <- as.matrix(fread(opt$pheno, data.table = F))
X <- as.matrix(BEDMatrix(opt$geno, simple_names = T))
t1_readXY <- Sys.time()

t0_maf <- Sys.time()
id <- colnames(X)

# Apply same filter as in GEMMA regarding MAF (default)
maf <- apply(X, 2, function(x){
    tbl <- table(x)
    maf <- sum(c(tbl[3], tbl[2]/2), na.rm = T) / sum(tbl, na.rm = T)
})
X <- X[, maf >= 0.01, drop = F]
id <- id[maf >= 0.01]
t1_maf = Sys.time()

transf <- opt$t
covariate_file <- opt$covariates

if(transf == "GAMMA"){
    t0_readRg <- Sys.time()
    Rg <- fread(opt$kinship, data.table = F)
    t1_readRg <- Sys.time()
  
    t0_GAMMA = Sys.time()
    vc <- colMeans(read.table(opt$vc))

    # GAMMA functions
    eigen_solve <- function(K, tol = 1e-12) {
        a <- eigen(K)$vectors
        b <- eigen(K)$values
        b[b < tol] <- tol
        b <- 1/sqrt(b)
        return(tcrossprod(a%*%diag(b), a))
    }
  
    # Transform phenotypes and genotypes
    Vg <- vc[1]; Ve <- vc[2] 
    Sigma <- Vg*Rg + Ve*diag(nrow(Y))
    U <- eigen_solve(Sigma)

    Y <- crossprod(U, Y)        
    X <- crossprod(U, X)
    t1_GAMMA = Sys.time()

} else if (transf == "PCA"){
  t0_cov <- Sys.time()    
  # Subset covariates
  k <- opt$number
  covariates <- fread(covariate_file, data.table = F)[, -c(1:2)]
  covariates <- covariates[, 1:k]
  t1_cov <- Sys.time()
}

## 1. Run mlm/manova

t0_mlm <- Sys.time()
if (transf != "PCA"){
   res <- apply(X, 2, function(x){tryCatch( {mlm(Y ~ x)$aov.tab[1,6]}, 
             error = function(e){return(NA)} ) })
} else {
   res <- apply(X, 2, function(x){tryCatch( {mlm(Y ~ ., data = data.frame(x, covariates))$aov.tab[1,6]},
             error = function(e){return(NA)} ) })
} 
res <- cbind.data.frame(id, res)
write.table(res, file = opt$mlm, col.names = F, row.names = F, quote = F, sep = "\t")
t1_mlm <- Sys.time()

if(!is.null(opt$manova)){
    t0_manova <- Sys.time()
    if (transf != "PCA"){
        res_manova <- apply(X, 2, function(x){tryCatch( {summary(manova(Y ~ x))$stats[1,6]}, 
                error = function(e){return(NA)} ) })
    } else {
        res_manova <- apply(X, 2, function(x){tryCatch( {summary(manova(Y ~ ., data = data.frame(x, covariates)))$stats[1,6]},
                error = function(e){return(NA)} )})
    }
    res_manova <- cbind.data.frame(id, res_manova)
    write.table(res_manova, file = opt$manova, col.names = F, row.names = F, quote = F, sep = "\t")
    t1_manova <- Sys.time()
}

# Runtime
if (opt$runtime){
    N <- nrow(Y)
    t_readXY <- as.numeric(difftime(t1_readXY, t0_readXY, units = "secs"))
    cat(sprintf("%s\t%s\tread_XY\t%s\n", N, transf, round(t_readXY)))
    t_maf <- as.numeric(difftime(t1_maf, t0_maf, units = "secs"))
    cat(sprintf("%s\t%s\tmaf\t%s\n", N, transf, round(t_maf)))
    if(transf == "GAMMA"){
        t_readRg <- as.numeric(difftime(t1_readRg, t0_readRg, units = "secs"))
        cat(sprintf("%s\t%s\tread_kinship\t%s\n", N, transf, round(t_readRg)))
        t_GAMMA <- as.numeric(difftime(t1_GAMMA, t0_GAMMA, units = "secs"))
        cat(sprintf("%s\t%s\tGAMMA\t%s\n", N, transf, round(t_GAMMA)))
    } else if (transf == "PCA"){
        t_cov <- as.numeric(difftime(t1_cov, t0_cov, units = "secs"))
        cat(sprintf("%s\t%s\tcov\t%s\n", N, transf, round(t_cov)))
    }
    t_mlm <- as.numeric(difftime(t1_mlm, t0_mlm, units = "secs"))
    cat(sprintf("%s\t%s\tmlm\t%s\n", N, transf, round(t_mlm)))
    if(!is.null(opt$manova)){
        t_manova <- as.numeric(difftime(t1_manova, t0_manova, units = "secs"))
        cat(sprintf("%s\t%s\tmanova\t%s\n", N, transf, round(t_manova)))
    }
}
