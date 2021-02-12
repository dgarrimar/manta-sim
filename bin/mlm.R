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
  make_option(c("-c", "--covariates"), type="character", default=NULL,
              help="Covariates (genotype PCs computed by plink2) [default %default]", metavar="character"),
  make_option(c("-k", "--number"), type="numeric", default=20,
              help="Number of PCs used to correct", metavar="numeric"),
  make_option(c("--mlm"), type="character", 
              help="Output (MLM p-values) file name)", metavar="character"),
  make_option(c("--manova"), type="character",
              help="Optional output (MANOVA p-values) file name", metavar="character"),
  make_option(c("--scale"), type="character", action = "store_true", default = FALSE,
              help="Scale response variables [default %default]", metavar="logical"),
  make_option(c("--runtime"), type="character", action = "store_true", default = FALSE,
              help="Report running time for each analysis step [default %default]", metavar="logical"),
  make_option(c("-i", "--id_replicate"), type="numeric", default=1, 
              help="Replicate id (used if --runtime)", metavar="numeric")
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

covariate_file <- opt$covariates

if (!is.null(covariate_file)){
  t0_cov <- Sys.time()    
  # Subset covariates
  k <- opt$number
  covariates <- fread(covariate_file, data.table = F)[, -c(1:2)]
  covariates <- covariates[, 1:k]
  t1_cov <- Sys.time()
}

## 1. Run mlm/manova

if(opt$scale){ # WARNING: Asymptotic null may not hold
    Y <- scale(Y)
}

t0_mlm <- Sys.time()
p <- ncol(X)
res <- rep(NA, p)
if (is.null(covariate_file)){
    for (snp in 1:p) {
        res[snp] <- tryCatch( {mlm(Y ~ ., data = data.frame(snp = X[, snp]), type = "I", subset = "snp")$aov.tab[1,6]}, error = function(e){return(NA)} )
    }
} else {
    for (snp in 1:p) {
        res[snp] <- tryCatch( {mlm(Y ~ ., data = data.frame(covariates, snp = X[, snp]), type = "I", subset = "snp")$aov.tab[1,6]}, 
                              error = function(e){return(NA)} )
    }
} 
res <- cbind.data.frame(id, res)
write.table(res, file = opt$mlm, col.names = F, row.names = F, quote = F, sep = "\t")
t1_mlm <- Sys.time()

if(!is.null(opt$manova)){
    t0_manova <- Sys.time()
    res_manova <- rep(NA, p)
    if (is.null(covariate_file)){
        for (snp in 1:p) {
            res_manova[snp] <- tryCatch( {summary(manova(Y ~ ., data.frame(snp = X[, snp])))$stats[1,6]}, error = function(e){return(NA)} )
        }
    } else {
        for (snp in 1:p) {
            res_manova[snp] <- tryCatch( {summary(manova(Y ~ ., data = data.frame(covariates, snp = X[, snp])))$stats["snp", 6]}, 
                                         error = function(e){return(NA)} )
        }
    }
    res_manova <- cbind.data.frame(id, res_manova)
    write.table(res_manova, file = opt$manova, col.names = F, row.names = F, quote = F, sep = "\t")
    t1_manova <- Sys.time()
}

# Runtime
if (opt$runtime){
    N <- nrow(Y)
    Q <- ncol(Y)
    R <- opt$id_replicate
    t_readXY <- as.numeric(difftime(t1_readXY, t0_readXY, units = "secs"))
    cat(sprintf("%s\t%s\t%s\t%s\tread_XY\t%s\n", N, Q, R, 'MLM', round(t_readXY)))
    t_maf <- as.numeric(difftime(t1_maf, t0_maf, units = "secs"))
    cat(sprintf("%s\t%s\t%s\t%s\tmaf\t%s\n", N, Q, R, 'MLM', round(t_maf)))
    if (!is.null(covariate_file)){
        t_cov <- as.numeric(difftime(t1_cov, t0_cov, units = "secs"))
        cat(sprintf("%s\t%s\t%s\t%s\tcov\t%s\n", N, Q, R, 'MLM', round(t_cov)))
    }
    t_mlm <- as.numeric(difftime(t1_mlm, t0_mlm, units = "secs"))
    cat(sprintf("%s\t%s\t%s\t%s\tmlm\t%s\n", N, Q, R, 'MLM', round(t_mlm)))
    if(!is.null(opt$manova)){
        t_manova <- as.numeric(difftime(t1_manova, t0_manova, units = "secs"))
        cat(sprintf("%s\t%s\t%s\t%s\tmanova\t%s\n", N, Q, R, 'MLM', round(t_manova)))
    }
}
