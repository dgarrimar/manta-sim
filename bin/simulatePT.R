#!/usr/bin/env Rscript

## Simulate phenotypes based on FPC's thesis

##  0. Load libraries, parse arguments

library(optparse)
library(MASS)
library(data.table)
library(BEDMatrix)
library(copula)

option_list = list(
  make_option(c("-s","--seed"), type="numeric", default=0,
              help="Set seed for random processes [default %default]", metavar="numeric"),
  make_option(c("-n","--nb_samples"), type="numeric", default=1000,
              help="Total number of samples [default %default]", metavar="numeric"),
  make_option(c("-q", "--nb_responses"), type="numeric", default=3,
              help="Number of dependent variables [default %default]", metavar="numeric"),
  make_option(c("--PTgen"), type="character", default="norm-0-1",
              help="Phenotype data generation: matrixNormal or copula [default %default]", metavar="character"),
  make_option(c("--geno"), type="character", default=NULL,
              help="genotype data (single variant) obtained by simulateGT.nf", metavar="character"),
  make_option(c("--kinship"), type="character", default=NULL,
              help="kinship matrix obtained by simulateGT.nf", metavar="character"),
  make_option(c("--hs2"), type="numeric", default=0,
              help="Average fraction of variance explained by causal variants across traits [default %default]", metavar="numeric"),
  make_option(c("--hg2"), type="numeric", default=0,
              help="Average fraction of variance explained by causal variants across traits [default %default]", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output (simulated phenotype) file name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$output) || is.null (opt$geno) || is.null (opt$kinship)){
  print_help(opt_parser)
  stop("Required I/O files must be supplied\n", call.=FALSE)
}

set.seed(opt$seed)

n <- opt$nb_samples
q <- opt$nb_responses                
PTgen <- opt$PTgen
geno <- opt$geno
kinship <- opt$kinship
hs2 <- opt$hs2
hg2 <- opt$hg2
outfile <- opt$output

rmatnorm_C <- function(M, U, V, tol = 1e-12){
  # Sample from Matrix Normal distribution via Cholesky
  # https://en.wikipedia.org/wiki/Matrix_normal_distribution#Drawing_values_from_the_distribution
  # When many eigenvalues are ~ 0, better use rmatnorm_K
  a <- nrow(M)
  b <- ncol(M)
  Z <- matrix(rnorm(a*b, 0, 1), a, b)
  L1 <- suppressWarnings(chol(U, pivot = T))
  L1 <- L1[, order(attr(L1, "pivot"))]
  L2 <- suppressWarnings(chol(V, pivot = T))
  L2 <- L2[, order(attr(L2, "pivot"))]
  return(M + crossprod(L1, Z) %*% L2)
}

cop <- function(n, sigma, distdef){
  
  # Get correlation matrix
  q <- nrow(sigma)
  D <- diag(q)*sqrt(diag(sigma))
  R <- solve(D)%*%sigma%*%solve(D)
  # https://math.stackexchange.com/questions/186959/correlation-matrix-from-covariance-matrix
  
  # Build and sample from copula
  distrib <- unlist(strsplit(distdef, "-"))[1]
  params <- as.numeric(unlist(strsplit(distdef, "-"))[-1])
  mar <- switch(distrib[1],
                "unif" = list(min = params[1], max = params[2]),
                "gamma" = list(shape = params[1], scale = params[2]),
                "beta" = list(shape1 = params[1], shape2 = params[2]),
                "t" = list(df = params[1]),
                "exp" = list(rate = params[1]),
                "lnorm" = list(meanlog = params[1], sdlog = params[2]),
                "binom" = list(size = params[1], p = params[2]),
                "norm" = list(mean = params[1], sd = params[2]))
  myCop <- normalCopula(param = P2p(R), dim = q, dispstr = "un")
  myMvd <- mvdc(copula = myCop, margins = rep(distrib[1], q),
                paramMargins = rep(list(mar), q))
  E <- rMvdc(n, myMvd)
  E <- t(apply(scale(E),1,function(e){e*sqrt(diag(sigma))}))
  return(E)
}

## 1. Effect of the SNP
  
if (hs2 != 0){
  
  # Load genotypes
  X <- as.matrix(BEDMatrix(geno, simple_names = T))
  X <- scale(X) # Standarized genotypes

  # Generate effects
  B <- matrix(sample(c(-1,1), replace = T, size = q), 1, q)
  XB <- X %*% B 
  XB <- XB / sqrt(mean(diag(cov(XB)))) * sqrt(hs2) # Rescale

} else {
  
  XB <- matrix(0, n, q)

}
  
## 2. Relatedness signal 

# Kinship 
Rg <- as.matrix(fread(kinship, data.table = FALSE, sep = "\t"))
     
# Shared genotype effect
A <- matrix(rnorm(q^2), q, q)
AAT <- tcrossprod(A)
G <- rmatnorm_C(M = matrix(0, n, q), U = Rg, V = AAT) 
G <- G / sqrt(mean(diag(cov(G)))) * sqrt(hg2) # Rescale

## 3. Residuals

B <- matrix(rnorm(q^2), q, q)
BBT <- tcrossprod(B)
E <- cop(n, BBT, PTgen)
   
E <- E / sqrt(mean(diag(cov(E)))) * sqrt( (1-hs2-hg2) ) # Rescale
   
## 4. Build Y
    
Y <- XB + G + E
    
## 5. Save
    
write.table(Y, file = outfile, col.names = F, row.names = F, sep = "\t")
