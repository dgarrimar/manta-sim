#!/usr/bin/env Rscript

## Simulate phenotypes based on FPC's thesis

##  0. Load libraries, parse arguments

library(optparse)
library(MASS)
library(data.table)

option_list = list(
  make_option(c("-r","--rpt"), type="numeric", default=0,
	      help="repeat id (seed) [default %default]", metavar="numeric"),
  make_option(c("-n","--nb_samples"), type="numeric", default=1000,
              help="Total number of samples [default %default]", metavar="numeric"),
  make_option(c("-q", "--nb_responses"), type="numeric", default=3,
              help="Number of dependent variables [default %default]", metavar="numeric"),
  make_option(c("-s", "--nb_causal"), type="numeric", default=0,
              help="Number of causal variants [default %default]", metavar="numeric"),
  make_option(c("--PTgen"), type="character", default="matrixNorm",
              help="Phenotype data generation: matrixNormal or copula [default %default]", metavar="character"),
  make_option(c("--geno"), type="character", default=NULL,
              help="genotype data obtained by simulateGT.nf", metavar="character"),
  make_option(c("--kinship"), type="character", default=NULL,
              help="kinship matrix obtained by simulateGT.nf", metavar="character"),
  make_option(c("--hs2"), type="numeric", default=0.01,
              help="Average fraction of variance explained by causal variants across traits [default %default]", metavar="numeric"),
  make_option(c("--hg2"), type="numeric", default=0.4,
              help="Average fraction of variance explained by causal variants across traits [default %default]", metavar="numeric"),
  make_option(c("--alphaG"), type="numeric", default=0.5,
              help="Fraction of signal from the relatedness contribution shared across traits [default %default]", metavar="numeric"),
  make_option(c("--lambda"), type="numeric", default=0.6,
              help="Fraction of structured noise [default %default]", metavar="numeric"),
  make_option(c("--alphaH"), type="numeric", default=0.3,
              help="Fraction of structured noise that is shared across traits [default %default]", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output (simulated phenotype) file name", metavar="character"),
  make_option(c("-i", "--id_file"), type="character", default=NULL,
              help="IDs of variants generated under H1, if any", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$output) || is.null (opt$geno) || is.null (opt$kinship) || is.null(opt$id_file)){
  print_help(opt_parser)
  stop("Required I/O files must be supplied\n", call.=FALSE)
}

set.seed(opt$rpt)

n <- opt$nb_samples
q <- opt$nb_responses                
s <- opt$nb_causal
PTgen <- opt$PTgen
geno <- opt$geno
kinship <- opt$kinship
hs2 <- opt$hs2
id_file <- opt$id_file
if(s == 0) {
  hs2 <- 0
} 
hg2 <- opt$hg2
alphaG <- opt$alphaG
lambda <- opt$lambda
alphaH <- opt$alphaH
outfile <- opt$output

if (PTgen == "matrixNorm"){
  
  rmatnorm_K <- function(M, U, V){
    # Sample from Matrix Normal distribution via Kronecker
    # Slower than via Cholesky
    n <- nrow(M)
    q <- ncol(M)
    X <- matrix(mvrnorm(mu = c(M), Sigma = kronecker(V, U)), n, q)
    return(X)
  }
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
  
  M <- matrix(0, n, q)   # Means

  ## 1. Effect of the SNP
  
   if (s != 0){
     # Load genotypes
     S <- fread(geno, data.table = FALSE, sep = ",")
     p <- nrow(S) # Total number of variants
          
     # Select causal SNP(s)
     sel <- sample(1:p, size = s)
     ids <- S[sel, 1]
     write.table(ids, file = id_file, col.names = F, row.names = F, quote = F)

     X <- scale(t(S[sel, -c(1:3), drop = F])) # Standarized genotypes
     
     # Generate effects
     B <- matrix(sample(c(-1,1), replace = T, size = s*q), s, q)
     
     # Multiply
     XB <- X %*% B 
     XB <- XB / sqrt(mean(diag(cov(XB)))) * sqrt(hs2) 
   } else {
     XB <- matrix(0, n, q)
     ignore <- file.create(id_file)
   }
  
  ## 2. Relatedness signal 
  
   ## Kinship (from standarized genotypes)
   Rg <- as.matrix(fread(kinship, data.table = FALSE, sep = "\t"))
     
   ## 2.1. Shared
    aaT <- tcrossprod(rnorm(q))
    Gs <- rmatnorm_K(M = M, U = Rg, V = aaT) 
    Gs <- Gs / sqrt(mean(diag(cov(Gs)))) * sqrt( hg2*alphaG ) # Rescale
   
   ## 2.2. Independent
    C <- diag(rnorm(q)^2)
    Gi <- rmatnorm_K(M = M, U = Rg, V = C)
    Gi <- Gi / sqrt(mean(diag(cov(Gi)))) * sqrt(hg2*(1-alphaG)) # Rescale
    
  ## 3. Residuals
   
   k <- 10
   LLT <- tcrossprod(matrix(rnorm(n*k), n, k))
    
   ## 3.1. Structured noise (shared)
  
    aHaHT <- tcrossprod(rnorm(q))
    Hs <- rmatnorm_K(M = M, U = LLT, V = aHaHT)
    Hs <- Hs / sqrt(mean(diag(cov(Hs)))) * sqrt( lambda*alphaH*(1-hs2-hg2) ) # Rescale 
    
   ## 3.2. Structured noise (independent)
  
    CH <- diag(rnorm(q)^2) 
    Hi <- rmatnorm_K(M = M, U = LLT, V = CH)
    Hi <- Hi / sqrt(mean(diag(cov(Hi)))) * sqrt( lambda*(1-alphaH)*(1-hs2-hg2) ) # Rescale
  
   ## 3.3. Unstructured noise
    D <- 10 # Number of products
    E <- 1
    for (d in 1:D){
      E <- E * matrix(rnorm(n*q), n, q)
    }
    E <- E / sqrt(mean(diag(cov(E)))) * sqrt( (1-lambda)*(1-hs2-hg2) )
   
  ## 4. Build Y
    
    Y <- XB + (Gs + Gi) + (Hs + Hi) + E
    
  ## 5. Save
    
    write.table(Y, file = outfile, col.names = F, row.names = F, sep = "\t")
    
} else {
  stop("Not implemented yet!")
}
