#!/usr/bin/env Rscript

## Simulate phenotypes 

##  0. Load libraries, parse arguments, define functions

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
  make_option(c("--PTgen"), type="character", default="mvnorm",
              help="Error in phenotype data generation: 'mvnorm', 'simplex', 'multinom' or copula definition (e.g. 'unif-0-1') [default %default]", metavar="character"),
  make_option(c("--geno"), type="character", default=NULL,
              help="genotype data (single variant) obtained by simulateGT.nf", metavar="character"),
  make_option(c("--kinship"), type="character", default=NULL,
              help="kinship matrix obtained by simulateGT.nf", metavar="character"),
  make_option(c("--hs2"), type="numeric", default=0,
              help="Average fraction of variance explained by causal variants across traits [default %default]", metavar="numeric"),
  make_option(c("--hg2"), type="numeric", default=0,
              help="Average fraction of variance explained by causal variants across traits [default %default]", metavar="numeric"),
  make_option(c("--varG"), type="character", default="random",
              help="[PTgen: mvnorm, copula] Genetic variances: 'equal', 'unequal', 'random' [default %default]", metavar="character"),
  make_option(c("--varE"), type="character", default="random",
              help="[PTgen: mvnorm, copula] Error variances: 'equal', 'unequal', 'random' [default %default]", metavar="character"),
  make_option(c("--vGr"), type="character", default=2,
              help="[PTgen: mvnorm, copula; varG: unequal] max/min ratio of genetic variances [default %default]", metavar="numeric"),
  make_option(c("--vEr"), type="character", default=2,
              help="[PTgen: mvnorm, copula; varE: unequal] max/min ratio of error variances [default %default]", metavar="numeric"),
  make_option(c("--corG"), type="numeric", default=0,
              help="[PTgen: 'mvnorm', 'copula' & varG not 'random'] Genetic correlations [default %default]", metavar="numeric"),
  make_option(c("--corE"), type="numeric", default=0,
              help="[PTgen: 'mvnorm', 'copula' & varG not 'random'] Error correlations [default %default]", metavar="numeric"),
  make_option(c("-p","--p_loc"), type="numeric", default=1,
              help="[PTgen: 'dirichlet', 'multinom'] parameter location [default %default]", 
              metavar="numeric"),
  make_option(c("-t", "--transf"), type="character", default="none",
              help="Transformation of response variables [default %default]", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output (simulated phenotype) file name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$output) || is.null (opt$geno) || is.null (opt$kinship)){
  print_help(opt_parser)
  stop("Required I/O files must be supplied\n", call.=FALSE)
}

n <- opt$nb_samples
q <- opt$nb_responses                
PTgen <- opt$PTgen
geno <- opt$geno
kinship <- opt$kinship
hs2 <- opt$hs2
hg2 <- opt$hg2
varG <- opt$varG
vGr <- opt$vGr
corG <- opt$corG
varE <- opt$varE
vEr <- opt$vEr
corE <- opt$corE
ploc <- opt$p_loc
transf <- opt$transf
outfile <- opt$output

rmatnorm_C <- function(M, U, V, tol = 1e-12){
  # Sample from Matrix Normal distribution via Cholesky
  # https://en.wikipedia.org/wiki/Matrix_normal_distribution#Drawing_values_from_the_distribution
  # Fast but requires U and V positive definite. 
  # Works well for positive semidefinite matrices via pivoting, when not many eigenvalues are ~ 0.
  # Otherwise, better use Kronecker product

  a <- nrow(M)
  b <- ncol(M)
  Z <- matrix(rnorm(a*b, 0, 1), a, b)
  L1 <- suppressWarnings(chol(U, pivot = T))
  L1 <- L1[, order(attr(L1, "pivot"))]
  L2 <- suppressWarnings(chol(V, pivot = T))
  L2 <- L2[, order(attr(L2, "pivot"))]

  return(M + crossprod(L1, Z) %*% L2)
}

getCov <- function(q, v, c, u, tol = 1e-10){

  if (v == 'random') {
    return(tcrossprod(matrix(rnorm(q^2), q, q)))
  } else if(v == 'equal'){
    vars <- rep(1, q)
  } else if (v == 'unequal'){
    # vars <- (q:1)/sum(q:1)
    vars <- seq(from = 1, to = u, length.out = q)
  } else {
    stop(sprintf("Unknown option: Var = '%s'.", v))
  }
  R <- matrix(c, nrow = q, ncol = q)
  diag(R) <- rep(1, q)
  sigma <- R * tcrossprod(sqrt(vars))

  if(any(eigen(sigma, only.values = T)$values < tol)){
    stop("Covariance matrix should be positive definite.")
  }

  return(sigma)
}

sim.copula <- function(n, sigma, distdef){

   # Obtain correlation matrix
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
   E <- rMvdc(n, myMvd) # This has sigma = R
   E <- t(apply(scale(E),1,function(e){e*sqrt(diag(sigma))})) 

   return(E)
}

rdirichlet <- function(n, alpha) {
   l <- length(alpha)
   x <- matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE)
   sm <- x%*%rep(1,l)
   return(x/as.vector(sm))
}

set.seed(opt$seed)

## 1. Effect of the SNP
  
if (hs2 != 0){
  
   # Load genotypes
   X <- as.matrix(BEDMatrix(geno, simple_names = T))

   # Generate effects
   B <- matrix(sample(c(-1,1), replace = T, size = q), 1, q)
   XB <- X %*% B 
   XB <- XB / sqrt(mean(diag(cov(XB)))) * sqrt(hs2) # Rescale

} 
  
## 2. Relatedness signal 

if (hg2 != 0){
   
   # Load kinship 
   Rg <- as.matrix(fread(kinship, data.table = FALSE, sep = "\t")) 
  
   # Genotype effect
   G <- rmatnorm_C(M = matrix(0, n, q), U = Rg, V = getCov(q, varG, corG, vGr)) 
   G <- G / sqrt(mean(diag(cov(G)))) * sqrt(hg2) # Rescale
} 

## 3. Residuals

sigma <- getCov(q, varE, corE, vEr)

if(PTgen == 'mvnorm'){
    E <- mvrnorm(n = n, mu = rep(0, q), Sigma = sigma) 
} else if (PTgen == 'dirichlet'){
    x <- c(ploc, rep(1, q-1))
    p <- x/sum(x)
    E <- rdirichlet(n, p)
} else if (PTgen == 'multinom'){
    x <- c(ploc, rep(1, q-1))
    p <- x/sum(x)
    N <- 1000
    E <- t(rmultinom(n, N, p))
}else {
    E <- sim.copula(n, sigma, PTgen)
}

E <- E / sqrt(mean(diag(cov(E)))) * sqrt( (1-hs2-hg2) ) # Rescale
   
## 4. Build Y (transform)

if(hs2 != 0 && hg2 != 0){
    Y <- XB + G + E
} else if (hs2 == 0 && hg2 != 0){
    Y <- G + E
} else if (hs2 != 0 && hg2 == 0){
    Y <- XB + E
} else if (hs2 == 0 && hg2 == 0){
    Y <- E
}

if (transf == "sqrt"){
    Y <- sqrt(Y)
} else if (transf == "log"){
    Y <- log(Y + 1)
} else if (transf == "none"){
    
} else {
    stop(sprintf("Unknown option: transf = '%s'.", transf))
}  

    
## 5. Save    

write.table(Y, file = outfile, col.names = F, row.names = F, sep = "\t")
