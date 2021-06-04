#!/usr/bin/env Rscript

## Simulate phenotypes 

##  0. Load libraries, parse arguments, define functions

library(optparse)
library(MASS)
library(data.table)
library(BEDMatrix)
library(copula)
library(mvnfast)

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
  make_option(c("--b"), type="character", default="equal",
              help="[PTgen: mvnorm, copula] Effect type: 'equal', 'unequal', 'block' [default %default]", metavar="character"),
  make_option(c("--ub"), type="character", default=2,
              help="[PTgen: mvnorm, copula; b: unequal] max/min effect ratio [default %default]", metavar="numeric"),
  make_option(c("-p","--p_loc"), type="numeric", default=1,
              help="[PTgen: 'simplex', 'multinom'] parameter location [default %default]", 
              metavar="numeric"),
  make_option(c("--hk"), type="numeric", default=1, help="[PTgen: 'mvnorm', varE: 'random'] heteroscedasticity (variances) [default %default]", 
              metavar="numeric"),
  make_option(c("--chk"), type="numeric", default=1, help="[PTgen: 'mvnorm', varE 'random'] heteroscedasticity (covariances) [default %default]",
              metavar="numeric"),
  make_option(c("-t", "--transf"), type="character", default="none",
              help="Transformation of response variables [default %default]", metavar="numeric"),
  make_option(c("-f", "--fx"), type="character", default=NULL,
              help="Path to helper functions", metavar="character"),
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
b <- opt$b
ub <- opt$ub
ploc <- opt$p_loc
hk <- opt$hk
chk <- opt$chk
transf <- opt$transf
outfile <- opt$output

source(sprintf("%s/fx.R", opt$fx))

set.seed(opt$seed)

## 1. Effect of the SNP
  
if (hs2 != 0 || hk != 1 || chk != 1){
  
   # Load genotypes
   X <- as.matrix(BEDMatrix(geno, simple_names = T))

   # Generate effects
   B <- getBeta(q, if(PTgen %in% c("simplex", "multinom")){PTgen}else{b})
   
   XB <- X %*% B 
   XB <- XB / sqrt(mean(diag(cov(XB)))) * sqrt(hs2) # Rescale

} 
  
## 2. Relatedness signal 

if (hg2 != 0){
   
   # Load kinship 
   Rg <- as.matrix(fread(kinship, data.table = FALSE, sep = "\t")) 
  
   # Genotype effect
   G <- rmatnorm_C(M = matrix(0, n, q), U = Rg, V = getCov(q, varG, corG, vGr, rep(1, q))) 
   G <- G / sqrt(mean(diag(cov(G)))) * sqrt(hg2) # Rescale
} 

## 3. Residuals

if(PTgen == 'mvnorm'){
  
    if ( hk != 1 || chk != 1 ) {
        E <- matrix(NA, n, q)
        tbl <- table(X)
        if (length(tbl) > 1){
            hkv <- c((1+hk)/2, hk)
            chkv <- c((1+chk)/2, chk)
            sigma0 <- tcrossprod(matrix(rnorm(q^2), q, q))   
            sigma1 <- sigma0*hkv[1]/chkv[1]
            sigma2 <- sigma0*hkv[2]/chkv[2]
            if (chk != 1) {
       	       diag(sigma2) <- diag(sigma1) <- diag(sigma0)	
            } 
            E[X == 0, ] <- mvrnorm(n = tbl['0'], mu = rep(0, q), Sigma = sigma0)
            E[X == 1, ] <- mvrnorm(n = tbl['1'], mu = rep(0, q), Sigma = sigma1)
            if (length(tbl) == 3){
                E[X == 2, ] <- mvrnorm(n = tbl['2'], mu = rep(0, q), Sigma = sigma2)
            }
        } else {
            E <- matrix(0, n, q) 
        } 
    } else {
        sigma <- getCov(q, varE, corE, vEr, if(hs2 != 0 && b == "block"){as.vector(B)}else{rep(1, q)})
        E <- mvrnorm(n = n, mu = rep(0, q), Sigma = sigma)
    }

} else if (PTgen == 'simplex'){
    
    tbl <- read.table(sprintf("%s/qlocstdev.%s.tsv", opt$fx, "norm"), h = T)
    colnames(tbl) <- c("Q", "L", "S")
    if(! q %in% unique(tbl$Q)){
        stop(sprintf("stdev not precomputed for q = %s", q))
    } 
    stdev <- subset(tbl, Q == q & L == ploc)$S
    
    x <- c(ploc, rep(1, q-1))
    y0 <- x/sum(x)
    E <- sim.simplex(q, n, y0, stdev, dist = "norm")
    # [NOTE] We are not checking here as intended for q <= 10
  
} else if (PTgen == 'multinom'){
  
    x <- c(ploc, rep(1, q-1))
    p <- x/sum(x)
    N <- 1000
    E <- t(rmultinom(n, N, p))

} else if ( grepl('t-\\d', PTgen) ) {
   
    sigma <- getCov(q, varE, corE, vEr, if(hs2 != 0 && b == "block"){as.vector(B)}else{rep(1, q)})
    E <- rmvt(n, rep(0, q), sigma, df = as.numeric(unlist(strsplit(PTgen, split = "-"))[2]))
    
} else {

    sigma <- getCov(q, varE, corE, vEr, if(hs2 != 0 && b == "block"){as.vector(B)}else{rep(1, q)})
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
