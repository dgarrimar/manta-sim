#!/usr/bin/env Rscript

## Simulate phenotypes

##  0. Parse arguments

library(optparse)

option_list = list(
  make_option(c("-s","--seed"), type="numeric", default=0,
              help="Set seed for random processes [default %default]", metavar="numeric"),
  make_option(c("-n","--nb_samples"), type="numeric", default=1000,
              help="Total number of samples [default %default]", metavar="numeric"),
  make_option(c("-q", "--nb_responses"), type="numeric", default=3,
              help="Number of dependent variables [default %default]", metavar="numeric"),
  make_option(c("-m", "--model"), type="character", default="mvnorm",
              help="H0/H1 generator model: 'mvnorm', 'simplex', 'multinom' or 'copula' [default %default]",
              metavar="character"),
  make_option(c("-d","--delta"), type="numeric", default=0,
              help="H1 generation parameter [default %default]", metavar="numeric"),
  make_option(c("-c","--correlation_y"), type="numeric", default=0,
              help="Correlation of the Y variables [default %default]", metavar="numeric"),
  make_option(c("-v","--variance_y"), type="character", default="equal",
              help="Variance of the Y variables: 'equal' or 'unequal' [default %default]", metavar="character"),
  make_option(c("-S","--stdev"), type="numeric", default=0.1,
              help="stdev for the 'simplex' generator model [default %default]", metavar="numeric"),
  make_option(c("-p","--position"), type="numeric", default=1,
              help="location of the 'simplex' generator model [default %default]", 
              metavar="numeric"),
  make_option(c("--p_dist"), type="character", default="norm",
              help="Distribution for 'simplex' generator model: norm, gamma or beta [default %default]", 
              metavar="character"),
  make_option(c("-D","--DistDef"), type="character", default="unif-0-1",
              help="Multivariate non-normal distribution definition [default %default]", 
              metavar="character"),
  make_option(c("-l","--lambda"), type="numeric", default=1000,
              help="lambda parameter (Poisson distribution) to generate size for 'multinom' generator model [default %default]", 
              metavar="numeric"),
  make_option(c("-H","--heterosk"), type="numeric", default= 1,
              help="Heteoskedasticity degree (level 1) [default %default]", metavar="numeric"),
  make_option(c("--geno"), type="character", default=NULL,
              help="genotype data (single variant) obtained by simulateGT.nf", metavar="character"),
  make_option(c("-t","--transf"), type="character", default="none",
              help="Data transformation: sqrt or None [default %default]", metavar="character"),
  make_option(c("-f", "--fx"), type="character", default=NULL,
              help="Path to helper functions", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output (simulated phenotype) file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("An output file must be supplied\n", call.=FALSE)
}

set.seed(opt$seed)

n <- opt$nb_samples
q <- opt$nb_responses                
Cor <- opt$correlation_y
Var <- opt$variance_y
dd <- opt$DistDef
lambda <- opt$lambda
stdev <- opt$stdev
loc <- opt$position
pdist <- opt$p_dist
hk <- opt$heterosk
delta <- opt$delta
modelSim <- opt$model
geno <- opt$geno
transf <- opt$transf
fx <- opt$fx
outfile <- opt$output

## 1. Load packages and functions

library(optparse)
library(MASS)
library(data.table)
library(BEDMatrix)
library(copula)

source(sprintf("%s/fx.R", opt$fx))

## 2. Define parameters 

if(modelSim == "simplex"){
  
  tbl <- read.table(sprintf("%s/qlocstdev.%s.tsv", fx, pdist), h = T)
  colnames(tbl) <- c("Q", "L", "S")
  
  if(! q %in% unique(tbl$Q)){
    stop(sprintf("stdev not precomputed for q = %s", q))
  } 
  
  stdev <- subset(tbl, Q == q & L == loc)$S
  
}

set.seed(opt$seed)

# Real SNP

ch <- as.matrix(BEDMatrix(geno, simple_names = T))[,1]

# Apply same filter as in GEMMA regarding MAF (default)
# Recode alleles

tbl <- table(ch)
if (length(tbl) == 1){
    maf <- 0
} else if (length(tbl) == 2){
    gt <- 0:2
    w <- gt[!gt %in% names(tbl)]
    if(w == 0){
        tbl <- c(0, tbl)
    } else if (w == 1){
        tbl <- c(tbl[1], 0, tbl[2])
    } else if (w == 2) {
        tbl <- c(tbl, 0)
    }
    names(tbl)[w+1] <- w
    maf <- (tbl['2'] + tbl['1']/2) / sum(tbl)
    maf <- min(maf, 1-maf)
    # Relabel
    ma <- names(which.min(table(ch))) # minor gt OBSERVED (could be het.)
    ch[ch == ma] <- 3 # Minor (mask)
    ch[ch != 3] <- 2  # Major
    ch[ch == 3] <- 1  # Minor (unmask)
} else {
  ma <- names(which.min(tbl[c('0', '2')])) # minor gt 
  maf <- (tbl[ma] + tbl['1']/2) / sum(tbl)
  Ma <- ifelse(ma == '0', '2', '0')
  ch[ch == 1] <- 3  # Het
  ch[ch == ma] <- 1 # Minor
  ch[ch == Ma] <- 2 # Major
} 

if(maf < 0.01) {

  write.table(matrix(0, n, q), file = outfile, col.names = F, row.names = F, sep = "\t")

} else{
  
  ## 3. Simulate
  
  ch <- as.factor(ch) 
  
  if (modelSim == "simplex") {
    
    Y <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, dist = pdist) 
    # [NOTE] We are not checking here as intended for q <= 10
    
  } else if (modelSim == "mvnorm") {
    
    Y <- Sim.mvnorm(ch, q, n, mu = rep(0, q), delta, hk, 1, Var, Cor)
    
  } else if (modelSim == "multinom") {
    
    N <- rpois(1, lambda)
    Y <- Sim.multinom(ch, q, n, N, delta, loc)
    
  } else if (modelSim == "copula") {
    
    Y <- Sim.copula(ch, q, n, mu = rep(0, q), delta, hk, Var, Cor, dd)
    
  } else {
    stop(sprintf("Unknown option: modelSim = '%s'.", modelSim))
  }
 
  if (transf == "sqrt"){
    Y <- sqrt(Y)
  } else if (transf == "log"){
    Y <- log(Y+1)
  } else if (transf == "none"){
    
  } else {
    stop(sprintf("Unknown option: transf = '%s'.", transf))
  }  

  ## 4. Save    
  write.table(Y, file = outfile, col.names = F, row.names = F, sep = "\t")
}

