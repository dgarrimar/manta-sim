#!/usr/bin/env Rscript

## Evaluation of asymptotic Anderson's test in complex models

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
outfile <- opt$output


## 1. Load packages and functions

library(optparse)
library(MASS)
library(data.table)
library(BEDMatrix)
library(copula)

source("/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/nf/bin/fx.R")

## 2. Define parameters 

if(modelSim == "simplex"){
  
  tbl <- read.table(sprintf("/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/nf/bin/qlocstdev.%s.tsv", pdist), h = T)
  colnames(tbl) <- c("Q", "L", "S")
  
  if(! q %in% unique(tbl$Q)){
    stop(sprintf("stdev not precomputed for q = %s", q))
  } 
  
  stdev <- subset(tbl, Q == q & L == loc)$S
  
}

set.seed(opt$seed)

# Real SNP
ch <- as.factor(as.matrix(BEDMatrix(geno, simple_names = T))[,1])

# Apply same filter as in GEMMA regarding MAF (default)
tbl <- table(ch)
maf <- sum(c(tbl[3], tbl[2]/2), na.rm = T) / sum(tbl, na.rm = T)

if(maf < 0.01) {
  write.table(matrix(0, n, q), file = outfile, col.names = F, row.names = F, sep = "\t")
} else{
  ## 3. Simulate
  
  if (modelSim == "simplex") {
    
    tol <- 0.01
    check <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = T, dist = pdist)
    wm <- which.max(check$exp)
    if( abs(check[wm, "obs"] - check[wm, "exp"]) > tol ) {
      stop("Deviation from expected centroid greater than tolerance.")
    }
    
    Y <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, dist = pdist)
    sd <- mean(apply(Y, 2, sd))
    
  } else if (modelSim == "mvnorm") {
    
    Y <- Sim.mvnorm(ch, q, n, mu = rep(0, q), delta, hk, hk2 = 1, Var, Cor)
    
  } else if (modelSim == "multinom") {
    
    N <- rpois(1, lambda)
    Y <- Sim.multinom(ch, q, n, N, delta, loc)
    
  } else if (modelSim == "copula") {
    
    Y <- Sim.copula(ch, q, n, mu = rep(0, q), delta, hk, Var, Cor, dd)
    
  } else {
    stop(sprintf("Unknown option: modelSim = '%s'.", modelSim))
  }
  
  ## 4. Save    
  write.table(Y, file = outfile, col.names = F, row.names = F, sep = "\t")
}

