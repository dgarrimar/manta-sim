#!/usr/bin/env Rscript

## Evaluation of asymptotic Anderson's test in complex models (running time)
## Model: Y ~ A + B + AB

##  0. Parse arguments

library(optparse)

option_list = list(
  make_option(c("-n", "--n"), type="numeric", default=100,
              help="Number of samples [default %default]", metavar="numeric"),
  make_option(c("-q", "--q"), type="numeric", default=3,
              help="Number of dependent variables [default %default]", metavar="numeric"),
  make_option(c("-r", "--replicates"), type="numeric", default=1,
              help="Number of replicates [default %default]", metavar="numeric"),
  make_option(c("-P", "--permutations"), type="numeric", default=1000,
              help="Number of permutations [default %default]", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("An output file must be supplied\n", call.=FALSE)
}

n <- opt$n
q <- opt$q
r <- opt$replicates
P <- opt$permutations
output <- opt$output

## 1. Load packages and functions

library(mlm)
library(vegan)
library(MASS)

source("/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/nf/bin/fx.R")

## 2. Define parameters 

labs <- label(a = 2, b = 3, n, u = 1, w = "B")
A <- labs[[1]]
B <- labs[[2]]

rt.df <- c()

## 3. Simulate
 
set.seed(opt$r)  
Y <- Sim.mvnorm(B, q, n, mu = rep(0, q), delta = 1, hk = 1, Var = 'equal', Cor = 0)

## 4. Measure time 
 
 # MLM
 t0_mlm <- Sys.time()
  MLM <- mlm(Y ~ A + B + A:B, type = "I")
 t1_mlm <- Sys.time() 

 # ADONIS 
 t0_adonis <- Sys.time()
  perm <- how(nperm = P)
  setBlocks(perm) <- B
  ADONIS <- adonis(dist(Y) ~ A + B + A:B, permutations = perm, parallel = 1)
 t1_adonis <- Sys.time()

 time <- data.frame(r, n, q, P,
                    mlm = as.numeric(difftime(t1_mlm, t0_mlm, units = "secs")),
                    adonis = as.numeric(difftime(t1_adonis, t0_adonis, units = "secs")))

## 4. Save output
 
 write.table(time, file = output, sep = "\t", col.names = F, row.names = F, quote = F)


