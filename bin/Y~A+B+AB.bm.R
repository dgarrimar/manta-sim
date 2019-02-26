#!/usr/bin/Rscript

## Evaluation of asymptotic Anderson's test in complex models
## Model: Y ~ A + B + AB, balanced

##  0. Parse arguments

library(optparse)


option_list = list(
    make_option(c("-q", "--q"), type="numeric", default=3,
                help="Number of dependent variables [default %default]", metavar="numeric"),
    make_option(c("-a","--a_levels"), type="numeric", default=3,
                help="Number of levels of factor A [default %default]", metavar="numeric"),
    make_option(c("-b","--b_levels"), type="numeric", default=3,
                help="Number of levels of factor B [default %default]", metavar="numeric"),
    make_option(c("-r","--replicates"), type="numeric", default=10,
                help="Number of replicates [default %default]", metavar="numeric"),
    make_option(c("-s","--simulations"), type="numeric", default=1e6,
                help="Number of simulations [default %default]", metavar="numeric"),
    make_option(c("-k","--k_chunk"), type="numeric", default=0,
                help="Chunk number [default %default]", metavar="numeric"),
    make_option(c("-c","--cores"), type="numeric", default=10,
                help="Number of cores [default %default]", metavar="numeric"),
    make_option(c("-m", "--model"), type="character", default="mvnorm",
                help="H0-generator model: 'dirichlet' or 'mvnorm' [default %default]",
                metavar="character"),
    make_option(c("-d","--delta"), type="numeric", default=0.1,
                help="H1 generation parameter [default %default]", metavar="numeric"),
    make_option(c("-t","--transf"), type="character", default="None",
                help="Data transformation: sqrt or None [default %default]", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$output)){
    print_help(opt_parser)
    stop("An output file must be supplied\n", call.=FALSE)
}

q <- opt$q
a <- opt$a_levels
b <- opt$b_levels
r <- opt$replicates
delta <- opt$delta
transf <- opt$transf
S <- round(opt$simulations)
modelSim <- opt$model
chunk <- opt$k_chunk
cores <- opt$cores
output <- opt$output

## 1. Load packages and functions

library(CompQuadForm)
library(car)
library(MCMCpack)
library(MASS)
library(vegan)

source("/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/nf/bin/fx.R")

## 2. Define parameters 

n <- a*b*r                            # Total n
A <- gl(a, b*r)                       # df = a-1
B <- gl(b, r, length = a*b*r)         # df = b-1

                                      # AB df = (a-1)*(b-1)
                                      # Error df = a*b*(r-1) = n - a*b

## 3. Simulate

if (chunk != 0) {
  
  tstats <- c()
  set.seed(chunk)
  
  for (i in 1:S){
    
    if (modelSim == "dirichlet") {
  
      h1 <- c(1, rep(0, q-1))
      
      alpha <- rep(1/q, q)
      Y <- rdirichlet(n, alpha)
  
      if (delta != 0){    
        s2d <- step2distance(alpha, h1, delta)
        if(is.na(s2d$d)){stop("delta => out of the simplex")}
        Y[B == 1,] <- rdirichlet(nrow(Y[B == 1,]), s2d$r)
      }
    } else if (modelSim == "mvnorm") {
      
      Y <- mvrnorm(n, mu = rep(1/q, q) , Sigma = diag(rep(1,q)))
      
      if (delta != 0){
        Y[B == 1, 1] <- Y[B == 1, 1] + delta
      }
    } else {
      NA
    }
    
    if (transf == "sqrt"){
      Y <- sqrt(Y)
    } 
    
    # lm and residuals
    Y <- scale(Y, center = T, scale = F)
    fit <- lm(Y ~ A + B + A:B)
    R <- fit$residuals
    
    # Sums of squares
    UU <- Anova(fit, type = "II") # SS type II 
    SS <- lapply(UU$SSP, function(x){sum(diag(x))})
    SSe <- sum(diag(UU$SSPE))
  
    # Df
    df.e <- fit$df.residual       # df.e <- (n-1) - sum(Df)
    Df <- table(fit$assign)[-1]
    names(Df) <- attributes(fit$terms)$term.labels
    
    # Statistic 
    f <- unlist(lapply(SS, function(x){x/SSe})) 
    Fs <- f/Df*df.e
    
    # Asymptotic test statistic (McArtor) 
    G <- crossprod(Y)
    eG <- eigen(G/n, only.values = T)$values # n or df.e ?
    
    # Asymptotic test statistic (Our proposal)
    e <- eigen(cov(R)*(n-1)/df.e, symmetric = T, only.values = T)$values
    
    tstats <- rbind(tstats, c(Fs, unlist(SS)))
    
  }
  
  write.table(tstats, file = output, sep = "\t", col.names = F, row.names = F, quote = F)
  
} else {

  set.seed(chunk)
  
  ######### This should be identical to the same piece of code above
  if (modelSim == "dirichlet") {
    
    h1 <- c(1, rep(0, q-1))
    
    alpha <- rep(1/q, q)
    Y <- rdirichlet(n, alpha)
    
    if (delta != 0){    
      s2d <- step2distance(alpha, h1, delta)
      if(is.na(s2d$d)){stop("delta => out of the simplex")}
      Y[B == 1,] <- rdirichlet(nrow(Y[B == 1,]), s2d$r)
    }
  } else if (modelSim == "mvnorm") {
    
    Y <- mvrnorm(n, mu = rep(1/q, q) , Sigma = diag(rep(1,q)))
    
    if (delta != 0){
      Y[B == 1, 1] <- Y[B == 1, 1] + delta
    }
  } else {
    NA
  }
  
  if (transf == "sqrt"){
    Y <- sqrt(Y)
  } 
  
  Y <- scale(Y, center = T, scale = F)
  #########
  
  d <- as.dist(interDist(Y))
  
  # Adonis permutations (raw data)
  ado <- adonis(d ~ A + B + A:B, permutations = S, parallel = cores)

  save(ado, file = sprintf("%s.RData",output))
}









