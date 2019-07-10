#!/usr/bin/Rscript

## Evaluation of asymptotic Anderson's test in complex models
## Model: Y ~ A + C

##  0. Parse arguments

library(optparse)


option_list = list(
  make_option(c("-q", "--q"), type="numeric", default=3,
              help="Number of dependent variables [default %default]", metavar="numeric"),
  make_option(c("-a","--a_levels"), type="numeric", default=2,
              help="Number of levels of factor A [default %default]", metavar="numeric"),
  make_option(c("-b","--b_levels"), type="numeric", default=3,
              help="Number of levels of factor B [default %default]", metavar="numeric"),
  make_option(c("-r", "--C_noise"), type="numeric", default=0,
              help="(If generator = simplex: noise added to the) correlation of covariate C with Y[,1] [default %default]", metavar="numeric"),
  make_option(c("--C_mean"), type="numeric", default=0,
              help="Mean of covariate C [default %default]", metavar="numeric"),
  make_option(c("--C_var"), type="numeric", default=1,
              help="Variance of covariate C [default %default]", metavar="numeric"),
  make_option(c("-u","--unbalance"), type="numeric", default=1,
              help="Unbalance degree (level 1) [default %default]", metavar="numeric"),
  make_option(c("-n","--no_samples"), type="numeric", default=100,
              help="Total number of samples [default %default]", metavar="numeric"),
  make_option(c("-S","--simulations"), type="numeric", default=1e3,
              help="Number of simulations [default %default]", metavar="numeric"),
  make_option(c("-m", "--model"), type="character", default="mvnorm",
              help="H0/H1 generator model: 'mvnorm', 'simplex' or 'multinom' [default %default]",
              metavar="character"),
  make_option(c("-d","--delta"), type="numeric", default=0,
              help="H1 generation parameter [default %default]", metavar="numeric"),
  make_option(c("-c","--correlation_y"), type="numeric", default=0,
              help="Correlation of the Y variables [default %default]", metavar="numeric"),
  make_option(c("-w","--which"), type="character", default = "C",
              help="Which factor changes: A or C", metavar = "character"),
  make_option(c("-v","--variance_y"), type="character", default="equal",
              help="Variance of the Y variables: 'equal' or 'unequal' [default %default]", metavar="character"),
  make_option(c("-s","--stdev"), type="numeric", default=0.1,
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
  make_option(c("-t","--transf"), type="character", default="none",
              help="Data transformation: sqrt or None [default %default]", metavar="character"),
  make_option(c("--adonis"), type="numeric", default=0,
              help="Should permutation test be performed? Specify number of permutations [default %default]",
              metavar="numeric"),
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
r <- opt$C_noise
C_mean <- opt$C_mean
C_var <- opt$C_var
n <- opt$no_samples
u <- opt$unbalance
Cor <- opt$correlation_y
Var <- opt$variance_y
dd <- opt$DistDef
lambda <- opt$lambda
stdev <- opt$stdev
loc <- opt$position
pdist <- opt$p_dist
hk <- opt$heterosk
delta <- opt$delta
transf <- opt$transf
S <- opt$simulations
modelSim <- opt$model
output <- opt$output
adonis <- opt$adonis

## 1. Load packages and functions

library(CompQuadForm)
library(car)
library(MCMCpack)
library(MASS)
library(plyr)
library(copula)

source("/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/nf/bin/fx.R")

## 2. Define parameters 
# which here is useless, changes are controlled using delta,r
labs <- label(a, 1, n, u, "A")
A <- labs[[1]]

if(modelSim == "simplex"){
  
  set.seed(1)
  rc <- c()
  
  tbl <- read.table(sprintf("/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/nf/bin/qlocstdev.%s.tsv", pdist), h = T)
  tbl2 <- read.table(sprintf("/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/nf/bin/qlocstdev2.%s.tsv", pdist), h = T)
  colnames(tbl) <- colnames(tbl2) <- c("Q", "L", "S")
  
  if(! q %in% unique(tbl$Q) || ! q %in% unique(tbl2$Q)){
    stop(sprintf("stdev not precomputed for q = %s", q))
  } else if (! loc %in% unique(tbl$L) || ! loc %in% unique(tbl$L)) {
    stop(sprintf("stdev not precomputed for p_loc = %s", loc))
  }
  
  stdev <- subset(tbl, Q == q & L == loc)$S
  stdev2 <- subset(tbl2, Q == q & L == loc)$S
  
  C.gen <- runif(n, min = -stdev2, max = stdev2)
  C.check <- runif(1e4, min = -stdev2, max = stdev2)
  C <- C.gen + runif(n, -r, r)
} else if (modelSim == "multinom"){
  set.seed(1)
  rc <- c()
  tbl2 <- read.table("/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/nf/bin/qlocstdev2.multinom.tsv", h = T)
  colnames(tbl2) <- c("Q", "L", "S")
  if(! q %in% unique(tbl2$Q)){
    stop(sprintf("stdev not precomputed for q = %s", q))
  } else if (! loc %in% unique(tbl2$L)) {
    stop(sprintf("stdev not precomputed for p_loc = %s", loc))
  }
  stdev2 <- subset(tbl2, Q == q & L == loc)$S
  C.gen <- runif(n, min = -stdev2, max = stdev2)
  C <- C.gen + runif(n, -r, r)
}

## 3. Simulate
pv.mt <- c()
for (i in 1:S){
  
  set.seed(i)
  
  if (modelSim == "simplex") {

    if(i == 1){ # Sanity check
      tol <- 0.01
      check <- Sim.simplex.C(A, C.check, q, n, loc, delta, hk, stdev, check = T, dist = pdist)
      wm <- which.max(check$exp)
      if( abs(check[wm, "obs"] - check[wm, "exp"]) > tol ) {
        stop("Deviation from expected centroid greater than tolerance.")
      }
    }
    
    Y <- Sim.simplex.C(A, C.gen, q, n, loc, delta, hk, stdev, dist = pdist)
    sd <- mean(apply(Y, 2, sd))
    rc <- c(rc, cor(Y[,1], C))
    
  } else if (modelSim == "mvnorm") {
    
    Yext <- Sim.mvnorm.C(A, q, n, mu = rep(0, q), r, delta, hk, Var, Cor, C_mean, C_var)
    Y <- Yext[, -1]
    C <- Yext[, 1]
    
  } else if (modelSim == "multinom") {
   
    N <- rpois(1, lambda)
    Y <- Sim.multinom.C(A, C, q, n, N, delta, loc)
    rc <- c(rc, cor(Y[,1], C))
    
  } else if (modelSim == "copula") {

    Yext <- Sim.copula.C(A, q, n, mu = rep(0, q), r, delta, hk, Var, Cor, dd, C_mean, C_var)
    Y <- Yext[, -1]
    C <- Yext[, 1]

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
  
  if (adonis != 0){
    library(vegan)
    ADONIS <- adonis(dist(Y) ~ A + B + A:B, permutations = adonis)$aov.tab[,6][1:3]
  }  
  
  # lm and residuals
  Y <- scale(Y, center = T, scale = F)
  fit <- lm(Y ~ C + A)
  R <- fit$residuals
  
  # Sums of squares
  UU <- Anova(fit, type = "II") # SS type II 
  SS <- lapply(UU$SSP, function(x){sum(diag(x))})
  SSe <- sum(diag(UU$SSPE))
  
  # Statistic 
  f <- unlist(lapply(SS, function(x){x/SSe})) 
  # We divide by SSe for numerical stability when running CompQuadForm::davies
  # Note that F = f/Df*df.e
  
  # Df
  df.e <- fit$df.residual       # df.e <- (n-1) - sum(Df)
  Df <- table(fit$assign)[-1]
  names(Df) <- attributes(fit$terms)$term.labels
  
  # eigendecomposition and pv calculation
  e <- eigen(cov(R)*(n-1)/df.e, symmetric = T, only.values = T)$values
  pv.acc <- mapply(pv.f, f = f, df.i = Df, MoreArgs = list(df.e = df.e, lambda = e))
  MANOVA <- tryCatch({summary(manova(fit))$stats[,6][1:2]}, 
                     error = function(e){return(rep(NA, 3))}) # MANOVA added for comparison, it may fail with lin. dep. variables
  pv.mt <- rbind(pv.mt, c(pv.acc[1,], MANOVA)) 
  
  if (adonis==0){
    pv.mt <- rbind(pv.mt, c(pv.acc[1,], MANOVA)) 
  } else {
    names(ADONIS) <- names(pv.acc[1,])
    pv.mt <- rbind(pv.mt, c(pv.acc[1,], MANOVA, ADONIS))
  }
}

if(modelSim == "mvnorm"){
  params <- c(a, C_mean, C_var, n, u, q, delta, r, hk, Var, Cor, transf)
} else if(modelSim == "copula"){
  params <- c(a, C_mean, C_var, n, u, q, delta, r, hk, Var, Cor, dd, transf)
} else if(modelSim == "simplex"){
  params <- c(a, C_mean, C_var, n, u, q, delta, mean(rc), hk, loc, pdist, sd, transf)
} else if(modelSim == "multinom"){
  params <- c(a, C_mean, C_var, n, u, q, delta, mean(rc), hk, loc, lambda, transf)
}

if (S == 1 && adonis != 0){
  result2write <- pv.mt
} else {
  result2write <- colMeans(pv.mt < 0.05)
}

write.table(t(c(params, result2write)), file = output, sep = "\t", col.names = F, row.names = F, quote = F)


