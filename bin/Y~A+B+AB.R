#!/usr/bin/Rscript

## Evaluation of asymptotic Anderson's test in complex models
## Model: Y ~ A + B + AB

##  0. Parse arguments

library(optparse)


option_list = list(
  make_option(c("-q", "--q"), type="numeric", default=3,
              help="Number of dependent variables [default %default]", metavar="numeric"),
  make_option(c("-a","--a_levels"), type="numeric", default=2,
              help="Number of levels of factor A [default %default]", metavar="numeric"),
  make_option(c("-b","--b_levels"), type="numeric", default=3,
              help="Number of levels of factor B [default %default]", metavar="numeric"),
  make_option(c("-u","--unbalance"), type="numeric", default=1,
              help="Unbalance degree (B level 1) [default %default]", metavar="numeric"),
  make_option(c("-n","--no_samples"), type="numeric", default=100,
              help="Total number of samples [default %default]", metavar="numeric"),
  make_option(c("-S","--simulations"), type="numeric", default=1e6,
              help="Number of simulations [default %default]", metavar="numeric"),
  make_option(c("-m", "--model"), type="character", default="mvnorm",
              help="H0-generator model: 'mvnorm' or 'simplex' [default %default]",
              metavar="character"),
  make_option(c("-d","--delta"), type="numeric", default=0,
              help="H1 generation parameter [default %default]", metavar="numeric"),
  make_option(c("-c","--correlation_y"), type="numeric", default=0,
              help="Correlation of the Y variables [default %default]", metavar="numeric"),
  make_option(c("-v","--variance_y"), type="character", default="equal",
              help="Variance of the Y variables: 'equal' or 'unequal' [default %default]", metavar="character"),
  make_option(c("-s","--stdev"), type="numeric", default=0.1,
              help="stdev for the 'simplex' generator model [default %default]", metavar="numeric"),
  make_option(c("-l","--location"), type="numeric", default=1,
              help="location of the 'simplex' generator model [default %default]", 
              metavar="numeric"),
  make_option(c("-H","--heterosk"), type="numeric", default= 1,
              help="Heteoskedasticity degree (B level 1) [default %default]", metavar="numeric"),
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
n <- opt$no_samples
u <- opt$unbalance
Cor <- opt$correlation_y
Var <- opt$variance_y
stdev <- opt$stdev
loc <- opt$location
hk <- opt$heterosk
delta <- opt$delta
transf <- opt$transf
S <- opt$simulations
modelSim <- opt$model
output <- opt$output

## 1. Load packages and functions

library(CompQuadForm)
library(car)
library(MCMCpack)
library(MASS)

source("/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simul/nf/bin/fx.R")

## 2. Define parameters 

labs <- label(a, b, n, u)
A <- labs[[1]]
B <- labs[[2]]

## 3. Simulate

pv.mt <- c()
for (i in 1:S){
  
  set.seed(i)
  
  if (modelSim == "simplex") {
    
    Y <- Sim.simplex(B, q, n, loc, delta, hk, stdev)
   
  } else if (modelSim == "mvnorm") {
    
    Y <- Sim.mvnorm(B, q, n, mu = rep(0, q), delta, hk, Var, Cor)
    
  } else {
    stop(sprintf("Unknown option: modelSim = '%s'.", modelSim))
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
  pv.mt <- rbind(pv.mt, c(pv.acc[1,],summary(manova(fit))$stats[,6][1:3])) # MANOVA added for comparison
}

if(modelSim == "mvnorm"){
  params <- c(a, b, n, u, q, delta, hk, Var, Cor, transf)
} else if(modelSim != "simplex"){
  params <- c(a, b, n, u, q, delta, hk, loc, stdev, transf)
}

result2write <- colMeans(pv.mt < 0.05)

write.table(t(c(params, result2write)), file = output, sep = "\t", col.names = F, row.names = F, quote = F)


