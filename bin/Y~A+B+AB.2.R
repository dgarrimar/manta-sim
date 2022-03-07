#!/usr/bin/env Rscript

## Evaluation of asymptotic PERMANOVA in complex models (II)
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
    make_option(c("-w","--which"), type="character", default = "B",
                help="Which factor changes: A, B or AB", metavar = "character"),
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
    make_option(c("-k","--k_chunk"), type="numeric", default=0,
                help="Chunk number [default %default]", metavar="numeric"),
    make_option(c("-x","--cpu"), type="numeric", default=10,
                help="Number of cores [default %default]", metavar="numeric"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output file name", metavar="character"),
    make_option(c("-f", "--fx"), type="character", default=NULL,
                help="Path to helper functions", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$output)) {
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
dd <- opt$DistDef
lambda <- opt$lambda
stdev <- opt$stdev
loc <- opt$position
pdist <- opt$p_dist
hk <- opt$heterosk
delta <- opt$delta
transf <- opt$transf
S <- round(opt$simulations)
modelSim <- opt$model
w <- opt$which
output <- opt$output
chunk <- opt$k_chunk
cores <- opt$cpu
fx <- opt$fx

## 1. Load packages and functions

library(CompQuadForm)
library(car)
library(MASS)
library(copula)
library(vegan)

source(sprintf("%s/fx.R", opt$fx))

## 2. Define parameters 

labs <- label(a, b, n, u, w)
A <- labs[[1]]
B <- labs[[2]]

if (w == "A") {
   ch <- A
} else if (w == "B") {
    ch <- B
} else if (w == "AB") {
    ch <- A:B
} else {
    stop(sprintf("Unknown factor: '%s'.", w))
}

if (modelSim == "simplex") {
  
    tbl <- read.table(sprintf("%s/qlocstdev.%s.tsv", fx, pdist), h = T)
    colnames(tbl) <- c("Q", "L", "S")
  
    if (! q %in% unique(tbl$Q)) {
        stop(sprintf("stdev not precomputed for q = %s", q))
    } 
  
    stdev <- subset(tbl, Q == q & L == loc)$S
  
}

## 3. Simulate

if (chunk != 0) {
  
    tstats <- c()
    set.seed(chunk)

    for (i in 1:S) {
    
        if (modelSim == "simplex") {
      
            if (i == 1) { # Sanity check
                tol <- 0.01
                check <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = T, dist = pdist)
                wm <- which.max(check$exp)
                if (abs(check[wm, "obs"] - check[wm, "exp"]) > tol) {
                    stop("Deviation from expected centroid greater than tolerance.")
                }
            }
      
            Y <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, dist = pdist)
            sd <- mean(apply(Y, 2, sd))
      
        } else if (modelSim == "mvnorm") {
      
            Y <- Sim.mvnorm(ch, q, n, mu = rep(0, q), delta, hk, Var, Cor)
      
        } else if (modelSim == "multinom") {
      
            # N <- rpois(1, lambda)
            Y <- Sim.multinom(ch, q, n, lambda, delta, loc)
      
        } else if (modelSim == "copula") {
      
            Y <- Sim.copula(ch, q, n, mu = rep(0, q), delta, hk, Var, Cor, dd)
      
        } else {
            stop(sprintf("Unknown option: modelSim = '%s'.", modelSim))
        }
    
        if (transf == "sqrt") {
            Y <- sqrt(Y)
        } else if (transf == "log") {
            Y <- log(Y+1)
        } else if (transf == "none") {
          
        } else {
            stop(sprintf("Unknown option: transf = '%s'.", transf))
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
        # We divide by SSe for numerical stability when running CompQuadForm::davies
        Fs <- f/Df*df.e 
        
        # Asymptotic test statistic (McArtor) 
        G <- crossprod(Y)
        eG <- eigen(G/n, only.values = T)$values # n or df.e ?
    
        # Asymptotic test statistic (Our proposal)
        e <- eigen(cov(R)*(n-1)/df.e, symmetric = T, only.values = T)$values
        
        tstats <- rbind(tstats, c(Fs, unlist(SS), SSe, e, eG))

    }
  
    write.table(tstats, file = output, sep = "\t", col.names = F, row.names = F, quote = F)
  
} else {
  
    set.seed(chunk)
    
    ######### This should be identical to the same piece of code above (except the check)
    if (modelSim == "simplex") {
    
        tol <- 0.01
        check <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = T, dist = pdist)
        wm <- which.max(check$exp)
        if (abs(check[wm, "obs"] - check[wm, "exp"]) > tol) {
            stop("Deviation from expected centroid greater than tolerance.")
        }
    
        Y <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, dist = pdist)
        sd <- mean(apply(Y, 2, sd))
    
    } else if (modelSim == "mvnorm") {
    
        Y <- Sim.mvnorm(ch, q, n, mu = rep(0, q), delta, hk, Var, Cor)
    
    } else if (modelSim == "multinom") {
    
        # N <- rpois(1, lambda)
        Y <- Sim.multinom(ch, q, n, lambda, delta, loc)
    
    } else if (modelSim == "copula") {
    
        Y <- Sim.copula(ch, q, n, mu = rep(0, q), delta, hk, Var, Cor, dd)
    
    } else {
        stop(sprintf("Unknown option: modelSim = '%s'.", modelSim))
    }
  
    if (transf == "sqrt") {
        Y <- sqrt(Y)
    } else if (transf == "log") {
        Y <- log(Y+1)
    } else if (transf == "none") {
    
    } else {
      stop(sprintf("Unknown option: transf = '%s'.", transf))
    }
  
    Y <- scale(Y, center = T, scale = F)
    
    #########
    
    d <- as.dist(interDist(Y))
  
    ## Adonis permutations (raw data + type I SS)
    # ado <- adonis(d ~ A + B + A:B, permutations = S, parallel = cores)

    ## Adonis permutations (strata + type I SS)
    perm <- how(nperm = S)
    setBlocks(perm) <- ch
    ado <- adonis(d ~ A + B + A:B, permutations = perm, parallel = cores) # WARNING: Adonis uses Type I SS, while we use Type II SS!
                                                                          # Only valid to study A:B under H0
    ## Adonis permutations (strata + type II SS)
    # ado <- list()
    # ado$f.perms <- matrix(NA, nrow = 1, ncol = 7)
  
    # Yp <- Y
    # for (p in 1:S) {
    #     print(p)
    #     for (i in 1:b) {Yp[B == i, ] <- Y[sample(which(B == i)),]}
    #     fit <- lm(Yp ~ A + B + A:B)
    #     R <- fit$residuals
    #     UU <- Anova(fit, type = "II") # SS type II
    #     SS <- lapply(UU$SSP, function(x){sum(diag(x))})
    #     SSe <- sum(diag(UU$SSPE))
    #     df.e <- fit$df.residual
    #     Df <- table(fit$assign)[-1]
    #     f <- unlist(lapply(SS, function(x){x/SSe}))
    #     ado$f.perms <- rbind(ado$f.perms, cbind(t(f/Df*df.e), t(as.numeric(SS)), t(SSe)))
    # }
    # ado$f.perms <- ado$f.perms[-1,]
    
    save(ado, file = sprintf("%s.RData",output))
}
