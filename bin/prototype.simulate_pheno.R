
## 0.1 Libraries and functions

 setwd("/nfs/users2/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/real_data/simreal-nf")

 library(MASS)

 rmatnorm2 <- function(M, U, V){
  # Sample from Matrix Normal distribution via Kronecker 
  n <- nrow(M)
  q <- ncol(M)
  X <- matrix(mvrnorm(mu = c(M), Sigma = kronecker(V, U)), n, q)
  return(X)
}

## 0.2 Parameters

 set.seed(12345)
 
 n <- 1e3               # Number of individuals
 q <- 3                 # Number of responses
 p <- 1e4               # Number of variants

 hs2 <- 0.05            # Average fraction of variance explained by causal variants across traits
 hg2 <- 0.8             # Average fraction of variance explained by relatedness across traits
 alphaG <- 0.5          # Fraction of signal from the relatedness contribution that is shared across traits
 lambda <- 0.5          # Fraction of structured noise
 alphaH <- 0.5          # Fraction of structured noise that is shared across traits
 k <- 10                # Number of unknown covariates
 
 s <- 5                 # Number of SNPs
 M <- matrix(0, n, q)   # Means
 
## 1. Effect of the SNP

 # Generate SNP matrix
 S <- matrix(NA, n, p)  
 mafs <- runif(p, min = 0, max = 0.5)
 for (snp in 1:p){
   maf <- mafs[snp]
   S[ ,snp] <- sample(size = n, c(2,1,0), replace = T, prob = c(maf^2, 2*maf*(1-maf), (1-maf)^2 ) )
 }  

 # Select SNP(s)
 # mafs[1:s]
 X <- S[, 1:s, drop = F]
 
 # Generate effects
 B <- matrix(sample(c(-1,1), replace = T, size = s*q), s, q)
 
 # Multiply
 XB <- scale(X, center = T, scale = F) %*% B 
 XB <- XB / sqrt(mean(diag(cov(XB)))) * sqrt(hs2)

## 2. Relatedness signal 

 # Rg <- tcrossprod(scale(S, center = T, scale = F))/p 
 write.table(matrix(0, n, p), file = "GEMMA/dummy.tsv", col.names = F, row.names = F, sep = "\t")
 GT <- cbind(
   paste0("rs", 1:p),
   matrix(sample(c('A','C','T','G'), replace = T, size = 2*p), p, 2),
   t(S))
 write.table(GT, file = "GEMMA/GT.tsv", col.names = F, row.names = F, quote = F, sep = ",")
 
 # Compute Kinship by GEMMA
 system('source ~/.conda4R && cd GEMMA && conda activate gemma && gemma -g GT.tsv -p dummy.tsv -gk -outdir . -o Rg')
 
 Rg <- as.matrix(read.table("GEMMA/Rg.cXX.txt"))
 
 ## 2.2. Shared
  aaT <- tcrossprod(rnorm(q))
  Gs <- rmatnorm2(M = M, U = Rg, V = aaT) 
  Gs <- Gs / sqrt(mean(diag(cov(Gs)))) * sqrt( hg2*alphaG ) # Rescale
 
 ## 2.3. Independent
  C <- diag(rnorm(q)^2)
  Gi <- rmatnorm2(M = M, U = Rg, V = C)
  Gi <- Gi / sqrt(mean(diag(cov(Gi)))) * sqrt(hg2*(1-alphaG)) # Rescale
  
## 3. Residuals

 LLT <- tcrossprod(matrix(rnorm(n*k), n, k))
  
 ## 3.1. Structured noise (shared)

  aHaHT <- tcrossprod(rnorm(q))
  Hs <- rmatnorm2(M = M, U = LLT, V = aHaHT)
  Hs <- Hs / sqrt(mean(diag(cov(Hs)))) * sqrt( lambda*alphaH*(1-hs2-hg2) ) # Rescale 
  
 ## 3.2. Structured noise (independent)

  CH <- diag(rnorm(q)^2) 
  Hi <- rmatnorm2(M = M, U = LLT, V = CH)
  Hi <- Hi / sqrt(mean(diag(cov(Hi)))) * sqrt( lambda*(1-alphaH)*(1-hs2-hg2) ) # Rescale

 ## 3.3. Unstructured noise
  D <- 10
  E <- 1
  for (d in 1:D){
    E <- E * matrix(rnorm(n*q), n, q)
  }
  E <- E / sqrt(mean(diag(cov(E)))) * sqrt( (1-lambda)*(1-hs2-hg2) )
 
## 4. Build Y
  
  Y <- XB + (Gs + Gi) + (Hs + Hi) + E
  
## 5. Save
  
  write.table(Y, file = "GEMMA/Y.tsv", col.names = F, row.names = F, sep = "\t")
  
## 6. Run GEMMA
  
  system('source ~/.conda4R && cd GEMMA && conda activate gemma && gemma -g GT.tsv -p Y.tsv -n 1 2 3 -k Rg.cXX.txt -lmm -outdir . -o res')

  res.gemma <- read.table("GEMMA/res.assoc.txt", as.is = T, h = T)[,c(2,17)]
  colnames(res.gemma) <- c("snp", "pv.gemma")

## 7. Run MLM

  library(mlm)

  chol_solve <- function(K, tol = 1e-12) {
    a = eigen(K)$vectors
    b = eigen(K)$values
    b[b < tol] = tol
    b = 1/sqrt(b)
    return(a%*%diag(b)%*%t(a))
  }
  rotate <- function(Y, sigma) {
    U <- chol_solve(sigma)
    tU <- t(U)
    UY <- tU%*%Y
    return(UY)
  }
 
  # Read variance components (median across phenotypes)
  # gemma -vc -p Y.tsv -k Rg.cXX.txt -n 3 -outdir . -o vc
  # gemma -p Y.tsv -g GT.tsv -k Rg.cXX.txt -lmm -n 3 -outdir . -o uni
  # for i in {1..3}; do gemma -vc 2 -p Y.tsv -k Rg.cXX.txt -n $i -outdir . -o vc > /dev/null ; grep -e "sigma2 estimates =" vc.log.txt | cut -d' ' -f 7,9; done
  
  # Transform
  Vg <- hg2 
  Ve <- 1 - hs2 - hg2 
  sigma <- Vg*Rg + Ve*diag(n)
  UY <- rotate(Y, sigma)		# Rotate genotypes and phenotypes
  US <- rotate(S, sigma)
  
  hmy <- 100
  pv.mlm2 <- pv.mlm <- rep(NA, hmy) # p
  for (i in 1:hmy){ # p
    print(i)
    if (length(table(S[,i])) < 3 || min(table(S[,i])) < 5 ) {next}
    pv.mlm[i] <- mlm(Y ~ S[, i])$aov.tab[1,6]
    pv.mlm2[i] <- mlm(UY ~ US[, i])$aov.tab[1,6]
  }

  res.mlm <- data.frame("snp" = paste0("rs",1:hmy), pv.mlm, pv.mlm2)
  res <- merge(res.gemma, res.mlm)
  res <- res[apply(res, 1, function(x){all(!is.na(x))}), ]
  res$snp <- as.numeric(gsub("rs", "", res$snp))
  res <- res[order(res$snp), ]

## 8. Plots
  library(reshape2)
  library(ggplot2)

  colnames(res) <- c("Position","GEMMA", "MLM", "MLM(transf)")
  res$color <- 0
  res$color[1:s] <- 1
  df.melt <- melt(res, id.vars = c("Position", "color"),
                  variable.name = "method", value.name = "pvalue")
  df.melt$FDR <- unlist(tapply(df.melt$pvalue, df.melt$method, function(x){p.adjust(x, "BH")}))

  ggplot(df.melt) +
    geom_bar(aes(x = Position, y = -log10(pvalue),
                 fill = as.factor(color)), stat = 'identity',
             position="identity", width = 1) +
    geom_hline(yintercept = -log10(0.05), col = "blue") +
    facet_wrap(~method, scales = "free") +
    scale_fill_manual(values = c("black", "red")) +
    guides(fill = F) +
    theme_classic(base_size = 22) +
    theme(strip.background = element_blank()) +
    labs(x = "Position", y = expression(-log[10](italic(p))))


  # source("~/bin/qqplot.R")
  # QQplot(as.numeric(res$GEMMA))
  # QQplot(as.numeric(res$MLM))
