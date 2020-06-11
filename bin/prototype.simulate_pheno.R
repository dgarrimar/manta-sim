
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
 
 n <- 1000              # Number of individuals
 q <- 3                 # Number of responses
 p <- 1000              # Number of variants

 hs2 <- 0.05            # Average fraction of variance explained by causal variants across traits
 hg2 <- 0.3             # Average fraction of variance explained by relatedness across traits
 alphaG <- 0.5          # Fraction of signal from the relatedness contribution that is shared across traits
 lambda <- 0.6          # Fraction of structured noise
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
  