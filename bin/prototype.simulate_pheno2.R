
## 0.1 Libraries and functions

setwd("/nfs/users2/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/real_data/simreal-nf")

library(MASS)

rmatnorm <- function(M, U, V, tol = 1e-12){
  # THIS FUNCTION MAY NOT WORK WHEN U/V is very rank deficient
  # M = mean of matrix
  # U = covariance matix of the columns
  # V = covariance matrix of the rows
  a <- nrow(M)
  b <- ncol(M)
  # Draw Z from MN(O, I, I)
  Z <- matrix(rnorm(a*b,0,1), a, b)
  # Cholesky decomposition of U and V (pivot to allow semidefinite matrices)
  # https://en.wikipedia.org/wiki/Matrix_normal_distribution#Drawing_values_from_the_distribution
  L1 <- suppressWarnings(chol(U, pivot = T))
  L1 <- L1[, order(attr(L1, "pivot"))]
  L2 <- suppressWarnings(chol(V, pivot = T))
  L2 <- L2[, order(attr(L2, "pivot"))]
  # Return draw from MN(M,U,V)
  return(M + crossprod(L1,Z) %*% L2)
}

## 0.2 Parameters

set.seed(12345)

n <- 1e3               # Number of individuals
q <- 10                # Number of responses
p <- 1e4               # Number of variants

hs2 <- 0
hg2 <- 0.5             
alphaG <- 1
lambda <- 0

M <- matrix(0, n, q)   # Means

## 1. Effect of the SNP

# Generate SNP matrix
S <- matrix(NA, n, p)  
mafs <- runif(p, min = 0, max = 0.5)
for (snp in 1:p){
  maf <- mafs[snp]
  S[ ,snp] <- sample(size = n, c(2,1,0), replace = T, prob = c(maf^2, 2*maf*(1-maf), (1-maf)^2 ) )
}  

## 2. Relatedness signal 

write.table(matrix(0, n, p), file = "GEMMA/dummy.tsv", col.names = F, row.names = F, sep = "\t")
GT <- cbind(
  paste0("rs", 1:p),
  matrix(sample(c('A','C','T','G'), replace = T, size = 2*p), p, 2),
  t(S))
write.table(GT, file = "GEMMA/GT.tsv", col.names = F, row.names = F, quote = F, sep = ",")
write.table(GT[,-c(1:3)], file = "GEMMA/GT2.tsv", col.names = F, row.names = F, quote = F, sep = " ")

# Compute Kinship by GEMMA
system('source ~/.conda4R && cd GEMMA && conda activate gemma && gemma -g GT.tsv -p dummy.tsv -gk 2 -outdir . -o Rg')

Rg <- as.matrix(read.table("GEMMA/Rg.sXX.txt"))

## 2.3. Independent
set.seed(2)
aaT <- tcrossprod(matrix(rnorm(q*q), q, q,))
Gs <- rmatnorm(M = M, U = Rg, V = aaT) 
Gs <- Gs / sqrt(mean(diag(cov(Gs)))) * sqrt( hg2*(alphaG ) )# Rescale

E <- matrix(rnorm(n*q), n, q)
E <- E / sqrt(mean(diag(cov(E)))) * sqrt( (1-lambda)*(1-hs2-hg2) )

## 4. Build Y

Y <- Gs + E

## 5. Save

write.table(Y, file = "GEMMA/Y.tsv", col.names = F, row.names = F, sep = "\t")
write.table(t(Y), file = "GEMMA/YT.tsv", col.names = F, row.names = F, sep = "\t")

## 6. Run GEMMA

system('source ~/.conda4R && cd GEMMA && conda activate gemma && gemma -g GT.tsv -p Y.tsv -n 1 2 3 4 5 6 7 8 9 10 -k Rg.sXX.txt -lmm -outdir . -o res')


