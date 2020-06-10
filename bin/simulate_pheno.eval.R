
## 0.1 Libraries and functions

 library(MASS)
 # library(MixMatrix) # R 3.5.2

 set.seed(12345)

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

 rmatnorm2 <- function(M, U, V){
  # TOO SLOW FOR LARGE MATRICES, BUT REQUIRED WHEN U or V is not positive definite
  n <- nrow(M)
  q <- ncol(M)
  X <- matrix(mvrnorm(mu = c(M), Sigma = kronecker(V, U)), n, q)
  return(X)
}

## 0.2 Parameters

 n <- 100
 q <- 5
 p <- 1e5

 M <- matrix(0, n, q) # Means

 times <- 1E5

 hs2 <- 0.006
 hg2 <- 0.4
 alpha <- 0.7

 sigmaV = sqrt(1)    # Row and
 # sigmaU = sqrt(pi) # column std.deviations

## 1. Effect of the SNP

 X <- matrix(NA, n, p)
 mafs <- runif(p,  min = 0, max = 0.5)
 for (snp in 1:p){
   maf <- mafs[snp]
   X[ ,snp] <- sample(size = n, c(2,1,0), replace = T, prob = c(maf^2, 2*maf*(1-maf), (1-maf)^2 ) )
 }  

 # Select SNP
 sel <- 1 
 maf <- mafs[sel]
 Xi <- X[,sel]
 
 # Generate effects
 bi <- sample(c(-1,1), replace = T, size = q)

 # Multiply
 XB <- tcrossprod(scale(Xi), bi) * sqrt(hs2)

## 2. Relatedness signal 

 # Row covariance
 Rg <- tcrossprod(scale(X, center = T, scale = F))/p 
 
 # Problems if not positive definite, try via pivot Cholesky (rmatnorm) or Kronecker (rmatnorm2)
 # Rg <- tcrossprod(matrix(rnorm(n*n, mean = 0, sd = sigmaU),n,n)) 

 ## 2.2. Shared

  # Column covariance
  a <- matrix(rnorm(q*q, mean = 0, sd = sigmaV), q, q)
  aaT <- tcrossprod(a)
  aaT <- aaT  / (q*sigmaV^2) * (hg2*alpha)   # Rescale (q*sigmaV^2 is the theoretical couterpart to mean(diag(aaT)) )
  Vs <- aaT/mean(diag(Rg)) # Rescale so that Vsprime is equal to aaT (n*sigmaU^2 is the theoretical counterpart to mean(diag(Rg)), X needs to be centered!)

  # Generation
  # Gs <- rmatrixnorm(mean = M, U = Rg, V = Vs, n = times, list = T)         # Generate via Cholesky/MixMatrix (requires pos def)
  Gs <- lapply(1:times, function(x){rmatnorm(M = M, U = Rg, V = Vs)})        # Generate via Kronecker
  Vsprime <- Reduce('+', lapply(Gs, cov))/times; median(Vsprime/aaT)                    # Check: Should be 1
  mean(sapply(Gs, function(A){diag(cov(A))})); hg2*alpha                                # Check: Should be centered at hg2*alpha
  Usprime <- Reduce('+', lapply(Gs, function(A){cov(t(A))}))/times; median(Usprime/Rg)  # Check: Should be mean(diag(Vs))
  mean(diag(Vs))
  
 ## 2.3. Independent

  # Column covariance
  C <- diag(rnorm(q, mean = 0, sd = sigmaV)^2)
  C <- C / sigmaV^2 * (hg2*(1-alpha))   # Rescale. Note here we do not need q as there is not crossproduct
  Vi <- C/mean(diag(Rg)) # Rescale (2)

 # Generation
  # Gi <- rmatrixnorm(mean = M, U = Rg, V = Vi, n = times, list = T)         # Generate via Cholesky/MixMatrix (requires pos def)
  Gi <- lapply(1:times, function(x){rmatnorm2(M = M, U = Rg, V = Vi)})       # Generate via Kronecker
  Viprime <- Reduce('+', lapply(Gi, cov))/times; median(diag(Viprime/C))                # Check: Should be 1
  mean(sapply(Gi, function(A){diag(cov(A))})); hg2*(1-alpha)                            # Check: Should be centered at hg2*(1-alpha)
  Uiprime <- Reduce('+', lapply(Gi, function(A){cov(t(A))}))/times; median(Uiprime/Rg)  # Check: Should be mean(diag(Vi))
  mean(diag(Vi))
  
 ## 3. Residuals

 lambda = 0.6
 alphaH = alpha

 # Row covariance
 k <- 10
 L <- matrix(rnorm(n*k, mean = 0, sd = 1), n, k)
 LLT <- tcrossprod(L)

 ## 3.1. Structured noise (shared)

  # Col covariance
  aH <- matrix(rnorm(q*q, mean = 0, sd = sigmaV), q, q)
  aHaHT <- tcrossprod(aH)
  aHaHT <- aHaHT  / (q*sigmaV^2) * (lambda*alphaH*(1-hs2-hg2))   # Rescale 
  VHs <- aHaHT/mean(diag(LLT))

  Hs <- lapply(1:times, function(x){rmatnorm2(M = M, U = LLT, V = VHs)})     # Generate via Kronecker
  VHsprime <- Reduce('+', lapply(Hs, cov))/times; median(VHsprime/aHaHT)                  # Check: Should be 1
  mean(sapply(Hs, function(A){diag(cov(A))})); lambda*alphaH*(1-hs2-hg2)                  # Check: Should be centered at lambda*alphaH*(1-hs2-hg2)
  UHsprime <- Reduce('+', lapply(Hs, function(A){cov(t(A))}))/times; median(UHsprime/LLT) # Check: Should be mean(diag(VHs))
  mean(diag(VHs))

 ## 3.2. Structured noise (independent)

  # Col covariance
  CH <- diag(rnorm(q, mean = 0, sd = sigmaV)^2)
  CH <- CH / sigmaV^2 * (lambda*(1-alphaH)*(1-hs2-hg2)) # Rescale
  VHi <- CH/mean(diag(LLT)) # Rescale (2)

  Hi <- lapply(1:times, function(x){rmatnorm2(M = M, U = LLT, V = VHi)})     # Generate via Kronecker
  VHiprime <- Reduce('+', lapply(Hi, cov))/times; median(diag(VHiprime/CH))               # Check: Should be 1
  mean(sapply(Hi, function(A){diag(cov(A))})); lambda*(1-alphaH)*(1-hs2-hg2)              # Check: Should be lambda*(1-alphaH)*(1-hs2-hg2)
  UHiprime <- Reduce('+', lapply(Hi, function(A){cov(t(A))}))/times; median(UHiprime/LLT) # Check: Should be mean(diag(VHi))
  mean(diag(VHi))
  
 ## 3.3. Unstructured noise

 E <- lapply(1:times, function(x){ matrix(rnorm(n*q, mean = 0, sd = 1), n, q) * sqrt((1-lambda) * (1-hs2-hg2)) })
 VEprime <- Reduce('+', lapply(E, cov))/times; median(diag(VEprime))                      # Check: Should be (1-lambda) * (1-hs2-hg2)
 (1-lambda) * (1-hs2-hg2)

