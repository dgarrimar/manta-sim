dM <- function (p, q) {  # Distance on the hyper-sphere , radius 1
  return(2*acos(sum(sqrt(p*q))))
}

dH <- function (p,q) {  # Hellinger distance
  a <- (sqrt(p) - sqrt(q))
  return(sqrt(sum(a*a)))
}

dE <- function (x1, x2) {
  a <- (x1-x2)
  b <- sqrt(sum(a*a))
  return(b)
}

interDist <- function(X) {
  n <- nrow(X)
  interdist <- matrix(data = 0, nrow = n, ncol = n)
  for(i in 2:n) {
    for(j in 1:(i-1)) {
      interdist[i, j] <- dE(X[i, ], X[j, ])
    }
  }
  return(interdist)
}

geodesic <- function (p, q, s) { # Geodesics on the hyper-sphere
  dm <- dM(p, q)
  r <- (sqrt(p)*cos(s/2)+(sqrt(q)-(sqrt(p)*cos(dm/2)))*sin(s/2)/sin(dm/2))^2  
  return(r)
}

step2distance <- function (p, q, step) { # Step, distance and H1 point
  if ( (p[1] + step <= 1) & (p[1] + step >= 0) & (step >= 0)) {
    dm <- dM(p,q)
    a    <- sqrt(p[1])   
    b    <- (sqrt(q[1])-a*cos(dm/2))/sin(dm/2)
    c    <- - sqrt(p[1] + step)
    termP <- (sqrt( (a^2)*(b^2)+(b^4)-(b^2)*(c^2))-a*c)/(a^2+b^2)
    termN <- - termP
    
    solP <- 2*acos(termP)
    solN <- 2*acos(termN)
    
    if (solP < 0) solP <- -solP
    if (solN < 0) solN <- -solN
    
    sol <- min(solP,solN)
    r <- geodesic(p,q,sol)
    dHell <- dH(p,r)
  } else{
    r <- rep(NA, 3)
    sol <- NA
    dHell <- NA
  }
  return(list("d" = sol, "r" = r , "dH" = dHell))
}

pv.f <- function(f, lambda, df.i, df.e, acc = 1e-14){
  
  pv.davies <- function(f, lambda, df.i, df.e, lim = 50000, acc = 1e-14){
    H <- c(rep(df.i, length(lambda)), rep(df.e, length(lambda)))
    pv <- suppressWarnings(davies(0, lambda = c(lambda, -f * lambda), h = H, lim = lim, acc = acc))
    if(pv$ifault != 0 || pv$Qq < 0 || pv$Qq > 1){
      return(pv)
    } else {
      return(pv$Qq)
    }
  }
  
  pv <- pv.davies(f = f, lambda = lambda, df.i = df.i, df.e = df.e, acc = acc)
  while (length(pv) > 1) {
    acc <- acc * 10
    pv  <- pv.davies(f = f, lambda = lambda, df.i = df.i, df.e = df.e, acc = acc)
  }
  if (pv < acc) {
    pv <- acc
  }
  return(c(pv, acc))
}

label <- function(a, b, n, u = 1, w = "B", plot = F, seed = 123) { # Get levels of the factors with unbalance u
  
  set.seed(seed)
  
  if (w == "A"){
    ua <- u
    ub <- 1
  } else if (w == "B") {
    ua <- 1
    ub <- u
  } else if (w == "AB") {
    ua <- ub <- u
  }
  
  xa <- c(ua, rep(1, a-1))
  pa <- round(xa/sum(xa)*n)
  r <- n - sum(pa)
  sel <- sample(1:length(pa), size = abs(r))
  if(r > 0){
    pa[sel] <- pa[sel] + 1
  } else if (r < 0) {
    pa[sel] <- pa[sel] - 1
  }
  A <- rep(1:a, times = pa)
  # A <- gl(a, n/a, length = n)
  
  B <- c()
  xb <- c(ub, rep(1, b-1))
  for (i in table(A)){
    pb <- round(xb/sum(xb)*i)
    r <- i - sum(pb)
    sel <- sample(1:length(pb), size = abs(r))
    if(r > 0){
      pb[sel] <- pb[sel] + 1
    } else if (r < 0) {
      pb[sel] <- pb[sel] - 1
    } 
    B <- c(B, rep(1:b, times = pb))
  }
  if(plot){
    pheatmap::pheatmap(cbind(A,B), cluster_cols = F, cluster_rows = F)
  }
  return(list(factor(A),factor(B)))
}

step2h1 <- function(p0, e, step) {
  # e should be one vertex of the simplex
  i <- which(e==1)
  if(p0[i] + step > 1 || p0[i] + step < 0) {
    stop("H1 out of the simplex.")
  } 
  p1 <- p0
  p1[i] <- p0[i] + step
  p1[-i] <- p0[-i] * (1 - step/(1-p0[i]))
  return(p1)
}

betasd <- function(a,b){
  sqrt((a*b) / ( (a+b)^2 * (a+b+1) ))
}

gammasd <- function(shape, scale){
  sqrt(shape*scale^2)
}

sim.simplex <- function(q, n, p0, stdev, tol = 1e-10, Y = NULL, dist = "norm"){
  
  if(is.null(Y)){
    Y <- t(matrix(p0, nrow = q, ncol = n))
  }
  
  elist <- list()
  for (i in 1:q){
    elist[[i]] <- rep(0, q)
    elist[[i]][i] <- 1 
  }
  
  for (i in 1:n){
    
    if(dist == "norm"){
      steps <- rnorm(q, mean = 0, sd = stdev) 
    } else if(dist == "gamma"){
      steps <- (rgamma(q, shape = 1, scale = 100) - 1)/gammasd(1, 100)*stdev
    } else if(dist == "beta"){
      steps <- (rbeta(q, shape1 = 0.5, shape2 = 0.5) - 0.5)/betasd(0.5, 0.5)*stdev
    } else {
      stop (sprintf("Unknown dist '%s'.", dist))
    }
    
    elist <- elist[sample(1:q)]
    for (j in 1:q) {
      k <- which(elist[[j]] == 1)
      if (steps[j] < -Y[i, k]){
        steps[j] <- -Y[i, k] + tol
        warning("One observation out of the simplex was corrected.")
      } else if(steps[j] > (1 - Y[i,k]) ){
        steps[j] <- (1 - Y[i, k]) - tol
        warning("One observation out of the simplex was corrected.")
      }
      Y[i, ] <- step2h1(Y[i, ], elist[[j]], steps[j])
    }
  }
  return(Y)
}

Sim.simplex <- function(B, q, n, loc, delta, hk, stdev, check = F, dist = "norm"){
  
  x <- c(loc, rep(1, q-1))
  y0 <- x/sum(x)
  
  if(check){
    e <- c(1, rep(0, q-1))
    y <- step2h1(y0, e, delta)
    Y <- sim.simplex(q, 1e4, y, stdev*hk, dist = dist)
    return(data.frame(exp = y, obs = colMeans(Y)))
  }
  
  Y <- sim.simplex(q, n, y0, stdev, dist = dist)
  
  if (delta != 0){
    
    # Generate e and H1 (+/- delta)
    e <- c(1, rep(0, q-1))
    y <- step2h1(y0, e, delta)
    yprime <- step2h1(y0, e, -delta)
    
    if(any(grepl(":", B))){ # Then we are changing the interaction
      levs <- levels(B)
      b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]   # Recover b
      B <- mapvalues(B, from = levs, to = 1:length(levs))             # Relevel
      Y[B == (1+b), ] <- sim.simplex(q, sum(B==(1+b)), yprime, stdev, dist = dist)
      Y[B == (2+b), ] <- sim.simplex(q, sum(B==(2+b)), y, stdev, dist = dist) 
    } 
    
    Y[B == 1, ] <- sim.simplex(q, sum(B==1), y, stdev*hk, dist = dist)   
    Y[B == 2, ] <- sim.simplex(q, sum(B==2), yprime, stdev, dist = dist)
    
  } else {
    
    if(any(grepl(":", B))){ # Then we are changing the interaction
      levs <- levels(B)
      b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]  # Recover b
      B <- mapvalues(B, from = levs, to = 1:length(levs))            # Relevel
    } 
    
    Y[B == 1, ] <- sim.simplex(q, sum(B==1), y0 , stdev*hk, dist = dist)
  }
  
  return(Y)
}

Sim.simplex.C <- function(B, C, q, n, loc, delta, hk, stdev, check = F, tol = 1e-10, dist = "norm"){
  
  if(check){
    n <- 1e4
    B <- rep(1, n)
  }
  
  x <- c(loc, rep(1, q-1))
  y0 <- x/sum(x)
  e <- c(1, rep(0, q-1))
  Y <- matrix(y0, nrow = n, ncol = q, byrow = T)
  
  if(delta != 0){
    Y[B==1, ] <- matrix(step2h1(y0, e, delta), nrow = sum(B==1), ncol = q, byrow = T)
    if(!check){
      Y[B==2, ] <- matrix(step2h1(y0, e, -delta), nrow = sum(B==2), ncol = q, byrow = T) 
    }
  }
  
  out <- 0
  for (k in 1:n){
    check2 <- Y[k, 1] + C[k]
    if (check2 < 0) {
      C[k] <- Y[k,1]
      out <- out + 1 
    } else if(check2 > 1) {
      C[k] <- 1 - Y[k,1]
      out <- out + 1 
    }
    Y[k,] <- step2h1(Y[k,], e, C[k])
  }
  if(out/n > 0.1) stop(">10% out")
  
  Y <- sim.simplex(q, n, p0, stdev, Y = Y, dist = dist)
  
  if(check){
    y <- step2h1(y0, e, delta)
    return(data.frame(exp = y, obs = colMeans(Y)))
  } else{
    return(Y)
  }
}

sim.mvnorm <- function(q, n, mu, v, c, tol = 1e-10){

  corr <- matrix(c, nrow = q, ncol = q)
  diag(corr) <- rep(1, q)
  sigma <- corr * tcrossprod(sqrt(v))
  
  if(any(eigen(sigma, only.values = T)$values < tol)){
    stop("Covariance matrix should be positive definite.")
  }
  
  Y <- mvrnorm(n, mu = mu, Sigma = sigma)
  
  return(Y)
}

sim.mvnorm.C <- function(q, n, mu, v, c, r, C_mean, C_var, tol = 1e-10){
  
  corr <- matrix(c, nrow = q, ncol = q)
  diag(corr) <- rep(1, q)
  
  mt <-c(1, r, rep(0, q-1))
  for (i in 1:ncol(corr)){
    if(i == 1){
      mt <- c(mt, r, corr[i,])
    } else {
      mt <- c(mt, 0, corr[i,])
    }
  }
  corr <- matrix(mt, ncol = q+1)
  sigma <- corr * tcrossprod(sqrt(c(C_var, v)))
  
  if(any(eigen(sigma, only.values = T)$values < tol)){
    stop("Covariance matrix should be positive definite.")
  }
  
  Yext <- mvrnorm(n, mu = c(C_mean, mu), Sigma = sigma)
  
  return(Yext) 
}

Sim.mvnorm <- function(B, q, n, mu, delta, hk, Var, Cor){
  
  if(Var == "equal"){
    vars <- rep(1, q)
  } else if (Var == "unequal"){
    vars <- (q:1)/sum(q:1)
  } else {
    stop(sprintf("Unknown option: Var = '%s'.", Var))
  }
  
  Y <- sim.mvnorm(q, n, mu, vars, Cor)
  
  if (delta != 0){
    
    if(any(grepl(":", B))){ # Then we are changing the interaction
      levs <- levels(B)
      b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]  # Recover b
      B <- mapvalues(B, from = levs, to = 1:length(levs))            # Relevel
      Y[B == (1+b), ] <- sim.mvnorm(q, sum(B==(1+b)), mu - delta, vars, Cor) 
      Y[B == (2+b), ] <- sim.mvnorm(q, sum(B==(2+b)), mu + delta, vars, Cor)
    } 
    
    Y[B == 1, ] <- sim.mvnorm(q, sum(B==1), mu + delta, vars*hk, Cor)   
    Y[B == 2, ] <- sim.mvnorm(q, sum(B==2), mu - delta, vars, Cor)      
    
  } else {
    
    if(any(grepl(":", B))){ # Then we are changing the interaction
      levs <- levels(B)
      b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]  # Recover b
      B <- mapvalues(B, from = levs, to = 1:length(levs))            # Relevel
    } 
    
    Y[B == 1, ] <- sim.mvnorm(q, sum(B==1), mu, vars*hk, Cor)
  }
  return(Y)
}

Sim.mvnorm.C <- function(B, q, n, mu, r, delta, hk, Var, Cor, C_mean, C_var){
  
  if(Var == "equal"){
    vars <- rep(1, q)
  } else if (Var == "unequal"){
    vars <- (q:1)/sum(q:1)
  } else {
    stop(sprintf("Unknown option: Var = '%s'.", Var))
  }
  
  Yext <- sim.mvnorm.C(q, n, mu, vars, Cor, r, C_mean, C_var)
  
  Y <- scale(Yext[,-1], center = T, scale = T)
  for (i in 1:q){Y[B != 1, i] <- Y[B != 1, i] * sqrt(vars[i]) + mu[i]}
  for (i in 1:q){Y[B == 1, i] <- Y[B == 1, i] * sqrt(vars[i]*hk) + mu[i]}
  
  if (delta != 0){
    Y[B == 1, ] <- Y[B == 1, ] + delta  
    Y[B == 2, ] <- Y[B == 2, ] - delta  
  } 
  
  Yext <- cbind(Yext[,1], Y)
  
  return(Yext)
}

Sim.multinom <- function(B, q, n, N, delta, loc) {
  
  x <- c(loc, rep(1, q-1))
  y0 <- x/sum(x)
  
  Y <- t(rmultinom(n, N, y0)) 
  
  if (delta != 0){
  
    # Generate e and H1 (+/- delta)
    e <- c(1, rep(0, q-1))
    y <- step2h1(y0, e, delta)
    yprime <- step2h1(y0, e, -delta)
    
    if(any(grepl(":", B))){ # Then we are changing the interaction
      levs <- levels(B)
      b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]   # Recover b
      B <- mapvalues(B, from = levs, to = 1:length(levs))             # Relevel
      Y[B == (1+b), ] <- t(rmultinom(sum(B==(1+b)), N, yprime)) 
      Y[B == (2+b), ] <- t(rmultinom(sum(B==(2+b)), N, y)) 
    } 
    
    Y[B == 1, ] <- t(rmultinom(sum(B==1), N, y))   
    Y[B == 2, ] <- t(rmultinom(sum(B==2), N, yprime))   
    
  } 
  return(Y)
}

sim.copula <- function(q, n, v, c, distdef, tol = 1e-10){
  
  corr <- matrix(c, nrow = q, ncol = q)
  diag(corr) <- rep(1, q)
  sigma <- corr * tcrossprod(sqrt(v))
  
  if(any(eigen(sigma, only.values = T)$values < tol)){
    stop("Covariance matrix should be positive definite.")
  }
  
  distrib <- unlist(strsplit(distdef, "-"))[1]
  params <- as.numeric(unlist(strsplit(distdef, "-"))[-1])
  
  mar <- switch(distrib[1],
                "unif" = list(min = params[1], max = params[2]),
                "gamma" = list(shape = params[1], scale = params[2]),
                "beta" = list(shape1 = params[1], shape2 = params[2]),
                "t" = list(df = params[1]),
                "exp" = list(rate = params[1]),
                "lnorm" = list(meanlog = params[1], sdlog = params[2]),
                "binom" = list(size = params[1], p = params[2]))

  myCop <- normalCopula(param = P2p(corr), dim = q, dispstr = "un")
  myMvd <- mvdc(copula = myCop, margins = rep(distrib[1], q), paramMargins = rep(list(mar), q))
  Y <- rMvdc(n, myMvd)
  return(Y)
}

sim.copula.C <- function(q, n, mu, v, c, distdef, r, C_mean, C_var, tol = 1e-10){
  
  corr <- matrix(c, nrow = q, ncol = q)
  diag(corr) <- rep(1, q)
  
  mt <-c(1, r, rep(0, q-1))
  for (i in 1:ncol(corr)){
    if(i == 1){
      mt <- c(mt, r, corr[i,])
    } else {
      mt <- c(mt, 0, corr[i,])
    }
  }
  corr <- matrix(mt, ncol = q+1)
  sigma <- corr * tcrossprod(sqrt(c(C_var, v)))
  
  if(any(eigen(sigma, only.values = T)$values < tol)){
    stop("Covariance matrix should be positive definite.")
  }
  
  distrib <- unlist(strsplit(distdef, "-"))[1]
  params <- as.numeric(unlist(strsplit(distdef, "-"))[-1])
  
  mar <- switch(distrib[1],
                "unif" = list(min = params[1], max = params[2]),
                "gamma" = list(shape = params[1], scale = params[2]),
                "beta" = list(shape1 = params[1], shape2 = params[2]),
                "t" = list(df = params[1]),
                "exp" = list(rate = params[1]),
                "lnorm" = list(meanlog = params[1], sdlog = params[2]),
                "binom" = list(size = params[1], p = params[2]))
  
  myCop <- normalCopula(param = P2p(corr), dim = q+1, dispstr = "un")
  myMvd <- mvdc(copula = myCop, margins = rep(distrib[1], q+1), paramMargins = rep(list(mar), q+1))
  Y <- rMvdc(n, myMvd)
  
  return(Y)
}

Sim.copula <- function(B, q, n, mu, delta, hk, Var, Cor, dd) {
  
  if(Var == "equal"){
    vars <- rep(1, q)
  } else if (Var == "unequal"){
    vars <- (q:1)/sum(q:1)
  } else {
    stop(sprintf("Unknown option: Var = '%s'.", Var))
  }
  
  Y <- sim.copula(q, n, vars, Cor, dd)

  Y <- scale(Y, center = T, scale = T)
  for (i in 1:q){Y[B != 1, i] <- Y[B != 1, i] * sqrt(vars[i]) + mu[i]}
  for (i in 1:q){Y[B == 1, i] <- Y[B == 1, i] * sqrt(vars[i]*hk) + mu[i]}
  
  if (delta != 0){
    if(any(grepl(":", B))){ # Then we are changing the interaction
      levs <- levels(B)
      b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]  # Recover b
      B <- mapvalues(B, from = levs, to = 1:length(levs))            # Relevel
      Y[B == (1+b), ] <-  Y[B == (1+b), ] - delta
      Y[B == (2+b), ] <- Y[B == (2+b), ] + delta
    } 
    
    Y[B == 1, ] <- Y[B == 1, ] + delta  
    Y[B == 2, ] <- Y[B == 2, ] - delta  
  } 
  
  return(Y)
}

Sim.copula.C <- function(B, q, n, mu, r, delta, hk, Var, Cor, dd, C_mean, C_var) {
  
  if(Var == "equal"){
    vars <- rep(1, q)
  } else if (Var == "unequal"){
    vars <- (q:1)/sum(q:1)
  } else {
    stop(sprintf("Unknown option: Var = '%s'.", Var))
  }
  
  Yext <- sim.copula.C(q, n, mu, vars, Cor, dd, r, C_mean, C_var)
  
  Y <- scale(Yext[,-1], center = T, scale = T)
  for (i in 1:q){Y[B != 1, i] <- Y[B != 1, i] * sqrt(vars[i]) + mu[i]}
  for (i in 1:q){Y[B == 1, i] <- Y[B == 1, i] * sqrt(vars[i]*hk) + mu[i]}
  
  if (delta != 0){
    Y[B == 1, ] <- Y[B == 1, ] + delta  
    Y[B == 2, ] <- Y[B == 2, ] - delta  
  } 
  
  Yext <- cbind(Yext[, 1], Y)
  
  return(Yext)
}

Sim.numcov <- function(x1, r = 0.5, m = 0, s = 1) {
  n     <- length(x1)
  theta <- acos(r)                              # corresponding angle
  x2    <- rnorm(n)                             # new random data
  X     <- cbind(x1, x2)                        # matrix
  Xctr  <- scale(X, center = TRUE, scale = F)   # center 
  
  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
  
  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
  # cor(x1, x)                                  # check                                   
  
  x <- as.numeric(scale(x, scale = T))          # scale x (it is already centered)
  x <- x * s + m                                # give custom mean and sd
  
  return(x)
}

