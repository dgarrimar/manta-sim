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

label <- function(a, b, n, u = 1, plot = F) { # Get levels of the factors with unbalance u
  
  # xa <- c(u, rep(1, a-1))
  # pa <- round(xa/sum(xa)*n)
  # pa[1] = pa[1] + n - sum(pa)
  # A <- rep(1:a, times = pa)
  A <- gl(a, n/a, length = n)
  
  B <- c()
  xb <- c(u, rep(1, b-1))
  for (i in table(A)){
    pb <- round(xb/sum(xb)*i)
    pb[1] = pb[1] + i - sum(pb)
    B <- c(B, rep(1:b, times = pb))
  }
  if(plot){
    pheatmap::pheatmap(cbind(A,B), cluster_cols = F, cluster_rows = F)
  }
  return(list(factor(A),factor(B)))
}

sim.props <- function(q, n, p0, stdev){
  p <- p0
  Y <- matrix(NA, nrow = n, ncol = q)
  for (i in 1:n){
    for (j in 1:q){
      e <- rep(0, q); e[j] <- 1
      d <- rnorm(1, mean = 0, sd = stdev)
      if(d > dM(p, e)){ d <- dM(p, e) }
      p <- geodesic(p, e, d)
    }
    Y[i, ] <- p
    p <- p0
  }
  return(Y)
}

Sim.props <- function(B, q, n, loc, delta, hk, stdev){
  
  x <- c(loc, rep(1, q-1))
  y0 <- x/sum(x)
  
  Y <- sim.props(q, n, y0, stdev)
  
  if (hk == 1){
    if (delta != 0){
      e <- c(1, rep(0, q-1))
      y <- step2distance(y0, e, delta)$r
      Y[B == 1, ] <- sim.props(q, sum(B==1), y, stdev)
    }
  } else {
    if (delta != 0){
      e <- c(1, rep(0, q-1))
      y <- step2distance(y0, e, delta)$r
      Y[B == 1, ] <- sim.props(q, sum(B==1), y, stdev*hk)
    } else {
      Y[B == 1, ] <- sim.props(q, sum(B==1), y0 , stdev*hk)
    }
  }
  return(Y)
}

sim.mvnorm <- function(q, n, mu, v, c){

  corr <- matrix(c, nrow = q, ncol = q)
  diag(corr) <- rep(1, q)
  sigma <- corr * tcrossprod(sqrt(v))
  
  Y <- mvrnorm(n, mu = mu, Sigma = sigma)
  
  return(Y)
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
  
  if (hk == 1){
    if (delta != 0){
      mu[1] <- mu[1] + delta
      Y[B == 1, ] <- sim.mvnorm(q, sum(B==1), mu, vars, Cor)
    }
  } else {
    if (delta != 0){
      mu[1] <- mu[1] + delta
      Y[B == 1, ] <- sim.mvnorm(q, sum(B==1), mu, vars*hk, Cor)
    } else {
      Y[B == 1, ] <- sim.mvnorm(q, sum(B==1), mu, vars*hk, Cor)
    }
  } 
  return(Y)
}
