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

