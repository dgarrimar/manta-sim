rmatnorm_C <- function(M, U, V, tol = 1e-12){
    # Sample from Matrix Normal distribution via Cholesky
    # https://en.wikipedia.org/wiki/Matrix_normal_distribution#Drawing_values_from_the_distribution
    # Fast but requires U and V positive definite. 
    # Works well for positive semidefinite matrices via pivoting, when not many eigenvalues are ~ 0.
    # Otherwise, better use Kronecker product
  
    a <- nrow(M)
    b <- ncol(M)
    Z <- matrix(rnorm(a*b, 0, 1), a, b)
    L1 <- suppressWarnings(chol(U, pivot = T))
    L1 <- L1[, order(attr(L1, "pivot"))]
    L2 <- suppressWarnings(chol(V, pivot = T))
    L2 <- L2[, order(attr(L2, "pivot"))]
  
    return(M + crossprod(L1, Z) %*% L2)
}

getBeta <- function(q, b = "equal", ub = 2, a = q, fr = 0.75){ 
  
    # Generate effects
    if(b == "equal"){ # a <= q
        B <- rep(0, q)
        B[1:a] <- 1
        return(matrix(B, 1, q))
    } else if (b == "unequal"){
        return(matrix(seq(from = 1, to = ub, length.out = q), 1, q))
    } else if (b == "block"){
        w <- max(1, round(fr*q))
        return(matrix(c(rep(1, ceiling(w/2)), rep(-1, floor(w/2)), rep(0, q-w)), 1, q))
    } else if (b == "simplex" || b == "multinom"){
        return(matrix(c(1, -rep(1/(q-1), q-1)), 1, q))
    }
}

getCov <- function(q, v, c, u, B, tol = 1e-10){
  
    if (v == 'random') {
        sigma <- tcrossprod(matrix(rnorm(q^2), q, q))
        return(sigma)
    } else if(v == 'equal'){
        vars <- rep(1, q)
    } else if (v == 'unequal'){
        # vars <- (q:1)/sum(q:1)
        vars <- seq(from = 1, to = u, length.out = q)
    } else {
        stop(sprintf("Unknown option: Var = '%s'.", v))
    }
  
    # R <- matrix(c, nrow = q, ncol = q)
    R <- tcrossprod(B) * c
    diag(R) <- rep(1, q)
    sigma <- R * tcrossprod(sqrt(vars))
  
    if(any(eigen(sigma, only.values = T)$values < tol)){
        stop("Covariance matrix should be positive definite.")
    }
  
    return(sigma)
}

sim.copula <- function(n, sigma, distdef){
  
    # Obtain correlation matrix
    q <- nrow(sigma)
    D <- diag(q)*sqrt(diag(sigma))
    R <- solve(D)%*%sigma%*%solve(D)
    # https://math.stackexchange.com/questions/186959/correlation-matrix-from-covariance-matrix
  
    # Build and sample from copula
    distrib <- unlist(strsplit(distdef, "-"))[1]
    params <- as.numeric(unlist(strsplit(distdef, "-"))[-1])
    mar <- switch(distrib[1],
                  "unif" = list(min = params[1], max = params[2]),
                  "gamma" = list(shape = params[1], scale = params[2]),
                  "beta" = list(shape1 = params[1], shape2 = params[2]),
                  "t" = list(df = params[1]),
                  "exp" = list(rate = params[1]),
                  "lnorm" = list(meanlog = params[1], sdlog = params[2]),
                  "binom" = list(size = params[1], p = params[2]),
                  "norm" = list(mean = params[1], sd = params[2]))
    myCop <- normalCopula(param = P2p(R), dim = q, dispstr = "un")
    myMvd <- mvdc(copula = myCop, margins = rep(distrib[1], q),
                paramMargins = rep(list(mar), q))
    E <- rMvdc(n, myMvd) # This has sigma = R
    E <- t(apply(scale(E),1,function(e){e*sqrt(diag(sigma))})) 
  
    return(E)
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
