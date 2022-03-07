step2h1 <- function(p0, e, step) {
    # e should be one vertex of the simplex
    i <- which(e == 1)
    if (p0[i] + step > 1 || p0[i] + step < 0) {
        stop("H1 out of the simplex.")
    } 
    p1 <- p0
    p1[i] <- p0[i] + step
    p1[-i] <- p0[-i] * (1 - step/(1 - p0[i]))
    return(p1)
}

betasd <- function(a,b) {
    sqrt((a*b) / ( (a+b)^2 * (a+b+1) ))
}

gammasd <- function(shape, scale) {
    sqrt(shape*scale^2)
}

sim.simplex <- function(q, n, p0, stdev, tol = 1e-10, Y = NULL, dist = "norm") {
  
    if (is.null(Y)) {
        Y <- t(matrix(p0, nrow = q, ncol = n))
    }
  
    elist <- list()
    for (i in 1:q) {
        elist[[i]] <- rep(0, q)
        elist[[i]][i] <- 1 
    }
  
    for (i in 1:n) {
    
        if (dist == "norm") {
            steps <- rnorm(q, mean = 0, sd = stdev) 
        } else if (dist == "gamma") {
            steps <- (rgamma(q, shape = 1, scale = 100) - 1)/gammasd(1, 100)*stdev
        } else if (dist == "beta") {
            steps <- (rbeta(q, shape1 = 0.5, shape2 = 0.5) - 0.5)/betasd(0.5, 0.5)*stdev
        } else {
            stop (sprintf("Unknown dist '%s'.", dist))
        }
    
        elist <- elist[sample(1:q)]
        for (j in 1:q) {
            k <- which(elist[[j]] == 1)
            if (steps[j] < -Y[i, k]) {
                steps[j] <- -Y[i, k] + tol
                warning("One observation out of the simplex was corrected.")
            } else if (steps[j] > (1 - Y[i,k])) {
                steps[j] <- (1 - Y[i, k]) - tol
                warning("One observation out of the simplex was corrected.")
            }
            Y[i, ] <- step2h1(Y[i, ], elist[[j]], steps[j])
        }
    }
    return(Y)
}

Sim.simplex <- function(B, q, n, loc, delta, hk, stdev, check = F, dist = "norm") {
  
    x <- c(loc, rep(1, q-1))
    y0 <- x/sum(x)
  
    if (check) {
        e <- c(1, rep(0, q-1))
        y <- step2h1(y0, e, delta)
        Y <- sim.simplex(q, 1e4, y, stdev*hk, dist = dist)
        return(data.frame(exp = y, obs = colMeans(Y)))
    }
  
    Y <- sim.simplex(q, n, y0, stdev, dist = dist)
  
    if (delta != 0) {
    
        # Generate e and H1 (+/- delta)
        e <- c(1, rep(0, q-1))
        y <- step2h1(y0, e, delta)
        yprime <- step2h1(y0, e, -delta)
    
        if (any(grepl(":", B))) { # Then we are changing the interaction
            levs <- levels(B)
            b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]   # Recover b
            B <- mapvalues(B, from = levs, to = 1:length(levs))             # Relevel
            Y[B == (1+b), ] <- sim.simplex(q, sum(B == (1+b)), yprime, stdev, dist = dist)
            Y[B == (2+b), ] <- sim.simplex(q, sum(B == (2+b)), y, stdev, dist = dist) 
        } 
    
        Y[B == 1, ] <- sim.simplex(q, sum(B == 1), y, stdev*hk, dist = dist)   
        Y[B == 2, ] <- sim.simplex(q, sum(B == 2), yprime, stdev, dist = dist)
    
    } else {
    
        if (any(grepl(":", B))) { # Then we are changing the interaction
            levs <- levels(B)
            b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]  # Recover b
            B <- mapvalues(B, from = levs, to = 1:length(levs))            # Relevel
        } 

        Y[B == 1, ] <- sim.simplex(q, sum(B == 1), y0 , stdev*hk, dist = dist)
    }

    return(Y)
}

sim.mvnorm <- function(q, n, mu, v, c, tol = 1e-10) {

    corr <- matrix(c, nrow = q, ncol = q)
    diag(corr) <- rep(1, q)
    sigma <- corr * tcrossprod(sqrt(v))
  
    if (any(eigen(sigma, only.values = T)$values < tol)) {
        stop("Covariance matrix should be positive definite.")
    }
  
    Y <- mvrnorm(n, mu = mu, Sigma = sigma)
  
    return(Y)
}

Sim.mvnorm <- function(B, q, n, mu, delta, hk, hk2, Var, Cor) {
  
    if (Var == "equal") {
        vars <- rep(1, q)
    } else if (Var == "unequal") {
        vars <- (q:1)/sum(q:1)
    } else {
        stop(sprintf("Unknown option: Var = '%s'.", Var))
    }
  
    Y <- sim.mvnorm(q, n, mu, vars, Cor)
  
    if (delta != 0) {
    
        if (any(grepl(":", B))) { # Then we are changing the interaction
            levs <- levels(B)
            b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]  # Recover b
            B <- mapvalues(B, from = levs, to = 1:length(levs))            # Relevel
            Y[B == (1+b), ] <- sim.mvnorm(q, sum(B == (1+b)), mu - delta, vars, Cor) 
            Y[B == (2+b), ] <- sim.mvnorm(q, sum(B == (2+b)), mu + delta, vars, Cor)
        } 
    
        Y[B == 1, ] <- sim.mvnorm(q, sum(B == 1), mu + delta, vars*hk, Cor*hk2)   
        Y[B == 2, ] <- sim.mvnorm(q, sum(B == 2), mu - delta, vars, Cor)      
    
    } else {
    
        if (any(grepl(":", B))) { # Then we are changing the interaction
            levs <- levels(B)
            b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]  # Recover b
            B <- mapvalues(B, from = levs, to = 1:length(levs))            # Relevel
        } 
    
        Y[B == 1, ] <- sim.mvnorm(q, sum(B == 1), mu, vars*hk, Cor*hk2)
    }
    return(Y)
}

Sim.multinom <- function(B, q, n, N, delta, loc) {
  
    x <- c(loc, rep(1, q-1))
    y0 <- x/sum(x)
  
    Y <- t(rmultinom(n, N, y0)) 
  
    if (delta != 0) {
  
        # Generate e and H1 (+/- delta)
        e <- c(1, rep(0, q-1))
        y <- step2h1(y0, e, delta)
        yprime <- step2h1(y0, e, -delta)
    
    if (any(grepl(":", B))) { # Then we are changing the interaction
        levs <- levels(B)
        b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]   # Recover b
        B <- mapvalues(B, from = levs, to = 1:length(levs))             # Relevel
        Y[B == (1+b), ] <- t(rmultinom(sum(B == (1+b)), N, yprime)) 
        Y[B == (2+b), ] <- t(rmultinom(sum(B == (2+b)), N, y)) 
    } 
    
    Y[B == 1, ] <- t(rmultinom(sum(B == 1), N, y))   
    Y[B == 2, ] <- t(rmultinom(sum(B == 2), N, yprime))   
    
    } 
    return(Y)
}

sim.copula <- function(q, n, v, c, distdef, tol = 1e-10) {
  
    corr <- matrix(c, nrow = q, ncol = q)
    diag(corr) <- rep(1, q)
    sigma <- corr * tcrossprod(sqrt(v))
  
    if (any(eigen(sigma, only.values = T)$values < tol)) {
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
                  "binom" = list(size = params[1], p = params[2])
    )

    myCop <- normalCopula(param = P2p(corr), dim = q, dispstr = "un")
    myMvd <- mvdc(copula = myCop, margins = rep(distrib[1], q), paramMargins = rep(list(mar), q))
    Y <- rMvdc(n, myMvd)
    return(Y)
}

Sim.copula <- function(B, q, n, mu, delta, hk, Var, Cor, dd) {
  
    if (Var == "equal") {
        vars <- rep(1, q)
    } else if (Var == "unequal") {
        vars <- (q:1)/sum(q:1)
    } else {
        stop(sprintf("Unknown option: Var = '%s'.", Var))
    }
  
    Y <- sim.copula(q, n, vars, Cor, dd)

    Y <- scale(Y, center = T, scale = T)
    for (i in 1:q) {Y[B != 1, i] <- Y[B != 1, i] * sqrt(vars[i]) + mu[i]}
    for (i in 1:q) {Y[B == 1, i] <- Y[B == 1, i] * sqrt(vars[i]*hk) + mu[i]}
  
    if (delta != 0) {
        if (any(grepl(":", B))) { # Then we are changing the interaction
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
