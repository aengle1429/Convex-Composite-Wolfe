#bsource("C:/Users/aengl/Dropbox/A.Engle/cvx_comp_R/nesterov_chebyshev/nesterov_functions.R")
BFGS.update <- function(yk, sk, invBk) {
  n <- length(yk)
  pk <- 1 / as.numeric(t(yk) %*% sk)
  mat1 <- diag(n) - pk * outer(as.vector(sk), as.vector(yk))
  return(mat1 %*% invBk %*% t(mat1) + pk * outer(as.vector(sk), as.vector(sk)))
}

weak.wolfe <- function(x, d, c1, c2) {
  alpha <- 0
  beta <- Inf
  t <- 1
  
  f.data <- smooth.nesterov(x)
  
  while (TRUE) {
    xnew <- x + t * d
    fnew.data <- smooth.nesterov(xnew)
    
    if (fnew.data[[1]] > f.data[[1]] + c1 * t * t(f.data[[2]]) %*% d) {
      beta <- t
    } else if (c2 * t(f.data[[2]]) %*% d > t(fnew.data[[2]]) %*% d) {
      alpha <- t  
    } else {
      return(t)
    }
    
    if (beta == Inf) {
      t <- 2*t
    } else {
      t <- (alpha + beta) / 2
    }
  }
}

nm.weak.wolfe <- function(x, d, c1, c2, fmax) {
  alpha <- 0
  beta <- Inf
  t <- 1
  
  f.data <- smooth.nesterov(x)
  
  while (TRUE) {
    xnew <- x + t * d
    fnew.data <- smooth.nesterov(xnew, eps)
    
    if (fnew.data[[1]] > fmax + c1 * t * t(f.data[[2]]) %*% d) {
      beta <- t
    } else if (c2 * t(f.data[[2]]) %*% d > t(fnew.data[[2]]) %*% d) {
      alpha <- t  
    } else {
      return(t)
    }
    
    if (beta == Inf) {
      t <- 2*t
    } else {
      t <- (alpha + beta) / 2
    }
  }
}

cat("\f") # clear console
set.seed(12)

n <- 4  # n must be at least 2

# xk <- rnorm(n, 0, 2)
xhat <- rep(1, n)
xhat[1] <- -1
xk = xhat

p <- 5  # nonmonotone parameter

cat("Initializing at xk = ", xk, ".\n")
invBk <- diag(n)
total_iter <- 50000

funvals <- rep(0, total_iter + p)
f0 <- smooth.nesterov(xk)[[1]]
funvals[1:p] <- f0

for(i in 1:total_iter) {
  cat("iteration is ", i, "\n")
  
  fk.data <- smooth.nesterov(xk)
  dk <- -invBk %*% fk.data[[2]]

  # fmax = max(funvals[seq(i, i + p - 1)])
  
  step <- weak.wolfe(xk, dk, 0.1, 0.8)
  # step <- bt(zk, zk, 0.1, eta)
  # step <- nm.weak.wolfe(zk, dk, 0.1, 0.8, fmax)
  cat("step is ", step, "\n")
  sk <- step * dk

  xnew <- xk + sk
  fnew.data <- smooth.nesterov(xnew)
  
  yk <- fnew.data[[2]] - fk.data[[2]]
  invBk <- BFGS.update(yk, sk, invBk)
  
  cat("xnew is ", xnew, "\n")
  
  funvals[i + p] <- fnew.data[[1]]
  cat("Function value is ", fnew.data[[1]], "\n")
  xk <- xnew
  # if (nesterov.cheby(zk[1 : n])[[1]] < 1e-8)
  # {
  #   stop("Terminated on objective < ", 1e-8, ", after ", i, " iterations.\n")
  # }
}
