source("C:/Users/aengl/Dropbox/A.Engle/cvx_comp_R/nesterov_chebyshev/nesterov_functions.R")
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
  
  f.data <- nesterov.cheby(x)
  
  while (TRUE) {
    xnew <- x + t * d
    fnew.data <- nesterov.cheby(xnew)
    
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


cvx.comp.wolfe <- function(x, d, c1, c2) {
  alpha <- 0
  beta <- Inf
  t <- 1
  
  f.data <- nesterov.cheby(x)
  
  Delta.f <- nesterov.cheby.Delta(x, d)
  while (TRUE) {
    xnew <- x + t * d
    fnew.data <- nesterov.cheby(xnew)
    Delta.f.new <- nesterov.cheby.Delta(xnew, d)
    
    if (fnew.data[[1]] > f.data[[1]] + c1 * t * Delta.f) {
      beta <- t
    } else if (c2 * Delta.f > Delta.f.new) {
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
n <- 2
xk <- rnorm(n, 0, 100)
cat("Initializing at xk = ", xk)
invBk <- diag(n)
for(i in 1:100) {
  cat("iteration is ", i, "\n")
  
  fk.data <- nesterov.cheby(xk)
  dk <- -invBk %*% fk.data[[2]]
  # cat("Delta is ", nesterov.cheby.Delta(xk, dk), "\n")
  
  # step <- cvx.comp.wolfe(xk, dk, 0.2, 0.8, w)
  step <- weak.wolfe(xk, dk, 0.0, 0.5)
  # cat("step is ", step, "\n")
  
  xnew <- xk + step * dk
  cat("xnew is ", xnew, "\n")
  
  fnew.data <- nesterov.cheby(xnew)
  
  cat("Function value is ", fnew.data[[1]], "\n")
  
  yk <- fnew.data[[2]] - fk.data[[2]]
  sk <- step * dk
  invBk <- BFGS.update(yk, sk, invBk)
  xk <- xnew
}
