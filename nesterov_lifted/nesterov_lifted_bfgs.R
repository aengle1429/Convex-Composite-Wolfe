#bsource("C:/Users/aengl/Dropbox/A.Engle/cvx_comp_R/nesterov_chebyshev/nesterov_functions.R")
BFGS.update <- function(yk, sk, invBk) {
  n <- length(yk)
  pk <- 1 / as.numeric(t(yk) %*% sk)
  mat1 <- diag(n) - pk * outer(as.vector(sk), as.vector(yk))
  return(mat1 %*% invBk %*% t(mat1) + pk * outer(as.vector(sk), as.vector(sk)))
}

weak.wolfe <- function(x, d, c1, c2, eps) {
  alpha <- 0
  beta <- Inf
  t <- 1
  
  f.data <- nesterov.lifted.helper(x, eps)
  
  while (TRUE) {
    xnew <- x + t * d
    fnew.data <- nesterov.lifted.helper(xnew, eps)
    
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

nm.weak.wolfe <- function(x, d, c1, c2, eps, fmax) {
  alpha <- 0
  beta <- Inf
  t <- 1
  
  f.data <- nesterov.lifted.helper(x, eps)
  
  while (TRUE) {
    xnew <- x + t * d
    fnew.data <- nesterov.lifted.helper(xnew, eps)
    
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
# set.seed(12)

n <- 10  # n must be at least 2
problem.size <- 2 * n - 1

xk <- rnorm(n, 0, 2)
lamk <- 20*rep(1, n - 1)  # have to go to something that looks like epsilon
eps <- 0.01

p <- 5  # nonmonotone parameter

zk <- c(xk, lamk)

cat("Initializing at xk = ", xk, " and lamk = ", lamk)
invBk <- diag(problem.size)
total_iter <- 50000

funvals <- rep(0, total_iter + p)
f0 <- nesterov.lifted.helper(zk, eps)[[1]]
funvals[1:p] <- f0

for(i in 1:total_iter) {
  cat("iteration is ", i, "\n")
  
  fk.data <- nesterov.lifted.helper(zk, eps)
  dk <- -invBk %*% fk.data[[2]]

  # fmax = max(funvals[seq(i, i + p - 1)])
  
  step <- weak.wolfe(zk, dk, 0.1, 0.8, eps)
  # step <- bt(zk, zk, 0.1, eta)
  # step <- nm.weak.wolfe(zk, dk, 0.1, 0.8, fmax)
  cat("step is ", step, "\n")
  sk <- step * dk

  znew <- zk + sk
  fnew.data <- nesterov.lifted.helper(znew, eps)
  
  yk <- fnew.data[[2]] - fk.data[[2]]
  invBk <- BFGS.update(yk, sk, invBk)
  
  cat("znew is ", znew, "\n")
  
  funvals[i + p] <- fnew.data[[1]]
  cat("Function value is ", fnew.data[[1]], "\n")
  cat("Derivative is ", norm(fnew.data[[2]], "F"), "\n")
  zk <- znew
  # if (nesterov.cheby(zk[1 : n])[[1]] < 1e-8)
  # {
  #   stop("Terminated on objective < ", 1e-8, ", after ", i, " iterations.\n")
  # }
}
