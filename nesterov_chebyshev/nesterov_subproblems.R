library(Matrix) # for sparse matrices and block diagonal ones.
library(Rcplex)

options(digits = 10)
nesterov.subprob <- function(x) {
  n <- length(x)
  # this function finds d from the GN update via
  # solving a quadratic program. notation is 0.5 x^T Q x + c^Tx s.t. Ax <= b
  Fx.data <- nesterov.cheby.F(x)
  Fx <- Fx.data[[1]]
  Fx.1 <- Fx[1]
  Fx.2 <- matrix(Fx[-1], nrow = n - 1)
  
  DFx <- Fx.data[[2]]
  DFx.1 <- DFx[1, ] # a vector
  DFx.2 <- DFx[-1, ] # a matrix unless n = 2.
  if (n == 2) {
    DFx.2 <- matrix(DFx.2, nrow = 1)
  }
  
  eye <- diag(n - 1)
  Atop <- cbind(DFx.2, -eye)
  Abot <- cbind(-DFx.2, -eye)
  
  A <- rbind(Atop, Abot)
  b <- rbind(-Fx.2, Fx.2)
  cvec <- c(0.5 * Fx.1 * DFx.1, rep(1, n - 1))
  Q <- bdiag(0.5 * outer(as.vector(DFx.1), as.vector(DFx.1)), 0 * diag(rep(1, n - 1))) # n = 2 gives errors for diag(0)
  ans <- Rcplex(cvec, A, b, Q, lb=-Inf, ub=Inf, control=list(trace=0))
  # cat("ans is ", ans[[1]], "\n")
  # cat("obj value is ", ans[[2]], "\n")
  # cat("status is ", ans[[3]], "\n")
  d <- matrix(ans[[1]][1:n], nrow = length(x)) # extract only (d1, d2)
  return(d)
}

cvx.comp.wolfe <- function(x, d, c1, c2) {
  alpha <- 0
  beta <- Inf
  t <- 1
  
  f.data <- nesterov.cheby(x)
  
  Delta.f <- nesterov.cheby.Delta(x, d)
  while (TRUE) {
    cat("Trying t = ", t, ".\n")
    xnew <- x + t * d
    fnew.data <- nesterov.cheby(xnew)
    Delta.f.new <- nesterov.cheby.Delta(xnew, d)
    
    # if (fnew.data[[1]] > f.data[[1]] - exp(log(c1) + log(t) + log(-Delta.f))) {
    if (fnew.data[[1]] > f.data[[1]] + c1 * t * Delta.f) {
      beta <- t
    } else if (c2 * Delta.f > Delta.f.new) {# || c2 * Delta.f - Delta.f.new < 1e-17) {
    # } else if (-exp(log(c2) + log(-Delta.f)) > Delta.f.new) {  
      alpha <- t
      cat("WWII fails with c2Delf-Delfnew=", c2 * Delta.f - Delta.f.new, "\n")
    } else {
      return(t)
    }
    
    if (beta == Inf) {
      t <- 2*t
    } else {
      t <- (alpha + beta) / 2
      # t <- exp(log(alpha + beta) - log(2))
    }
  }
}

cvx.comp.bt <- function(x, d, c1) {
  t <- 1
  
  f.data <- nesterov.cheby(x)
  
  Delta.f <- nesterov.cheby.Delta(x, d)
  while (TRUE) {
    # cat("Trying t = ", t, ".\n")
    xnew <- x + t * d
    fnew.data <- nesterov.cheby(xnew)
    Delta.f.new <- nesterov.cheby.Delta(xnew, d)
    
    if (fnew.data[[1]] > f.data[[1]] + c1*t*Delta.f) {
      t <- 0.7 * t
    } else {
      cat("Final step size is ", t, ".\n")
      return(t)
    }
  }
}

cat("\f") # clear console
# set.seed(10)
set.seed(8)
n <- 10 # num vars
N <- 10000 # num iters
xk <- rnorm(n, 1, 0.1)
cat("Initializing at xk = ", xk)
p <- 5 # number of funvals to consider
funvals <- rep(0, N + p)
xks <- matrix(0, nrow = n, ncol = N)

f0 <- nesterov.cheby(xk)[[1]]
funvals[1:p] <- f0
for(i in 1:N) {
  cat("iteration is ", i, "\n")

  dk <- nesterov.subprob(xk)
  cat("\n\n\n")
  
  cat("Delta is ", nesterov.cheby.Delta(xk, dk), "\n")
  
  fmax = max(funvals[seq(i, i + p - 1)])
  
  # step <- nm.cvx.comp.wolfe(xk, dk, 0.000001, 0.00008, fmax)
  
  step <- cvx.comp.wolfe(xk, dk, 0.1, 0.8)
  
  # step <- cvx.comp.bt(xk, dk, 0.5)
  cat("step is ", step, "\n")
  
  xnew <- xk + step * dk
  cat("xnew is ", xnew, "\n")
  
  fnew.data <- nesterov.cheby(xnew)
  
  cat("Function value is ", fnew.data[[1]], "\n")
  
  xk <- xnew
}
