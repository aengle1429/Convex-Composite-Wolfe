library(Matrix) # for sparse matrices and block diagonal ones.
library(Rcplex)

options(digits = 10)

nesterov.newton.subprob <- function(x, xprev, dprev, deltak) {
  # Solves the full Hessian subproblems where
  # lambda^(k+1) in sd h(F(x^k) + DF(x^k)d^k). 
  n <- length(x)
  
  # constructs the Lagrange multiplier
  
  Fxprev.data <- nesterov.cheby.F(xprev)
  Fxprev <- Fxprev.data[[1]]
  Fxprev.1 <- Fxprev[1]
  Fxprev.2 <- matrix(Fxprev[-1], nrow = n - 1)
  
  DFxprev <- Fxprev.data[[2]]
  DFxprev.1 <- DFxprev[1, ] # a vector
  DFxprev.2 <- DFxprev[-1, ] # a matrix unless n = 2.
  if (n == 2) {
    DFxprev.2 <- matrix(DFxprev.2, nrow = 1)
  }
  
  lam.data <- nesterov.cheby.h(Fxprev + DFxprev %*% dprev)
  lam <- lam.data[[2]]
  
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
  Hk <- diag(c(-4 * lam[-1], 0))
  
  eye <- diag(n - 1)
  Atop <- cbind(DFx.2, -eye)
  Abot <- cbind(-DFx.2, -eye)
  
  A <- rbind(Atop, Abot)
  b <- rbind(-Fx.2, Fx.2)
  cvec <- c(0.5 * Fx.1 * DFx.1, rep(1, n - 1))
  Q11 <- 0.5 * outer(as.vector(DFx.1), as.vector(DFx.1))#   + Hk
  lam.min <- min(eigen(Q11)$values)
  if (lam.min <= 1e-6) {
    Q11 <- Q11 + (abs(lam.min) + 1e-4) * diag(n)
  }
  Q22 <- 0 * diag(rep(1, n - 1))
  Q <- bdiag(Q11, Q22) # n = 2 gives errors for diag(0)
  QC = list(
    QC = list(Q = list(bdiag(diag(n), 0 * diag(rep(1, n - 1)))), L = NULL), dir = "L", b = deltak
  )
  # ans <- Rcplex_solve_QCP(cvec, A, b, Q, QC, lb=-Inf, ub=Inf, control=list(trace=0))
  ans <- Rcplex(cvec, A, b, Q, lb=-Inf, ub=Inf, control=list(trace=0))
  d <- matrix(ans[[1]][1:n], nrow = length(x)) # extract only (d1, d2)
  return(d)
}

cat("\f") # clear console
n <- 100

set.seed(1428)

deltak <- 1  # initial trust region radius

xk <- rnorm(n, 1, 100)
xprev <- xk

# xk <- c(1.000012352,1.000049409,1.00019764,1.000790636,1.003163796,1.012675201,1.051022127,1.209295021,1.924788896,6.409624589)

dk <- rep(0, n)
dprev <- dk
cat("Initializing at xk = ", xk)

xks <- matrix(0, nrow = n, ncol = N)

f0 <- nesterov.cheby(xk)[[1]]

num_iters = 20
for(i in 1:num_iters) {
  cat("iteration is ", i, "\n")
  
  dprev <- dk
  dk <- nesterov.newton.subprob(xk, xprev, dprev, deltak)
  cat("Delta is ", nesterov.cheby.Delta(xk, dk), "\n")
  xprev <- xk
  # step <- cvx.comp.wolfe(xk, dk, 0.0001, 0.9)
  step <- 1
  cat("step is ", step, "\n")
  xk <- xk + step * dk
  cat("xnew is ", xk, "\n")
  fnew.data <- nesterov.cheby(xk)
  cat("Function value is ", fnew.data[[1]], "\n")
  cat("\n\n\n")
}
