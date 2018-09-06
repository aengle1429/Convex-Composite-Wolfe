library(Matrix)
library(Rcplex)
nesterov.SR1.subprob <- function(x, B) {
  # Solves the full Hessian subproblems where
  # lambda^(k+1) in sd h(F(x^k) + DF(x^k)d^k). 
  n <- length(x)
  
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
  Q11 <- 0.5 * outer(as.vector(DFx.1), as.vector(DFx.1)) + B
  lam.min <- min(eigen(Q11)$values)
  if (lam.min <= 1e-6) {
    Q11 <- Q11 + (abs(lam.min) + 1e-4) * diag(n)
  }
  Q22 <- 0 * diag(rep(1, n - 1))
  Q <- bdiag(Q11, Q22) # n = 2 gives errors for diag(0)
  # QC = list(
  #   QC = list(Q = list(bdiag(diag(n), 0 * diag(rep(1, n - 1)))), L = NULL), dir = "L", b = deltak
  # )
  # ans <- Rcplex_solve_QCP(cvec, A, b, Q, QC, lb=-Inf, ub=Inf, control=list(trace=0))
  ans <- Rcplex(cvec, A, b, Q, lb=-Inf, ub=Inf, control=list(trace=0))
  d <- matrix(ans[[1]][1:n], nrow = length(x)) # extract only (d1, d2)
  return(d)
}


nesterov.bfgs.subprob <- function(x, B) {
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
  Q11 <- 0.5 * outer(as.vector(DFx.1), as.vector(DFx.1)) + B 
  Q <- bdiag(Q11, 0 * diag(rep(1, n - 1))) # n = 2 gives errors for diag(0)
  ans <- Rcplex(cvec, A, b, Q, lb=-Inf, ub=Inf, control=list(trace=0))
  # cat("ans is ", ans[[1]], "\n")
  # cat("obj value is ", ans[[2]], "\n")
  # cat("status is ", ans[[3]], "\n")
  d <- matrix(ans[[1]][1:n], nrow = length(x)) # extract only (d1, d2)
  return(d)
}

cat("\f") # clear console
n <- 10
xk <- rnorm(n, 0, 20)
# xk <- rnorm(n, 1, 1)
# xk<-c(1.000012352,1.000049409,1.00019764,1.000790636,1.003163796,1.012675201,1.051022127,1.209295021,1.924788896,6.409624589)
Bk <- diag(n) / 1000000
eta <- 0.5
alpha <- 0.5
cat("Initializing at xk = ", xk)
N <- 20000
p <- 5 # number of funvals to consider
funvals <- rep(0, N + p)
xks <- matrix(0, nrow = n, ncol = N)

f0 <- nesterov.cheby(xk)[[1]]
funvals[1:p] <- f0
for(i in 1:N) {
  cat("iteration is ", i, "\n")
  Fk.data <- nesterov.cheby.F(xk)
  hFk.data <- nesterov.cheby.h(Fk.data[[1]])
  # dk <- nesterov.bfgs.subprob(xk, Bk)
  dk <- nesterov.SR1.subprob(xk, Bk)
  cat("\n\n\n")
  Delta <- nesterov.cheby.Delta(xk, dk)
  cat("Delta is ", Delta, "\n")
  stopifnot(Delta < 0)
  
  tk <- cvx.comp.wolfe(xk, dk, 0.0001, 0.9)
  # tk <- 1
  cat("tk is ", tk, "\n")
  
  xnew <- xk + tk * dk

  hLin.data <- nesterov.cheby.h(Fk.data[[1]] + Fk.data[[2]] %*% (tk * dk))
  
  Fnew.data <- nesterov.cheby.F(xnew)
  if (tk > 1) {
    hFnew.data <- nesterov.cheby.h(Fnew.data[[1]] + Fnew.data[[2]] %*% dk)
  }
  else {
    hFnew.data <- nesterov.cheby.h(Fnew.data[[1]] + Fnew.data[[2]] %*% dk)
  }
  
  # yk <- t(Fnew.data[[2]]) %*% hFnew.data[[2]]  - t(Fk.data[[2]]) %*% hFk.data[[2]]  
  yk <- (t(Fnew.data[[2]]) - t(Fk.data[[2]])) %*% hLin.data[[2]]
  
  # Bk <- DFP.update(yk, dk, tk, Bk)
  Bk <- damped.BFGS.update(yk, dk, tk, Bk)
  # Bk <- SR1.update(yk, dk, tk, Bk)
  cat("Condition number of Bk is ", kappa(Bk), ".\n")
  
  # cat("Eigenvalues of Bk are ", eigen(Bk)$values, ".\n")
  fnew.data <- nesterov.cheby(xnew)
  
  cat("Function value is ", fnew.data[[1]], "\n")

  xk <- xnew
  cat("xk is ", xk, "\n")
  xks[, i] <- xk
}
