library(Matrix)
library(Rcplex)
cat("\f") # clear console
n <- 10
xk <- rnorm(n, 1, 1)
Hk <- diag(n)
deltak <- 10
gam1 <- 0.5
gam2 <- 0.75
gam3 <- 1.1

beta1 <- 0.3
beta2 <- 0.6
beta3 <- 0.9

cat("Initializing at ", xk, "\n")

nesterov.jim.TRS <- function(x, H, delta) {
  n <- length(x)
  # this function finds d from the GN update via a trust region QCQP
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
  Q11 <- 0.5 * outer(as.vector(DFx.1), as.vector(DFx.1)) + H
  Q <- bdiag(Q11, 0 * diag(rep(1, n - 1))) # n = 2 gives errors for diag(0)
  QC = list(
    QC = list(Q = list(bdiag(diag(n), 0 * diag(rep(1, n - 1)))), L = NULL), dir = "L", b = delta
  )
  ans <- Rcplex_solve_QCP(cvec, A, b, Q, QC, lb=-Inf, ub=Inf, control=list(trace=0))
  d <- matrix(ans[[1]][1:n], nrow = length(x)) # extract only (d1, d2)
  return(d)
}

cat("Initializing at xk = ", xk)

num_iter <- 3000

for(i in 1:num_iter) {
  cat("iteration is ", i, "\n")
  Fk.data <- nesterov.cheby.F(xk)
  hFk.data <- nesterov.cheby.h(Fk.data[[1]])
  dk <- nesterov.jim.TRS(xk, Hk, deltak)
  Delta <- nesterov.cheby.Delta(xk, dk)
  cat("Delta is ", Delta, "\n")
  
  stopifnot(Delta < 0)
  
  tk <- cvx.comp.wolfe(xk, dk, 0.1, 0.8)
  # cat("tk is ", tk, "\n")
  dk <- tk * dk
  
  xnew <- xk + dk
  
  cat("xnew is ", xnew, "\n")
  
  Fnew.data <- nesterov.cheby.F(xnew)
  fnew.data <- nesterov.cheby(xnew)
  
  
  
  rk <- (fnew.data[[1]] - hFk.data[[1]]) / (Delta + 0.5 * as.numeric(crossprod(dk, as.vector(Hk %*% dk))))
  
  if (rk > beta3) {
    deltak <- 0.5 * (deltak + gam3 * deltak)
  } else if (beta2 <= rk && rk <= beta3) {
    deltak <- deltak
  } else {
    deltak <- 0.5 * deltak * (gam1 + gam2)
  }
  
  if (rk < beta1) {
    xnew <- xk  # do not update
    Hk <- Hk  # can not choose Hk
  } else {
    xnew <- xk + dk
    Hk <- Hk  # can choose Hk
  }
  
  cat("\n\n\n")
  
  # tk <- cvx.comp.wolfe(xk, dk, 0.1, 0.8)
  # tk <- 1
  # cat("tk is ", tk, "\n")
  
  # hFnew.data <- nesterov.cheby.h(Fnew.data[[1]] + Fnew.data[[2]] %*% dk)
  
  # yk <- t(Fnew.data[[2]]) %*% hFnew.data[[2]]  - t(Fk.data[[2]]) %*% hFk.data[[2]]
  
  Fnew.data <- nesterov.cheby.F(xnew)
  fnew.data <- nesterov.cheby(xnew)
  
  cat("Function value is ", fnew.data[[1]], "\n")
  
    
  xk <- xnew
  cat("xk is ", xk, "\n")
}
