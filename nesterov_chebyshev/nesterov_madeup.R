library(Matrix)
library(Rcplex)

nesterov.prox.subprob.bfgs <- function(eta, alpha, t, x, B) {
  n <- length(x)
  # this function finds d from the GN update via minimizing based at x with backtracking
  # solve strongly cvx QP   Fx.data <- nesterov.cheby.F(x)
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
  while (TRUE) {
    Q11 <- 0.5 * outer(as.vector(DFx.1), as.vector(DFx.1)) + 1 / (alpha * t) * diag(n) + B
    Q <- bdiag(Q11, 0 * diag(rep(1, n - 1))) # n = 2 gives errors for diag(0)
    
    ans <- Rcplex(cvec, A, b, Q, lb=-Inf, ub=Inf, control=list(trace=0))
    d <- matrix(ans[[1]][1:length(x)], nrow = length(x)) # extract only (d1, d2)
    if (nesterov.cheby(x + d)[[1]] <= norm(d, type="F")^2 / (2 * t) + nesterov.cheby.h(Fx.data[[1]] + Fx.data[[2]] %*% d)[[1]]) {
      break
    } else {
      t <- eta * t
      cat("BACKTRACKING!\n")
    }
  }
  return(list(1 / (alpha * t), t, d))
}

nesterov.prox.subprob.momentum <- function(x, v, t, alpha) { # no need for v when g=0
  n <- length(x)
  Fx.data <- nesterov.cheby.F(x)
  # linearize at x (y in the paper), update v!
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
  b <- rbind(-Fx.2, Fx.2) / alpha
  cvec <- c(0.5 * Fx.1 * DFx.1 / alpha, rep(1, n - 1))
  Q11 <- 0.5 * outer(as.vector(DFx.1), as.vector(DFx.1)) / alpha + diag(n) / t
  Q <- bdiag(Q11, 0 * diag(rep(1, n - 1))) # n = 2 gives errors for diag(0)
  ans <- Rcplex(cvec, A, b, Q, lb=-Inf, ub=Inf, control=list(trace=0))
  # cat("ans is ", ans[[1]], "\n")
  # cat("obj value is ", ans[[2]], "\n")
  # cat("status is ", ans[[3]], "\n")
  d <- matrix(ans[[1]][1:length(x)], nrow = length(x)) # extract only (d1, d2)
  return(d)
}

cat("\f")
n <- 4
# set.seed(8)
eta <- 0.5
alpha <- 0.5
xk <- rnorm(n, 0, 150)
Bk <- diag(n)

# xk <- c(0.4680186805, -0.5619170299, -0.3684985032, -0.7284177063, 0.0611847097)
# xk <- c(-0.9841175122, 0.9369745743, 0.7558427059, 0.142596392, -0.959332538, 0.8406378179, 0.4133438813, -0.6582936735, -0.1332988838, -0.96446285)
tk <- 1
for(k in 1:30000) {
  # cat("iteration is ", k, "\n")
  # get this from solving the subproblem based on y
  bt.data <- nesterov.prox.subprob.bfgs(eta, alpha, tk, xk, Bk)
  muk <- bt.data[[1]]
  tk <- bt.data[[2]]
  dk <- bt.data[[3]]
  step <- cvx.comp.wolfe(xk, d1k, 0.000001, 0.8)
  xk <- xk + step*dk # this is the S mapping
  cat("xk is ", xk, ".\nproxgrad norm is ", norm(dk, type="F")^2 * muk, "\n")
  cat("Delta is ", nesterov.cheby.Delta(yk, dk), "\n")
  fnew <- nesterov.cheby(xk)
  cat("function value is ", fnew[[1]], "\n\n\n")
  
  Fnew.data <- nesterov.cheby.F(xk)
  if(tk > 1) {
    hFnew.data <- nesterov.cheby.h(Fnew.data[[1]] + Fnew.data[[2]] %*% dk)
  }
  else {
    hFnew.data <- nesterov.cheby.h(Fnew.data[[1]] + Fnew.data[[2]] %*% dk * tk)
  }
  
  yk <- t(Fnew.data[[2]]) %*% hFnew.data[[2]] - t(Fk.data[[2]]) %*% hFk.data[[2]]
  Bk <- BFGS.update(yk, dk, step, Bk)
  
}

