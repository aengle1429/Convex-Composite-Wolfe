library(Matrix)
library(Rcplex)

nesterov.prox.subprob <- function(eta, alpha, t, x) {
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
    Q11 <- 0.5 * outer(as.vector(DFx.1), as.vector(DFx.1)) + 1 / (alpha * t) * diag(n)
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
n <- 10
# set.seed(8)
eta <- 0.5
alpha <- 0.5
xk <- rnorm(n, 1, 10)
vk <- rnorm(n, 1, 30)
# xk <- rep(1, n)
# xk <- c(1.000012352,1.000049409,1.00019764,1.000790636,1.003163796,1.012675201,1.051022127,1.209295021,1.924788896,6.409624589)
tk <- 1
N <- 5000000
funvals <- rep(0, N)
xks <- matrix(0, nrow = n, ncol = N)
proxgrads <- rep(0, N)
for(k in 1:N) {
  cat("iteration is ", k, "\n")
  ak <- 2 / (k + 1)
  yk <- ak * vk + (1 - ak) * xk
  # get this from solving the subproblem based on y
  bt.data <- nesterov.prox.subprob(eta, alpha, tk, yk)
  muk <- bt.data[[1]]
  tk <- bt.data[[2]]
  # cat("tk is ", tk, "\n")
  d1k <- bt.data[[3]]
  
  ## adding nonsense...
  xk <- yk + d1k # this is the S mapping
  # step = cvx.comp.wolfe(xk, d1k, 0.1, 0.8)
  # xk <- yk + step * d1k
  ##
  
  xks[, k] <- xk
  d2k <- nesterov.prox.subprob.momentum(yk, vk, 1 / (muk * ak), ak)
  vk <- vk + d2k # this is the bivariate S mapping
  cat("xk is ", xk, ".\nproxgrad norm is ", norm(d1k, type="F")^2 * muk, "\n")
  cat("Delta is ", nesterov.cheby.Delta(yk, d1k), "\n")
  fnew <- nesterov.cheby(xk)
  cat("function value is ", fnew[[1]], "\n\n\n")
  funvals[k] <- fnew[[1]]
}



