library(Matrix)
library(Rcplex)
BFGS.update <- function(yk, dk, tk, Bk) {
  yksk <- as.numeric(crossprod(yk, tk*dk))
  stopifnot(yksk > 0)
  n <- length(yk)
  vec1 <- as.vector(Bk %*% dk)
  
  # yhat <- yk / sqrt(yksk)
  # shat <- sk / sqrt(as.numeric(crossprod(sk, vec1)))
  
  # vec1hat <- Bk %*% shat
  # return(Bk - tcrossprod(vec1hat, vec1hat) +  tcrossprod(yhat, yhat))
  return(Bk - vec1 %o% vec1 / as.numeric(crossprod(dk, vec1)) +  as.vector(yk)%o%as.vector(yk) / yksk)
}

DFP.update <- function(yk, dk, tk, Bk) {
  yksk <- as.numeric(crossprod(yk, tk*dk))
  stopifnot(yksk > 0)
  n <- length(yk)
  
  mat1 <- diag(n) - tcrossprod(yk, dk) / as.numeric(crossprod(yk, dk))
  return(mat1 %*% Bk %*% t(mat1) + as.vector(yk)%o%as.vector(yk) / yksk)
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

cvx.comp.wolfe <- function(x, d, c1, c2) {
  alpha <- 0
  beta <- Inf
  t <- 1
  
  f.data <- nesterov.cheby(x)
  
  Delta.f <- nesterov.cheby.Delta(x, d)
  while (TRUE) {
    # cat("Trying t = ", t, ".\n")
    xnew <- x + t * d
    fnew.data <- nesterov.cheby(xnew)
    Delta.f.new <- nesterov.cheby.Delta(xnew, d)
    
    if (fnew.data[[1]] > f.data[[1]] + c1 * t * Delta.f) {
      beta <- t
    } else if (c2 * Delta.f > Delta.f.new) {
      alpha <- t  
      # cat("Have sufficient decrease.\n")
    } else {
      # cat("Final step is ", t, ".\n")
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
# set.seed(10)
# set.seed(8)
n <- 4
xk <- rnorm(n, 0, 150)
# xk <- rnorm(n, 1, 1)
Bk <- diag(n)
eta <- 0.5
alpha <- 0.5
L <- 5 # limited memory thing
# seq(i - L + 1, i)
cat("Initializing at xk = ", xk)
memory_yk <- matrix(0, nrow = n, ncol = L)
memory_sk <- matrix(0, nrow = n, ncol = L)
for(i in 1:20000) {
  cat("iteration is ", i, "\n")
  Fk.data <- nesterov.cheby.F(xk)
  hFk.data <- nesterov.cheby.h(Fk.data[[1]])
  dk <- nesterov.bfgs.subprob(xk, Bk)
  
  cat("\n\n\n")
  Delta <- nesterov.cheby.Delta(xk, dk)
  cat("Delta is ", Delta, "\n")
  stopifnot(Delta < 0)
  
  tk <- cvx.comp.wolfe(xk, dk, 0.00000001, 0.00008)
  cat("tk is ", tk, "\n")
  
  xnew <- xk + tk * dk
  # cat("xnew is ", xnew, "\n")
  
  Fnew.data <- nesterov.cheby.F(xnew)
  if(tk > 1) {
    hFnew.data <- nesterov.cheby.h(Fnew.data[[1]] + Fnew.data[[2]] %*% dk)
  }
  else {
    hFnew.data <- nesterov.cheby.h(Fnew.data[[1]] + Fnew.data[[2]] %*% dk * tk)
  }
  yk <- t(Fnew.data[[2]]) %*% hFnew.data[[2]] - t(Fk.data[[2]]) %*% hFk.data[[2]]
  
  idx = (i %% L) + 1
  memory_sk[, idx] <- tk*dk
  memory_yk[, idx] <- yk
  
  Bk <- BFGS.update(yk, dk, tk, Bk)
  cat("Lambda min is ", eigen(Bk, only.values = TRUE)$values)
  # Bk <- DFP.update(yk, dk, 1, Bk)
  cat("Condition number of Bk is ", kappa(Bk), ".\n")
  
  # cat("Eigenvalues of Bk are ", eigen(Bk)$values, ".\n")
  fnew.data <- nesterov.cheby(xnew)
  
  cat("Function value is ", fnew.data[[1]], "\n")
  
  xk <- xnew
}
cat(xk)

