nesterov.cheby.h <- function(y) {
  # 0.25 * y_1^2 + |y_2| + ... + |y_n|
  # along with its derivative, where it makes sense
  n <- length(y)
  term1 <- 0.25 * y[1]^2
  term2 <- sum(abs(y[-1]))
  h.y <- term1 + term2
  Dh.y <- matrix(NA, nrow=n, ncol=1)
  Dh.y[1, 1] <- 0.5*y[1]
  for(i in 2:n) {
    if (y[i] > 0) {
      Dh.y[i, 1] <- 1 
    } else if (y[i] < 0) {
      Dh.y[i, 1] <- -1
    } else {
      Dh.y[i, 1] <- runif(1, min = -1, max = 1)  #TODO fix to return a specific subgrad
    }
  }
  return(list(h.y, Dh.y))
}

nesterov.cheby.F <- function(x) {
  # inner function F for the nesterov-chebyshev-rosenbrock
  n <- length(x)
  DF.x <- diag(n)
  for(i in 2:n) {
    DF.x[i, i - 1] <- -4 * x[i - 1]
  }
  F.x <- c(x[1] - 1, x[-1] - 2 * x[-n]^2 + 1)
  return(list(F.x, DF.x))
}

nesterov.cheby <- function(x) {
  # Convex composite h o F of the
  # nonsmooth Nesterov function from Overton-Lewis paper
  n <- length(x)
  F.data <- nesterov.cheby.F(x)
  F.x <- F.data[[1]]
  DF.x <- F.data[[2]]
  h.data <- nesterov.cheby.h(F.x)
  h.x <- h.data[[1]]
  Dh.x <- t(DF.x) %*% h.data[[2]]
  Hh.x <- diag(c(-4 * h.data[[2]][-1], 0))
  return(list(h.x, Dh.x, Hh.x))
}

nesterov.cheby.Delta <- function(x, d) {
  # continuous approx to directional derivative in the cvx-comp case:
  # h(F(x) + F'(x) d) - h(F(x))
  F.data <- nesterov.cheby.F(x)
  f.new <- nesterov.cheby.h(F.data[[1]] + F.data[[2]] %*% d)
  f.old <- nesterov.cheby(x)
  return(f.new[[1]] - f.old[[1]])
}

phi.n <- function(x, lam) {
  # experimental function using cvx-comp Lagrangian
  # TODO change name from phi
  F.data <- nesterov.cheby.F(x)
  return(crossprod(F.data[[1]], lam))
}

nesterov.cheby.Delta.new <- function(x, d, lam) {
  # experimental Delta function using cvx-comp Lagrangian
  F.data <- nesterov.cheby.F(x)
  return(crossprod(lam, F.data[[2]] %*% d))
}

