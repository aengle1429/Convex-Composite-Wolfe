smooth.nesterov.h <- function(y) {
  n <- length(y)
  term1 <- 0.25 * y[1]^2
  term2 <- sum(y[-1]^2)
  h.y <- term1 + term2
  Dh.y <- matrix(NA, nrow=n, ncol=1)
  Dh.y[1, 1] <- 0.5*y[1]
  for(i in 2:n) {
    Dh.y[i, 1] <- 2 * y[i]
  }
  return(list(h.y, Dh.y))
}

nesterov.cheby.F <- function(x) {
  n <- length(x)
  DF.x <- diag(n)
  for(i in 2:n) {
    DF.x[i, i - 1] <- -4 * x[i - 1]
  }
  F.x <- c(x[1] - 1, x[-1] - 2 * x[-n]^2 + 1)  # replace smooth with 0.25 * (x_1 - 0.5)^2
  return(list(F.x, DF.x))
}

smooth.nesterov <- function(x) {
  F.data <- nesterov.cheby.F(x)
  F.x <- F.data[[1]]
  DF.x <- F.data[[2]]
  h.data <- smooth.nesterov.h(F.x)
  h.x <- h.data[[1]]
  Dh.x <- t(DF.x) %*% h.data[[2]]
  return(list(h.x, Dh.x))
}

smooth.nesterov.Delta <- function(x, d) {
  # h(F(x) + F'(x) d) - h(F(x))
  F.data <- nesterov.cheby.F(x)
  f.new <- smooth.nesterov.h(F.data[[1]] + F.data[[2]] %*% d)
  f.old <- smooth.nesterov(x)
  return(f.new[[1]] - f.old[[1]])
}

