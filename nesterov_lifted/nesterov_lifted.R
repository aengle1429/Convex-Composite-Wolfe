# this function is the analytic approx on the interior of its domain
psi <- function(y, lambda, eps) {
  # Psi(y, lam) is the quadratic over linear approximation to min |x|
  # Args:
  #   y: primal variable
  #   lambda: dual variable
  #   eps: regularization parameter
  # Returns:
  #   psi(y,lambda) and its derivative on intr(dom(psi)), NA else
  if (length(y) != length(lambda) || length(y) > 1) {
    stop("Error: incorrect argument dimensions.")
  }
  if (lambda > 0) {
    psi.y <- 0.5 * ((eps + y^2) / lambda + lambda)
    Dpsi.y <- y / lambda
    Dpsi.lam <- 0.5 * (-(eps + y^2) / lambda^2 + 1)
  } else if (y == lambda && y == 0) {
    psi.y <- 0
    Dpsi.y <- NA
    Dpsi.lam <- NA
  } else {
    psi.y <- Inf
    Dpsi.y <- NA
    Dpsi.lam <- NA
  }
  return(list(psi.y, c(Dpsi.y, Dpsi.lam)))
}

psi.helper <- function(data) {
  return(psi(data[1], data[2], data[3]))
}

Fi <- function(x, lambda, i) {
  # Computes the ith F matrix which is to be precomposed with psi.
  #
  # Args:
  #   x: The x vector
  #   lambda: The vector of lambads
  #   i: the index to determine which Fi. Makes sense for 1 <= i <= length(lambda)
  #
  # Returns:
  #   F_i(x,lambda) and the Jacobian DF_i(x,lambda)
  n <- length(x)
  if ((n - 1) != length(lambda)) {
    stop("Error: x vector and lambda vector do not have the same length.")
  }
  if ( 1 > i || i > (n - 1)) {
    stop("Error: index i is out of range.")
  }
  Fi.xl <- c(x[i + 1] - 2*x[i]^2 + 1, lambda[i])
  row1 <- rep(0, 2 * n - 1)
  row1[i] <- -4*x[i]
  row1[i + 1] <- 1
  row2 <- rep(0, 2 * n - 1)
  row2[n + i] <- 1
  DFi.xl <- rbind(row1, row2)
  return(list(Fi.xl, DFi.xl))
}

nesterov.lifted <- function(x, lambda, eps) {
  n <- length(x)
  if (n != length(lambda) + 1) {
    stop("Error argument dimensions do not match.")
  }
  stopifnot(n == (length(lambda) + 1))
  f1 <- 0.25 * (x[1] - 1)^2
  # obj <- f1
  obj <- f1 + x[1]^2 * x[2]^2
  gradient <- matrix(rep(0, 2 * n - 1), ncol = 1)
  # gradient[1] <- 0.5 * (x[1] - 1)
  gradient[1] <- 0.5 * (x[1] - 1) + 2 * x[2]^2 * x[1]
  gradient[2] <- 2 * x[1]^2 * x[2]

  # compute the n - 1 compositions for obj and grad
  for (i in 1 : (n - 1)) {
    Fi.data <- Fi(x, lambda, i)
    psi.data <- psi.helper(c(Fi.data[[1]], eps))
    obj <- obj + psi.data[[1]]
    gradient <- gradient + t(Fi.data[[2]]) %*% psi.data[[2]]
  }
  
  return(list(obj, gradient))
}

nesterov.lifted.helper <- function(xk.lamk, eps) {
  # Unpacks the xk and lambdak into separate vectors.
  m <- length(xk.lamk)  # this equals 2 n - 1 for some n.
  n <- (1 + m) / 2  # this particular n.
  return(nesterov.lifted(xk.lamk[1:n], xk.lamk[(n + 1) : m], eps))
}

