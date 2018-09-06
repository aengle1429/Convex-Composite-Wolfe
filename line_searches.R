#
# Implements classical and convex-composite backtracking and weak Wolfe
# line searches applied to the Nesterov Chebyshev-Rosenbrock Functions
#
cvx.comp.wolfe <- function(x, d, c1, c2) {
  # Line search developed in Burke-Engle 2018
  alpha <- 0
  beta <- Inf
  t <- 1
  
  f.data <- nesterov.cheby(x)
  
  Delta.f <- nesterov.cheby.Delta(x, d)
  while (TRUE) {
    # cat("Trying t = ", t, ".\n")
    xnew <- x + t * d
    fnew.data <- nesterov.cheby(xnew)
    # Delta.f.new <- nesterov.cheby.Delta(xnew, d)
    Delta.f.new <- nesterov.cheby.Delta(xnew, -1 * d)
    
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
      t <- 2 * t
    } else {
      t <- (alpha + beta) / 2
    }
  }
}

cvx.comp.bt <- function(x, d, c1) {
  # convex-composite backtracking scheme developed by Burke '85
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

weak.wolfe <- function(x, d, c1, c2) {
  # smooth weak Wolfe bisection
  # described in Overton-Lewis
  alpha <- 0
  beta <- Inf
  t <- 1
  
  f.data <- nesterov.cheby(x)
  
  while (TRUE) {
    xnew <- x + t * d
    fnew.data <- nesterov.cheby(xnew)
    
    if (fnew.data[[1]] > f.data[[1]] + c1 * t * t(f.data[[2]]) %*% d) {
      beta <- t
    } else if (c2 * t(f.data[[2]]) %*% d > t(fnew.data[[2]]) %*% d) {
      alpha <- t  
    } else {
      return(t)
    }
    
    if (beta == Inf) {
      t <- 2*t
    } else {
      t <- (alpha + beta) / 2
    }
  }
}

cvx.comp.wolfe.new <- function(x, d, c1, c2, lam) {
  # Experimental line search procedure based on cvx-comp Lagrangian
  alpha <- 0
  beta <- Inf
  t <- 1
  
  Delta.f <- nesterov.cheby.Delta.new(x, d, lam)
  while (TRUE) {
    # cat("Trying t = ", t, ".\n")
    xnew <- x + t * d
    Delta.f.new <- nesterov.cheby.Delta.new(xnew, d, lam)
    
    if (phi.n(xnew, lam) > phi.n(x, lam) + c1 * t * Delta.f) {
      beta <- t
    } else if (c2 * Delta.f > Delta.f.new) {
      alpha <- t  
      cat("Step is ", t, ".\n")
      cat("Have sufficient decrease.\n")
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