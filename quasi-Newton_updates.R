# 
# Functions for computing updates to quasi-Newton update methods
# including the BFGS, SR1, DFP, and damped versions
# 
damped.BFGS.update <- function(yk, dk, tk, Bk) {
  yksk <- as.numeric(crossprod(yk, tk * dk))
  vec1 <- as.vector(Bk %*% dk)
  skBksk <- tk^2 * as.numeric(crossprod(dk, vec1))
  # cat("Error here... ", yksk, "\n", skbksk, "\n")
  if (yksk >= 0.2 * skBksk) {
    cat("Not updating thetak.\n")
    thetak <- 1
  } else {
    thetak <- (0.8 * skBksk) / (skBksk - yksk)  
  }
  
  rk <- thetak * yk + (1 - thetak) * tk * vec1
  return(Bk - tk^2 * vec1 %o% vec1 / skBksk +  as.vector(rk)%o%as.vector(rk) / as.numeric(crossprod(tk * dk, rk)))
}

damped.DFP.update <- function(yk, dk, tk, Bk) {
  yksk <- as.numeric(crossprod(yk, tk * dk))
  vec1 <- as.vector(Bk %*% dk)
  skBksk <- tk^2 * as.numeric(crossprod(dk, vec1))
  # cat("Error here... ", yksk, "\n", skbksk, "\n")
  if (yksk >= 0.2 * skBksk) {
    cat("Not updating thetak.\n")
    thetak <- 1
  } else {
    thetak <- (0.8 * skBksk) / (skBksk - yksk)  
  }
  
  rk <- thetak * yk + (1 - thetak) * tk * vec1
  return(Bk - tk^2 * vec1 %o% vec1 / skBksk +  as.vector(rk)%o%as.vector(rk) / as.numeric(crossprod(tk * dk, rk)))
}


SR1.update <- function(yk, dk, tk, Bk) {
  vec1 <- as.vector(yk - Bk %*% (tk * dk))
  if (abs(as.numeric(crossprod(vec1, dk * tk))) >= 1e-8 * sqrt(sum((dk * tk)^2)) * sqrt(sum(vec1^2))) {
    return(Bk + vec1 %o% vec1 / as.numeric(crossprod(vec1, dk * tk)))
  } else {
    return(Bk)
  }
}

BFGS.update <- function(yk, dk, tk, Bk) {
  yksk <- as.numeric(crossprod(yk, tk*dk))
  stopifnot(yksk > 0)
  n <- length(yk)
  vec1 <- as.vector(Bk %*% dk)
  return(Bk - vec1 %o% vec1 / as.numeric(crossprod(dk, vec1)) +  as.vector(yk)%o%as.vector(yk) / yksk)
}

DFP.update <- function(yk, dk, tk, Bk) {
  yksk <- as.numeric(crossprod(yk, tk*dk))
  n <- length(yk)
  
  mat1 <- diag(n) - tcrossprod(yk, dk) / as.numeric(crossprod(yk, dk))
  B <- mat1 %*% Bk %*% t(mat1) + tcrossprod(yk, yk) / yksk
  return((B + t(B))/ 2)
}
