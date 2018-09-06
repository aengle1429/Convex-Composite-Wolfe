source("C:/Users/aengl/Dropbox/A.Engle/cvx_comp_R/nesterov_chebyshev/nesterov_functions.R")
cat("\f") # clear console
# set.seed(12)
n <- 10
xk <- rnorm(n, 0, 1000)
cat("Initializing at xk = ", xk)
invBk <- diag(n)
for(i in 1:50) {
  cat("iteration is ", i, "\n")
  
  fk.data <- nesterov.cheby(xk)
  dk <- -invBk %*% fk.data[[2]]
  cat("Delta is ", nesterov.cheby.Delta(xk, dk), "\n")
  
  # step <- cvx.comp.wolfe(xk, dk, 0.2, 0.8, w)
  step <- weak.wolfe(xk, dk, 0.0, 0.5)
  # cat("step is ", step, "\n")
  
  xnew <- xk + step * dk
  cat("xnew is ", xnew, "\n")
  
  fnew.data <- nesterov.cheby(xnew)
  
  cat("Function value is ", fnew.data[[1]], "\n")
  
  yk <- fnew.data[[2]] - fk.data[[2]]
  sk <- step * dk
  invBk <- BFGS.update(yk, sk, invBk)
  xk <- xnew
}
