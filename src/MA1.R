init_ma1_model <- function(D = 1000L) {
  simulate <- function(theta) {
    stopifnot(length(theta) == 2)
    th <- theta[1]              
    ls <- theta[2]              
    if (abs(th) >= 1) stop("theta must be in (-1,1)")
    sg <- exp(ls)
    
    eps <- rnorm(D + 1L, mean = 0, sd = sg)
    y   <- eps[-1] + th * eps[-(D + 1L)]  
    matrix(y, ncol = 1)
  }
  
  rprior <- function(N = 1L, hyperparams = NULL) {
    cbind(
      theta    = runif(N, -0.99999, 0.99999),
      logsigma = rnorm(N, 0, 1)
    )
  }
  
  dprior <- function(theta, hyperparams = NULL, sum = FALSE) {
    theta <- as.matrix(theta)
    if (ncol(theta) != 2)
      stop("theta must have two columns (theta, logsigma)")
    th <- theta[, 1]
    ls <- theta[, 2]
    in_support   <- (th > -1 & th < 1)
    logpdf_theta <- ifelse(in_support, log(0.5), -Inf)    
    logpdf_lsg   <- dnorm(ls, mean = 0, sd = 1, log = TRUE)
    logpdf <- logpdf_theta + logpdf_lsg
    if (sum) sum(logpdf) else logpdf
  }
  
  list(simulate = simulate,
       rprior   = rprior,
       dprior   = dprior,
       info     = list(D = D))
}
