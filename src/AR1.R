init_ar1_model <- function(D = 200L) {
  simulate <- function(theta) {
    stopifnot(length(theta) == 2)
    phi      <- theta[1]       
    log_sig  <- theta[2]
    if (abs(phi) >= 1) stop("phi must be in (-1,1)")
    sigma    <- exp(log_sig)
    x <- numeric(D)
    x[1] <- rnorm(1, mean = 0, sd = sigma / sqrt(1 - phi^2))
    w <- rnorm(D - 1L, mean = 0, sd = 1)
    for (t in 1:(D - 1L)) {
      x[t + 1] <- phi * x[t] + sigma * w[t]
    }
    matrix(x, ncol = 1)
  }
  
  rprior <- function(N = 1L, hyperparams = NULL) {
    cbind(
      phi       = runif(N, -1, 1),
      logsigma  = rnorm(N, 0, 1)
    )
  }
  
  dprior <- function(theta, hyperparams = NULL, sum = FALSE) {
    theta <- as.matrix(theta)
    if (ncol(theta) != 2)
      stop("theta must have two columns (phi, logsigma)")
    phi      <- theta[, 1]
    logsigma <- theta[, 2]
    in_support <- (phi > -1 & phi < 1)
    logpdf_phi <- ifelse(in_support, log(0.5), -Inf)     
    logpdf_lsg <- dnorm(logsigma, mean = 0, sd = 1, log = TRUE)
    logpdf <- logpdf_phi + logpdf_lsg
    if (sum) sum(logpdf) else logpdf
  }
  
  list(simulate = simulate,
       rprior   = rprior,
       dprior   = dprior,
       info     = list(D = D))
}
