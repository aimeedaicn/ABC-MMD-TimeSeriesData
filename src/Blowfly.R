init_ecological_model <- function(T = 180L) {
  
  simulate <- function(theta) {
    stopifnot(length(theta) == 6)
    logP       <- theta[1]
    logN0      <- theta[2]
    log_sd     <- theta[3]
    log_sp     <- theta[4]
    log_tau    <- theta[5]
    log_delta  <- theta[6]
    
    P       <- exp(logP)
    N0      <- exp(logN0)
    sigma_d <- exp(log_sd)
    sigma_p <- exp(log_sp)
    tau     <- max(1L, as.integer(round(exp(log_tau))))  
    delta   <- exp(log_delta)
    
    N <- numeric(T)
    N[1:max(2L, tau)] <- N0
    
    if (T >= tau + 1) {
      for (t in (tau):(T - 1L)) {
        e_t   <- rgamma(1L, shape = 1 / (sigma_p^2), scale = sigma_p^2)
        eps_t <- rgamma(1L, shape = 1 / (sigma_d^2), scale = sigma_d^2)
        N[t + 1L] <- P * N[t - tau + 1L] * exp(-N[t - tau + 1L] / N0) * e_t +
          N[t] * exp(-delta * eps_t)
      }
    }
    matrix(N, ncol = 1)
  }
  
  rprior <- function(N = 1L, hyperparams = NULL) {
    cbind(
      logP        = rnorm(N,  2, sqrt(4)),
      logN0       = rnorm(N,  5, sqrt(0.25)),
      log_sigma_d = rnorm(N, -0.5, 1),
      log_sigma_p = rnorm(N, -0.5, 1),
      log_tau     = rnorm(N,  2, 1),
      log_delta   = rnorm(N, -1, sqrt(0.16))
    )
  }
  
  dprior <- function(theta, hyperparams = NULL, sum = FALSE) {
    theta <- as.matrix(theta)
    if (ncol(theta) != 6) stop("theta must have 6 columns (logP, logN0, log_sigma_d, log_sigma_p, log_tau, log_delta)")
    lp <- dnorm(theta[,1],  2, sqrt(4),      log=TRUE) +
      dnorm(theta[,2],  5, sqrt(0.25),   log=TRUE) +
      dnorm(theta[,3], -0.5, 1,          log=TRUE) +
      dnorm(theta[,4], -0.5, 1,          log=TRUE) +
      dnorm(theta[,5],  2, 1,            log=TRUE) +
      dnorm(theta[,6], -1, sqrt(0.16),   log=TRUE)
    if (sum) sum(lp) else lp
  }
  
  list(simulate = simulate,
       rprior   = rprior,
       dprior   = dprior,
       info     = list(T = T))
}
