init_cosine_model <- function(n        = 100L,    
                              omega_lo = 0, omega_hi = 1/10,
                              phi_lo   = 0, phi_hi = 2*pi) {
  simulate <- function(theta) {
    stopifnot(length(theta) == 4)
    omega      <- theta[1]
    phi        <- theta[2]
    log_sigma  <- theta[3]
    log_A      <- theta[4]
    sigma <- exp(log_sigma)
    A     <- exp(log_A)
    t_seq <- seq_len(n)
    mu    <- A * cos(2*pi*omega * t_seq + phi)
    y     <- mu + sigma * rnorm(n)
    matrix(y, ncol = 1)
  }
  
  rprior <- function(N = 1L, hyperparams = NULL) {
    cbind(
      omega     = runif(N, omega_lo, omega_hi),
      phi       = runif(N, phi_lo,   phi_hi),
      logsigma  = rnorm(N, 0, 1),   
      logA      = rnorm(N, 0, 1)   
    )
  }
  
  dprior <- function(theta, hyperparams = NULL, sum = FALSE) {
    theta <- as.matrix(theta)
    if (ncol(theta) != 4) stop("theta must have 4 columns (omega, phi, log_sigma, log_A)")
    in_support <- (theta[,1] >= omega_lo & theta[,1] <= omega_hi) &
      (theta[,2] >= phi_lo   & theta[,2] <= phi_hi)
    logp_uni <- ifelse(in_support,
                       log(1/(omega_hi-omega_lo)) + log(1/(phi_hi-phi_lo)),
                       -Inf)
    logp_norm <- dnorm(theta[,3], 0, 1, log = TRUE) + dnorm(theta[,4], 0, 1, log = TRUE)
    out <- logp_uni + logp_norm
    if (sum) sum(out) else out
  }
  
  list(simulate = simulate,
       rprior   = rprior,
       dprior   = dprior,
       info     = list(n = n,
                       omega_range = c(omega_lo, omega_hi),
                       phi_range   = c(phi_lo,   phi_hi)))
}
