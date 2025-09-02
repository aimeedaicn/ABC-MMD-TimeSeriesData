init_ou_model <- function(D       = 50L,      # number of time steps
                          delta_t = 0.2,
                          x0      = 10,
                          sigma_w = 0.5) {
  # Trajectory simulator
  simulate <- function(theta) {
    stopifnot(length(theta) == 2)
    theta1 <- theta[1]      # must be in (0,1)
    theta2 <- theta[2]      # must be in (-2,2)
    
    x <- numeric(D)
    x[1] <- x0
    
    eps <- rnorm(D - 1L, mean = 0, sd = sqrt(delta_t))
    
    for (t in 1:(D - 1L)) {
      drift  <- theta1 * (exp(theta2) - x[t]) * delta_t
      diff   <- sigma_w * eps[t]
      x[t+1] <- x[t] + drift + diff
    }
    matrix(x, ncol = 1)
  }
  
  # Prior sampler 
  rprior <- function(N = 1L, hyperparams = NULL) {
    cbind(
      theta1 = runif(N, 0, 1),
      theta2 = runif(N, -2, 2)
    )
  }
  
  # Prior log-density 
  dprior <- function(theta, hyperparams = NULL, sum = FALSE) {
    theta <- as.matrix(theta)
    if (ncol(theta) != 2)
      stop("theta must have two columns (theta1, theta2)")
    
    in_support <- (theta[, 1] > 0  & theta[, 1] < 1) &
      (theta[, 2] > -2 & theta[, 2] < 2)
    
    logpdf <- ifelse(in_support,
                     log(1 / 1) + log(1 / 4),   # Uniform PDF: 1/(b-a)
                     -Inf)
    if (sum) sum(logpdf) else logpdf
  }
  
  # loglik_ou <- function(theta, x) {
  #   theta1 <- theta[1]; theta2 <- theta[2]
  #   resid  <- x[-1] - x[-length(x)] -
  #     theta1 * (exp(theta2) - x[-length(x)]) * delta_t
  #   sigma2 <- sigma_w^2 * delta_t
  #   -0.5 * length(resid) * log(2 * pi * sigma2) -
  #     sum(resid^2) / (2 * sigma2)
  # }
  
  # Return list 
  list(simulate = simulate,
       rprior   = rprior,
       dprior   = dprior,
       info     = list(D = D, delta_t = delta_t, x0 = x0, sigma_w = sigma_w))
}

