suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(gridExtra)
})

sample_from_grid_posterior <- function(theta1_g, theta2_g, post_mat, n = 5000) {
  stopifnot(is.matrix(post_mat))
  P <- post_mat / sum(post_mat)
  idx <- sample.int(length(P), size = n, replace = TRUE, prob = as.vector(P))
  i1 <- ((idx - 1L) %% nrow(post_mat)) + 1L
  i2 <- ((idx - 1L) %/% nrow(post_mat)) + 1L
  cbind(theta1 = theta1_g[i1], theta2 = theta2_g[i2])
}

# MMD
median_heuristic_sigma <- function(Z) {
  Z <- as.matrix(Z)
  D <- as.matrix(dist(Z))           
  median(D[lower.tri(D)], na.rm = TRUE)
}

rbf_kernel <- function(X, Y, sigma) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  XX <- rowSums(X * X); YY <- rowSums(Y * Y)
  sq <- outer(XX, YY, "+") - 2 * (X %*% t(Y))
  exp(-sq / (2 * sigma * sigma))
}

mmd_u_unbiased <- function(X, Y, sigma) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  m <- nrow(X); n <- nrow(Y)
  Kxx <- rbf_kernel(X, X, sigma)
  Kyy <- rbf_kernel(Y, Y, sigma)
  Kxy <- rbf_kernel(X, Y, sigma)
  
  diag(Kxx) <- 0
  diag(Kyy) <- 0
  term_xx <- sum(Kxx) / (m * (m - 1))
  term_yy <- sum(Kyy) / (n * (n - 1))
  term_xy <- (2 / (m * n)) * sum(Kxy)
  term_xx + term_yy - term_xy
}

# sliced 1-Wasserstein distance
oneD_w1 <- function(a, b) {
  na <- length(a); nb <- length(b)
  if (na == nb) {
    aa <- sort(a); bb <- sort(b)
    mean(abs(aa - bb))
  } else {
    Q <- min(na, nb)
    probs <- (seq_len(Q) - 0.5) / Q
    qa <- as.numeric(quantile(a, probs = probs, type = 1))
    qb <- as.numeric(quantile(b, probs = probs, type = 1))
    mean(abs(qa - qb))
  }
}

sliced_w1 <- function(X, Y, n_projections = 200L, seed = 123) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  d <- ncol(X); stopifnot(ncol(Y) == d)
  set.seed(seed)
  U <- matrix(rnorm(d * n_projections), nrow = d)
  U <- sweep(U, 2, sqrt(colSums(U^2)), "/")  
  vals <- vapply(seq_len(ncol(U)), function(j) {
    u <- U[, j, drop = FALSE]
    oneD_w1(as.numeric(X %*% u), as.numeric(Y %*% u))
  }, numeric(1))
  mean(vals)
}

assess_posteriors <- function(true_samples,
                              samples_by_method,
                              n_projections = 200L,
                              sigma = c("median", "fixed"),
                              sigma_value = NULL) {
  sigma <- match.arg(sigma)
  rows <- list()
  for (meth in names(samples_by_method)) {
    obj <- samples_by_method[[meth]]
    reps <- if (is.list(obj)) obj else list(obj)
    for (r in seq_along(reps)) {
      X <- as.matrix(reps[[r]])
      Y <- as.matrix(true_samples)
      sig <- if (sigma == "median") {
        median_heuristic_sigma(rbind(X, Y))
      } else {
        if (is.null(sigma_value) || !is.finite(sigma_value) || sigma_value <= 0)
          stop("Provide a positive sigma_value when sigma='fixed'.")
        sigma_value
      }
      w1  <- sliced_w1(X, Y, n_projections = n_projections)
      mmd <- mmd_u_unbiased(X, Y, sigma = sig)
      rows[[length(rows) + 1L]] <- tibble(
        method = meth, rep = r, W1 = w1, MMD2 = mmd
      )
    }
  }
  bind_rows(rows)
}

plot_assessment_boxplots <- function(metrics_df, out_file = NULL) {
  df_long <- metrics_df |>
    tidyr::pivot_longer(cols = c(W1, MMD2), names_to = "metric", values_to = "value")
  
  p_w1 <- ggplot(df_long %>% filter(metric == "W1"),
                 aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    labs(x = "", y = expression(W[1](pi[ABC], pi))) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1),
          panel.border = element_rect(colour = "black", fill = NA))
  
  p_mmd <- ggplot(df_long %>% filter(metric == "MMD2"),
                  aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    labs(x = "", y = expression(MMD[U](pi[ABC], pi))) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1),
          panel.border = element_rect(colour = "black", fill = NA))
  
  if (!is.null(out_file)) {
    g <- gridExtra::grid.arrange(p_w1, p_mmd, nrow = 1)
    ggsave(out_file, g, width = 12, height = 5)
  }
  list(p_w1 = p_w1, p_mmd = p_mmd)
}