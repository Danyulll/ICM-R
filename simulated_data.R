###############################################################################
# SIMULATED DATA GENERATION AND ANOMALY TYPES
# ---------------------------------------------------------------------------
# This file contains all functions related to data generation, anomaly types,
# and dataset creation for the semi-supervised DP-ICM-GP model.
###############################################################################

# --------------------------- Data generation (ICM) ------------------------------

gen_icm_curves <- function(N = 50, P = 128, M = 3, t = seq(0, 1, len = P),
                           B = matrix(c(1, .7, .5, .7, 1, .4, .5, .4, 1), 3, 3),
                           eta = rep(0.02, 3),
                           kern = kernels[[1]], par = list(l_scale = 0.2)) {
  # --- NEW: make sure B is a proper MxM matrix ---
  if (!is.matrix(B)) {
    if (length(B) == 1) {
      B <- matrix(B, M, M)
    } else if (length(B) == M) {
      B <- diag(as.numeric(B), M, M)
    } else if (length(B) == M * M) {
      B <- matrix(as.numeric(B), M, M, byrow = FALSE)
    } else {
      stop(sprintf("B must be coercible to %dx%d; got length=%d", M, M, length(B)))
    }
  }
  if (nrow(B) != ncol(B) || nrow(B) != M) {
    stop(sprintf("B must be %dx%d; got %dx%d", M, M, nrow(B), ncol(B)))
  }
  # ------------------------------------------------

  Kx <- kern$fun(t, par)
  K <- kronecker_icm(B, Kx, eta)
  Y <- vector("list", N)
  for (i in 1:N) {
    yvec <- as.numeric(MASS::mvrnorm(1, mu = rep(0, P * M), Sigma = K + diag(1e-6, P * M)))
    Y[[i]] <- matrix(yvec, nrow = P, ncol = M)
  }
  Y
}

# ======================= Anomaly generators & dataset ========================

.pick_win <- function(P, min_len = max(4, ceiling(P * 0.05)), max_len = max(6, ceiling(P * 0.20))) {
  L <- sample(min_len:max_len, 1)
  s <- sample(1:(P - L + 1), 1)
  list(s = s, e = s + L - 1, L = L)
}

.pick_subset <- function(M, min_k = 1) {
  k <- sample(min_k:M, 1)
  sort(sample(1:M, k))
}

.hann <- function(n) {
  if (n <= 1) {
    return(rep(1, n))
  }
  0.5 - 0.5 * cos(2 * pi * (0:(n - 1)) / (n - 1))
}


anom_spikes <- function(y, k = 5, amp = NULL) { # Increased k from 3 to 5
  P <- nrow(y)
  M <- ncol(y)
  if (is.null(amp)) amp <- 5 * sd(as.vector(y)) # Increased from 3x to 5x
  for (u in 1:k) {
    t0 <- sample(1:P, 1)
    m <- sample(1:M, 1)
    y[t0, m] <- y[t0, m] + rnorm(1, 0, amp)
  }
  y
}

anom_burst <- function(y, amp = NULL, f = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  if (is.null(amp)) amp <- 1.5 * sd(as.vector(y))
  if (is.null(f)) f <- runif(1, 4, 12)
  tt <- seq(0, 1, length.out = w$L)
  msel <- sample(1:M, sample(1:M, 1))
  for (m in msel) y[w$s:w$e, m] <- y[w$s:w$e, m] + amp * sin(2 * pi * f * tt)
  y
}

anom_varburst <- function(y, factor = 5) { # Increased from 3 to 5
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  add <- matrix(rnorm(w$L * M, 0, factor * sd(as.vector(y))), nrow = w$L, ncol = M)
  y[w$s:w$e, ] <- y[w$s:w$e, ] + add
  y
}

anom_step <- function(y, delta = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  if (is.null(delta)) delta <- rnorm(M, 0, 4 * apply(y, 2, sd)) # Increased from 2x to 4x
  y[w$s:w$e, ] <- sweep(y[w$s:w$e, ], 2, delta, "+")
  y
}

anom_ramp <- function(y, slope = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  if (is.null(slope)) slope <- rnorm(M, 0, 3 * apply(y, 2, sd) / w$L)
  ramp <- outer(seq(0, 1, length.out = w$L), slope, "*")
  y[w$s:w$e, ] <- y[w$s:w$e, ] + ramp
  y
}

anom_freqchange <- function(y, amp = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  if (is.null(amp)) amp <- apply(y, 2, sd)
  f_new <- runif(1, 6, 14)
  tt <- seq(0, 1, length.out = w$L)
  msel <- sample(1:M, sample(1:M, 1))
  for (m in msel) y[w$s:w$e, m] <- y[w$s:w$e, m] + amp[m] * sin(2 * pi * f_new * tt)
  y
}

anom_phase <- function(y, shift = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  if (is.null(shift)) shift <- runif(1, 0.2, 0.8)
  Ls <- floor(shift * w$L)
  for (m in 1:M) {
    seg <- y[w$s:w$e, m]
    if (Ls > 0 && Ls < length(seg)) y[w$s:w$e, m] <- c(seg[(Ls + 1):length(seg)], seg[1:Ls])
  }
  y
}

anom_corrbreak <- function(y, theta = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  if (M < 2) {
    return(y)
  }
  w <- .pick_win(P)
  if (is.null(theta)) theta <- runif(1, pi / 6, pi / 3)
  R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  idx <- 1:min(2, M)
  y[w$s:w$e, idx] <- as.matrix(y[w$s:w$e, idx]) %*% R
  y
}

anom_warp <- function(y, factor = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  if (is.null(factor)) factor <- runif(1, 0.5, 1.5)
  x <- seq(0, 1, length.out = w$L)
  x_new <- pmin(1, pmax(0, x^factor))
  for (m in 1:M) {
    seg <- y[w$s:w$e, m]
    y[w$s:w$e, m] <- approx(x, seg, xout = x_new, rule = 2)$y
  }
  y
}

anom_dropout <- function(y) {
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  y[w$s:w$e, ] <- matrix(y[max(1, w$s - 1), ], nrow = w$L, ncol = M, byrow = TRUE)
  y
}

anom_swap <- function(y) {
  P <- nrow(y)
  M <- ncol(y)
  if (M < 2) {
    return(y)
  }
  w <- .pick_win(P)
  ch <- sample(1:M, 2)
  tmp <- y[w$s:w$e, ch[1]]
  y[w$s:w$e, ch[1]] <- y[w$s:w$e, ch[2]]
  y[w$s:w$e, ch[2]] <- tmp
  y
}

anom_global_scale <- function(y, factor = NULL) {
  if (is.null(factor)) factor <- runif(1, 1.5, 2.2)
  y * factor
}

# -------- Vertical translation anomalies (pure additive offsets) --------------

# (a) Global constant offset (all channels, all times)
anom_vglob <- function(y, delta = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  if (is.null(delta)) delta <- 2.0 * sd(as.vector(y)) # scale knob
  y + delta
}

# (b) Per-channel constant offsets (all times)
anom_vchan <- function(y, delta = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  if (is.null(delta)) delta <- rnorm(M, 0, 2.0 * apply(y, 2, sd))
  sweep(y, 2, delta, "+")
}

# (c) Group-shared constant offset (subset of channels)
anom_vgroup <- function(y, delta = NULL, min_k = 1) {
  P <- nrow(y)
  M <- ncol(y)
  S <- .pick_subset(M, min_k = min_k)
  if (is.null(delta)) delta <- 2.0 * mean(apply(y[, S, drop = FALSE], 2, sd))
  y[, S] <- y[, S, drop = FALSE] + delta
  y
}

# (d) Time-localized step (plateau) per channel
anom_vstep <- function(y, delta = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  if (is.null(delta)) delta <- rnorm(M, 0, 2.0 * apply(y, 2, sd))
  y[w$s:w$e, ] <- sweep(y[w$s:w$e, , drop = FALSE], 2, delta, "+")
  y
}

# (e) Windowed smooth bump (Hann window) per channel
anom_vbump <- function(y, delta = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  w <- .pick_win(P)
  if (is.null(delta)) delta <- rnorm(M, 0, 2.0 * apply(y, 2, sd))
  win <- .hann(w$L)
  for (m in 1:M) y[w$s:w$e, m] <- y[w$s:w$e, m] + delta[m] * win
  y
}

# (f) Correlated multichannel offset (fixed direction u)
anom_vdir <- function(y, delta = NULL, u = NULL) {
  P <- nrow(y)
  M <- ncol(y)
  if (is.null(u)) {
    u <- rnorm(M)
    u <- u / sqrt(sum(u^2) + 1e-12)
  } else {
    u <- as.numeric(u)
    u <- u / sqrt(sum(u^2) + 1e-12)
  }
  if (is.null(delta)) delta <- 2.0 * sd(as.vector(y))
  y + matrix(delta * rep(u, each = P), nrow = P, ncol = M)
}


ANOM_REG <- list(
  spikes = anom_spikes,
  burst = anom_burst,
  varb = anom_varburst,
  step = anom_step,
  ramp = anom_ramp,
  fchg = anom_freqchange,
  phase = anom_phase,
  corr = anom_corrbreak,
  warp = anom_warp,
  drop = anom_dropout,
  swap = anom_swap,
  gsc = anom_global_scale, vglob = anom_vglob,
  vchan = anom_vchan,
  vgroup = anom_vgroup,
  vstep = anom_vstep,
  vbump = anom_vbump,
  vdir = anom_vdir
)

make_anomaly_dataset <- function(anom_type = NULL, N = 200, P = 128, M = 3, t = seq(0, 1, len = P),
                                 base_B = diag(M), base_eta = rep(0.02, M),
                                 base_kern = kernels[[1]], base_par = list(l_scale = 0.2),
                                 anom_rates = NULL, frac_anom = 0.05, seed = 1) {
  set.seed(seed)
  Y <- gen_icm_curves(N = N, P = P, M = M, t = t, B = base_B, eta = base_eta, kern = base_kern, par = base_par)
  y_all <- Y
  ylab <- rep(0L, N) # 0=normal, >0 anomalous type id
  n_anom <- floor(frac_anom * N)
  idx_anom <- if (n_anom > 0) sort(sample(1:N, n_anom)) else integer(0)

  types <- names(ANOM_REG)
  if (!is.null(anom_type)) {
    # If specific anomaly type is requested, use only that type
    if (!anom_type %in% types) {
      stop(sprintf("Unknown anomaly type: %s. Available types: %s", anom_type, paste(types, collapse = ", ")))
    }
    w <- rep(1e-6, length(types))
    names(w) <- types
    w[anom_type] <- 1.0
  } else if (is.null(anom_rates)) {
    w <- rep(1 / length(types), length(types))
    names(w) <- types
  } else {
    stopifnot(all(names(anom_rates) %in% types))
    w <- rep(1e-6, length(types))
    names(w) <- types
    w[names(anom_rates)] <- anom_rates
    w <- w / sum(w)
  }
  for (i in idx_anom) {
    ty <- sample(types, 1, prob = w)
    y_all[[i]] <- ANOM_REG[[ty]](y_all[[i]])
    ylab[i] <- match(ty, types)
  }
  true_anom <- as.integer(seq_len(N) %in% idx_anom)
  list(
    Y = y_all,
    labels = ylab,
    true_anom = true_anom,
    types = types,
    idx_anom = idx_anom,
    t = t # <-- include time grid
  )
}

# ------------------------ Multi-cluster anomaly detection logic ------------

# Function to determine normal vs anomaly clusters when there are multiple clusters
determine_normal_anomaly_clusters <- function(z_hat, reveal_idx = integer(0), min_clusters_for_combination = 3) {
  unique_clusters <- unique(z_hat)
  n_clusters <- length(unique_clusters)

  # If we have fewer than min_clusters_for_combination clusters, use original logic
  if (n_clusters < min_clusters_for_combination) {
    if (length(reveal_idx) > 0) {
      revealed_clusters <- z_hat[reveal_idx]
      tab_revealed <- table(revealed_clusters)
      normal_cluster <- as.integer(names(which.max(tab_revealed)))
    } else {
      tabz <- table(z_hat)
      normal_cluster <- as.integer(names(which.max(tabz)))
    }
    pred_anom <- as.integer(z_hat != normal_cluster)
    return(list(
      pred_anom = pred_anom, normal_cluster = normal_cluster,
      anomaly_clusters = setdiff(unique_clusters, normal_cluster)
    ))
  }

  # Multi-cluster logic: cluster with most revealed normals is normal, others are anomaly
  if (length(reveal_idx) > 0) {
    revealed_clusters <- z_hat[reveal_idx]
    tab_revealed <- table(revealed_clusters)
    normal_cluster <- as.integer(names(which.max(tab_revealed)))
  } else {
    # If no revealed normals, use the largest cluster as normal
    tabz <- table(z_hat)
    normal_cluster <- as.integer(names(which.max(tabz)))
  }

  # All other clusters are considered anomalies
  anomaly_clusters <- setdiff(unique_clusters, normal_cluster)
  pred_anom <- as.integer(z_hat != normal_cluster)

  return(list(
    pred_anom = pred_anom, normal_cluster = normal_cluster,
    anomaly_clusters = anomaly_clusters
  ))
}

# Plotting functions are now in plotting_functions.R
