###############################################################################
# MCMC SAMPLING FUNCTIONS AND MODEL COMPONENTS
# ---------------------------------------------------------------------------
# This file contains all MCMC sampling functions, model components, utilities,
# and diagnostic functions for the semi-supervised DP-ICM-GP model.
###############################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(invgamma)
  library(waveslim)
  library(coda)
  library(ROCR)
})

# --------------------------- Utilities ---------------------------------------

ensure_dyadic_J <- function(P, J) {
  if (is.null(J)) J <- log2(P)
  J_int <- as.integer(round(J))
  if (abs(P - 2^J_int) > .Machine$double.eps * max(1, P)) {
    stop(sprintf("P must be 2^J (dyadic). Got P=%s, J≈%.6f (rounded J=%d gives 2^J=%s).",
                 P, J, J_int, 2^J_int))
  }
  J_int
}

logit   <- function(x) log(x/(1-x))
invlogit<- function(z) 1/(1+exp(-z))

chol_loglik <- function(y, Sigma, jitter=1e-6) {
  p <- nrow(Sigma)
  ok <- FALSE; tries <- 0; L <- NULL; J <- jitter
  while(!ok && tries < 7) {
    tries <- tries + 1
    out <- try(chol(Sigma + diag(J, p)), silent=TRUE)
    if(!inherits(out, "try-error")) { L <- out; ok <- TRUE } else J <- J*10
  }
  if(!ok) return(list(ll = -Inf, L = NULL, jitter = NA_real_))
  # evaluate *with the same PD matrix* that succeeded
  Sigma_pd <- tcrossprod(L)  # equals Sigma + diag(J, p)
  ll <- mvtnorm::dmvnorm(y, sigma = Sigma_pd, log = TRUE)
  list(ll = ll, L = L, jitter = J)
}

kronecker_icm <- function(B, Kx, eta) {
  # Coerce B to a proper square matrix if possible
  if (!is.matrix(B)) {
    # infer M from eta if available
    M_guess <- if (is.numeric(eta) && length(eta) > 0) length(eta) else NULL
    if (!is.null(M_guess)) {
      if (length(B) == 1) {
        B <- matrix(B, M_guess, M_guess)
      } else if (length(B) == M_guess) {
        B <- diag(as.numeric(B), M_guess, M_guess)
      } else if (length(B) == M_guess * M_guess) {
        B <- matrix(as.numeric(B), M_guess, M_guess, byrow = FALSE)
      }
    }
  }
  if (!is.matrix(B) || nrow(B) != ncol(B)) stop("B must be square.")
  if (!is.matrix(Kx) || nrow(Kx) != ncol(Kx)) stop("Kx must be square.")
  M <- nrow(B); P <- nrow(Kx)
  if (length(eta) != M) stop(sprintf("eta length (%d) must match nrow(B)=%d.", length(eta), M))
  K <- kronecker(B, Kx)
  K <- K + kronecker(diag(eta, M, M), diag(P))
  K
}



pack_L <- function(L) L[lower.tri(L, diag=TRUE)]
unpack_L <- function(theta, m) {
  L <- matrix(0, m, m)
  L[lower.tri(L, diag=TRUE)] <- theta
  diag(L) <- abs(diag(L)) + 1e-8
  L
}

# --------------------------- Kernel families (self-contained) ---------------

make_kernels <- function() {
  # helpers live INSIDE the closure -> serialized to workers
  k_sqexp <- function(t, l_scale) {
    D2 <- as.matrix(dist(t))^2
    exp(-0.5 * D2 / (l_scale^2))
  }
  k_mat32 <- function(t, l_scale) {
    D <- as.matrix(dist(t)); r <- D / l_scale; a <- sqrt(3) * r
    (1 + a) * exp(-a)
  }
  k_mat52 <- function(t, l_scale) {
    D <- as.matrix(dist(t)); r <- D / l_scale; a <- sqrt(5) * r
    (1 + a + 5*r^2/3) * exp(-a)
  }
  k_periodic <- function(t, l_scale, period) {
    D <- as.matrix(dist(t))
    exp( - 2 * sin(pi * D / period)^2 / (l_scale^2) )
  }
  
  list(
    list(
      name="SE", fun=function(t, par) k_sqexp(t, par$l_scale),
      pnames=c("l_scale"),
      prior = function(par) stats::dgamma(par$l_scale, 2, 2, log=TRUE),
      pstar = function() list(l_scale = stats::rgamma(1, 2, 2)),
      prop_sd = list(l_scale=0.20)
    ),
    list(
      name="Mat32", fun=function(t, par) k_mat32(t, par$l_scale),
      pnames=c("l_scale"),
      prior = function(par) stats::dgamma(par$l_scale, 2, 2, log=TRUE),
      pstar = function() list(l_scale = stats::rgamma(1, 2, 2)),
      prop_sd = list(l_scale=0.20)
    ),
    list(
      name="Mat52", fun=function(t, par) k_mat52(t, par$l_scale),
      pnames=c("l_scale"),
      prior = function(par) stats::dgamma(par$l_scale, 2, 2, log=TRUE),
      pstar = function() list(l_scale = stats::rgamma(1, 2, 2)),
      prop_sd = list(l_scale=0.20)
    ),
    list(
      name="Periodic", fun=function(t, par) k_periodic(t, par$l_scale, par$period),
      pnames=c("l_scale","period"),
      prior = function(par) stats::dgamma(par$l_scale, 3, 2, log=TRUE) + stats::dbeta(par$period, 5, 5, log=TRUE),
      pstar = function() list(l_scale = stats::rgamma(1, 3, 2), period = stats::rbeta(1, 5, 5)),
      prop_sd = list(l_scale=0.20, period=0.20)
    )
  )
}

# ------------------------ Wavelet wrappers (waveslim) ------------------------

wt_inverse_mat_keep <- function(y_mat, wtf_list, gamma_list) {
  P <- nrow(y_mat); M <- ncol(y_mat)
  y_norm <- matrix(0, P, M)
  for (m in 1:M) {
    cvec <- wtf_list[[m]]$coeff
    gvec <- gamma_list[[m]]
    cvec[gvec == 0] <- 0
    y_norm[,m] <- wt_inverse_1d(cvec, wtf_list[[m]]$map)
  }
  y_norm
}

wt_stack_channel <- function(Y_list, wf="la8", J=NULL, boundary="periodic") {
  M <- ncol(Y_list[[1]]); N <- length(Y_list)
  tmp <- wt_forward_mat(Y_list[[1]], wf=wf, J=J, boundary=boundary)
  ncoeff <- length(tmp[[1]]$coeff)
  D_arr <- array(NA_real_, dim=c(ncoeff, N, M))
  for (i in 1:N) {
    wtf <- wt_forward_mat(Y_list[[i]], wf=wf, J=J, boundary=boundary)
    for (m in 1:M) D_arr[,i,m] <- wtf[[m]]$coeff
  }
  list(D_arr=D_arr, maps=tmp)  # maps per channel (length M)
}

# ------------------- Besov schedules (optional, simple) ----------------------

besov_pi_schedule <- function(J, b = 0.5, c2 = 0.8) {
  setNames(pmin(1, c2 * 2^(-b * (1:J))), paste0("d", 1:J))
}

besov_g_hyper <- function(J, a = 0.5, base_shape = 2.0, base_rate = 2.0) {
  data.frame(
    level = paste0("d", 1:J),
    shape = base_shape + 0.25 * (1:J),
    rate  = base_rate  * 2^(a * (1:J)),
    row.names = paste0("d", 1:J)
  )
}

# ---------------- Ray–Mallick wavelet block with sigma_i^2 (hetero) ----------

update_wavelet_params_rm <- function(Y_list, wf, J, boundary,
                                     wpar, sigma2_i,
                                     a_pi=1, b_pi=1,
                                     g_hyp=NULL,
                                     a_sig=2.5, b_sig=0.02,
                                     use_besov_pi=FALSE, pi_sched=NULL) {
  stopifnot(!is.null(J))
  M <- ncol(Y_list[[1]]); N <- length(Y_list)
  
  st  <- wt_stack_channel(Y_list, wf=wf, J=J, boundary=boundary)
  D   <- st$D_arr
  maps<- st$maps
  
  ncoeff <- dim(D)[1]
  lev_names <- names(maps[[1]]$map$idx)
  det_names <- lev_names[grepl("^d", lev_names)]
  s_name    <- lev_names[grepl("^s", lev_names)]
  
  if (is.null(wpar$pi_level))  wpar$pi_level <- setNames(rep(0.5, length(det_names)), det_names)
  if (is.null(wpar$g_level))   wpar$g_level  <- setNames(rep(2.0, length(det_names)), det_names)
  if (is.null(wpar$gamma_ch))  wpar$gamma_ch <- lapply(1:M, function(m) rbinom(ncoeff,1,0.2))
  if (length(sigma2_i) != N)   sigma2_i <- rep(ifelse(is.finite(mean(sigma2_i)), mean(sigma2_i), 0.05), N)
  
  if (use_besov_pi && !is.null(pi_sched)) {
    wpar$pi_level[names(pi_sched)] <- pi_sched
  }
  if (!is.null(g_hyp)) {
    if (!all(names(wpar$g_level) %in% rownames(g_hyp))) stop("g_hyp must index detail levels.")
  }
  
  # 1) Update gamma
  for (m in 1:M) {
    Dm <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
    gam <- wpar$gamma_ch[[m]]
    for (lev in det_names) {
      ids <- maps[[m]]$map$idx[[lev]]
      if (length(ids) == 0) next
      pi_j <- wpar$pi_level[[lev]]
      g_j  <- wpar$g_level[[lev]]
     v_spike <- matrix(rep(sigma2_i, each=length(ids)), nrow=length(ids), ncol=N)
     v_slab  <- (1 + g_j) * v_spike
      Dsub    <- Dm[ids, , drop=FALSE]
     log_like_spike <- -0.5 * rowSums(log(2*pi*v_spike) + (Dsub^2)/v_spike)
     log_like_slab  <- -0.5 * rowSums(log(2*pi*v_slab ) + (Dsub^2)/v_slab )
      logit_val <- log(pi_j) + log_like_slab - (log(1-pi_j) + log_like_spike)
      p1        <- plogis(pmax(pmin(logit_val, 35), -35))
      gam[ids]  <- rbinom(length(ids), 1, p1)
    }
    if (length(s_name) == 1) {
      ids_s <- maps[[m]]$map$idx[[s_name]]
      if (length(ids_s) > 0) gam[ids_s] <- 1L
    }
    wpar$gamma_ch[[m]] <- gam
  }
  
  # 2) Update g_level
  for (lev in det_names) {
    ss_over_sigma <- 0
    n_sel_total   <- 0
    for (m in 1:M) {
      ids <- maps[[m]]$map$idx[[lev]]
      if (length(ids) == 0) next
      Dm  <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
      sel <- wpar$gamma_ch[[m]][ids] == 1
      if (any(sel)) {
        Did <- Dm[ids[sel], , drop=FALSE]
        ss_over_sigma <- ss_over_sigma + sum( sweep(Did^2, 2, sigma2_i, "/") )
        n_sel_total   <- n_sel_total + nrow(Did)
      }
    }
    shape0 <- if (!is.null(g_hyp)) g_hyp[lev, "shape"] else 2.0
    rate0  <- if (!is.null(g_hyp)) g_hyp[lev, "rate"]  else 2.0
    shape_post <- shape0 + 0.5 * n_sel_total
    rate_post  <- rate0  + 0.5 * ss_over_sigma
    wpar$g_level[[lev]] <- invgamma::rinvgamma(1, shape=shape_post, rate=rate_post)
    
  }
  
  # 3) Update pi_level
  for (lev in det_names) {
    n1 <- 0; n0 <- 0
    for (m in 1:M) {
      ids <- maps[[m]]$map$idx[[lev]]
      if (length(ids) == 0) next
      gm <- wpar$gamma_ch[[m]][ids]
      n1 <- n1 + sum(gm == 1)
      n0 <- n0 + sum(gm == 0)
    }
    if (!use_besov_pi) wpar$pi_level[[lev]] <- rbeta(1, a_pi + n1, b_pi + n0)
  }
  
  # 4) Update sigma2_i
  for (i in 1:N) {
    ss <- 0
    n_eff <- 0
    for (m in 1:M) {
      Dm <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
      gm <- wpar$gamma_ch[[m]]
      for (lev in det_names) {
        ids <- maps[[m]]$map$idx[[lev]]
        if (length(ids) == 0) next
        dvec <- Dm[ids, i]
        g_j  <- wpar$g_level[[lev]]
        sel  <- gm[ids] == 1
        if (any(sel)) {
          ss <- ss + sum(dvec[sel]^2   / (1 + g_j))
          n_eff <- n_eff + sum(sel)
        }
        if (any(!sel)) {
         ss <- ss + sum(dvec[!sel]^2  / 1.0)
          n_eff <- n_eff + sum(!sel)
        }
      }
      if (length(s_name) == 1) {
        ids_s <- maps[[m]]$map$idx[[s_name]]
        if (length(ids_s) > 0) {
          dvec <- Dm[ids_s, i]
          ss   <- ss + sum(dvec^2)
          n_eff <- n_eff + length(dvec)
        }
      }
    }
    sigma2_i[i] <- invgamma::rinvgamma(1, shape = a_sig + 0.5 * n_eff, rate = b_sig + 0.5 * ss)
  }
  
  list(wpar=wpar, sigma2_i=sigma2_i, maps=maps)
}

reconstruct_normals_for_cluster <- function(Y_list, wf, J, boundary, wpar) {
  N <- length(Y_list); P <- nrow(Y_list[[1]]); M <- ncol(Y_list[[1]])
  out <- vector("list", N)
  for (i in 1:N) {
    wtf <- wt_forward_mat(Y_list[[i]], wf=wf, J=J, boundary=boundary)
    gam <- wpar$gamma_ch
    out[[i]] <- wt_inverse_mat_keep(Y_list[[i]], wtf, gam)
  }
  out
}

# ================== Cluster-level Ray–Mallick wavelet updater =================

update_cluster_wavelet_params_rm <- function(
  Y_list, wf, J, boundary,
  wpar, sigma2,
  a_pi = 1, b_pi = 1,         # Beta prior for pi_level
  g_hyp = NULL,               # Optional Gamma(a,b) on g_level per detail level
  a_sig = 2.5, b_sig = 0.02   # IG prior on sigma2
) {
  stopifnot(length(Y_list) > 0)
  M <- ncol(Y_list[[1]])
  N <- length(Y_list)

  # Stack wavelet coefficients across samples i and channels m
  st  <- wt_stack_channel(Y_list, wf=wf, J=J, boundary=boundary)
  D   <- st$D_arr                       # [ncoeff x N x M]
  maps<- st$maps                        # per-channel maps (length M)
  ncoeff <- dim(D)[1]

  lev_names <- names(maps[[1]]$map$idx)
  det_names <- lev_names[grepl("^d", lev_names)]
  s_name    <- lev_names[grepl("^s", lev_names)]

  # Initialize missing wpar slots
  if (is.null(wpar$pi_level))  wpar$pi_level <- setNames(rep(0.5, length(det_names)), det_names)
  if (is.null(wpar$g_level))   wpar$g_level  <- setNames(rep(2.0,  length(det_names)), det_names)
  if (is.null(wpar$gamma_ch))  wpar$gamma_ch <- lapply(1:M, function(m) rbinom(ncoeff, 1, 0.2))

  # ---------- 1) Update gamma (detail levels) with spike vs slab evidence -----
  for (m in 1:M) {
    Dm  <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)  # [ncoeff x N]
    gam <- wpar$gamma_ch[[m]]

    for (lev in det_names) {
      ids <- maps[[m]]$map$idx[[lev]]
      if (length(ids) == 0) next
      pi_j <- wpar$pi_level[[lev]]
      g_j  <- wpar$g_level[[lev]]

      # Marginal under spike: D ~ N(0, sigma2)
      # Marginal under slab (beta integrated out): D ~ N(0, (1+g_j) * sigma2)
      v_spike <- sigma2
      v_slab  <- (1 + g_j) * sigma2
      Dsub    <- Dm[ids, , drop=FALSE]   # [n_ids x N]

      ll_spike <- -0.5 * rowSums( log(2*pi*v_spike) + (Dsub^2)/v_spike )
      ll_slab  <- -0.5 * rowSums( log(2*pi*v_slab ) + (Dsub^2)/v_slab  )

      logit_val <- log(pi_j) + ll_slab - (log(1 - pi_j) + ll_spike)
      p1        <- plogis(pmax(pmin(logit_val, 35), -35))
      gam[ids]  <- rbinom(length(ids), 1, p1)
    }

    # Force scaling band on (standard wavelet practice)
    if (length(s_name) == 1) {
      ids_s <- maps[[m]]$map$idx[[s_name]]
      if (length(ids_s) > 0) gam[ids_s] <- 1L
    }
    wpar$gamma_ch[[m]] <- gam
  }

  # ---------- 2) Update g_level (detail-specific slab inflation) --------------
  for (lev in det_names) {
    # Conjugate Gamma prior on g_j; use simple moment update via marginal score
    # We adopt the same IG-moment style as your previous block but with common sigma2.
    shape0 <- if (!is.null(g_hyp)) g_hyp[lev, "shape"] else 2.0
    rate0  <- if (!is.null(g_hyp)) g_hyp[lev, "rate"]  else 2.0

    # Use slab responsibility: only active rows contribute with 'slab precision'
    ss_over_sigma <- 0
    n_sel_total   <- 0
    for (m in 1:M) {
      ids <- maps[[m]]$map$idx[[lev]]
      if (length(ids) == 0) next
      Dm  <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
      sel <- wpar$gamma_ch[[m]][ids] == 1
      if (any(sel)) {
        Did <- Dm[ids[sel], , drop=FALSE]
        ss_over_sigma <- ss_over_sigma + sum(Did^2) / sigma2
        n_sel_total   <- n_sel_total + nrow(Did)
      }
    }
    shape_post <- shape0 + 0.5 * n_sel_total
    rate_post  <- rate0  + 0.5 * ss_over_sigma
    # Sample g_j (keeps your stochastic flavour); for deterministic EB, use rate/shape
    wpar$g_level[[lev]] <- invgamma::rinvgamma(1, shape=shape_post, rate=rate_post)
  }

  # ---------- 3) Update pi_level (detail inclusion rates) ---------------------
  for (lev in det_names) {
    n1 <- 0; n0 <- 0
    for (m in 1:M) {
      ids <- maps[[m]]$map$idx[[lev]]
      if (length(ids) == 0) next
      gm <- wpar$gamma_ch[[m]][ids]
      n1 <- n1 + sum(gm == 1)
      n0 <- n0 + sum(gm == 0)
    }
    wpar$pi_level[[lev]] <- rbeta(1, a_pi + n1, b_pi + n0)
  }

  # ---------- 4) Update beta_ch (posterior mean for active coefficients) ------
  # beta_{r,m} | D, sigma2, gamma=1 ~ N( (n*mean)/(n + 1/g_j), sigma2/(n + 1/g_j) )
  beta_ch <- lapply(1:M, function(m) numeric(ncoeff))
  for (m in 1:M) {
    Dm  <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
    gam <- wpar$gamma_ch[[m]]
    b   <- numeric(ncoeff)

    for (lev in det_names) {
      ids <- maps[[m]]$map$idx[[lev]]
      if (length(ids) == 0) next

      g_j <- wpar$g_level[[lev]]
      n   <- N
      Dbar<- rowMeans(Dm[ids, , drop=FALSE])              # [|ids|]

      # posterior mean for active; zero for inactive
      shrink <- (n) / (n + 1 / g_j)
      b_active <- shrink * Dbar

      is_on <- (gam[ids] == 1L)
      if (any(is_on)) b[ids[is_on]] <- b_active[is_on]
      # off remain 0
    }

    # Keep scaling band unshrunk (mean equals average of scaling coeffs)
    if (length(s_name) == 1) {
      ids_s <- maps[[m]]$map$idx[[s_name]]
      if (length(ids_s) > 0) b[ids_s] <- rowMeans(Dm[ids_s, , drop=FALSE])
    }
    beta_ch[[m]] <- b
  }

  # ---------- 5) Update sigma2 (cluster noise, common across all coeffs) ------
  # D_{r,i,m} - beta_{r,m} ~ N(0, sigma2)
  ss <- 0
  for (m in 1:M) {
    Dm <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
    resid <- sweep(Dm, 1, beta_ch[[m]], FUN = "-")
    ss <- ss + sum(resid^2)
  }
  n_eff <- ncoeff * N * M
  sigma2 <- invgamma::rinvgamma(1, shape = a_sig + 0.5 * n_eff, rate = b_sig + 0.5 * ss)

  list(wpar = wpar, beta_ch = beta_ch, sigma2 = sigma2, maps = maps)
}

# ------------------------ Likelihood & MH updates ----------------------------

# ------------------------ Likelihood with tau_B -------------------------------

loglik_icm <- function(y_norm, t, L, eta, kern_cfg, kp, tau_B = 1.0, mu = NULL) {
  # Build B with amplitude tau_B (unit-trace shape for stability)
  M <- length(eta)
  Bshape <- tcrossprod(L)
  trB <- sum(diag(Bshape))
  if (trB > 0) Bshape <- Bshape * (M / trB)
  B  <- tau_B * Bshape

  Kx <- kern_cfg$fun(t, kp)
  K  <- kronecker_icm(B, Kx, eta)

  Yc <- if (is.null(mu)) y_norm else (y_norm - mu)
  yv <- as.numeric(Yc)
  chol_loglik(yv, K, jitter=1e-6)$ll
}

# ------------------------ MH update for tau_B --------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

mh_update_tauB <- function(tau_B, step, t, Y_list, L, eta, kern_cfg, kp, acc) {
  cur <- tau_B
  prp <- rlnorm(1, log(cur), step)  # log-normal RW

  ll_cur <- sum(sapply(Y_list, loglik_icm, t=t, L=L, eta=eta,
                       kern_cfg=kern_cfg, kp=kp, tau_B=cur))
  ll_prp <- sum(sapply(Y_list, loglik_icm, t=t, L=L, eta=eta,
                       kern_cfg=kern_cfg, kp=kp, tau_B=prp))

  # Robust prior: half-Cauchy on sqrt(tau_B)  =>  p(tau_B) ∝ 1/(1+tau_B)
  lp_cur <- -log1p(cur)
  lp_prp <- -log1p(prp)

  q_cgpr <- dlnorm(cur, log(prp), step, log=TRUE)
  q_prgc <- dlnorm(prp, log(cur), step, log=TRUE)

  a <- (ll_prp + lp_prp + q_cgpr) - (ll_cur + lp_cur + q_prgc)
  if (is.finite(a) && log(runif(1)) < a) { tau_B <- prp; acc$tauB <- (acc$tauB %||% 0L) + 1L }
  list(tau_B=tau_B, acc=acc)
}

# ------------------------ MH updates (now carry tau_B) -----------------------

mh_update_kernel <- function(kern_cfg, kp, t, Ynorm_list, L, eta, tau_B, acc) {
  for(pn in kern_cfg$pnames) {
    cur <- kp[[pn]]
    if (pn == "period") {
      z_cur <- logit(pmin(pmax(cur,1e-6), 1-1e-6))
      z_prp <- rnorm(1, z_cur, kern_cfg$prop_sd[[pn]])
      prp   <- invlogit(z_prp)
      kp_prop <- kp; kp_prop[[pn]] <- prp
      ll_cur <- sum(sapply(Ynorm_list, loglik_icm, t=t, L=L, eta=eta, kern_cfg=kern_cfg, kp=kp,       tau_B=tau_B))
      ll_prp <- sum(sapply(Ynorm_list, loglik_icm, t=t, L=L, eta=eta, kern_cfg=kern_cfg, kp=kp_prop, tau_B=tau_B))
      lp_cur <- kern_cfg$prior(kp)
      lp_prp <- kern_cfg$prior(kp_prop)
      a <- (ll_prp + lp_prp) - (ll_cur + lp_cur)
      if (is.finite(a) && log(runif(1)) < a) { kp <- kp_prop; acc$kern <- acc$kern + 1 }
    } else {
      prp   <- rlnorm(1, log(cur), kern_cfg$prop_sd[[pn]])
      kp_prop <- kp; kp_prop[[pn]] <- prp
      ll_cur <- sum(sapply(Ynorm_list, loglik_icm, t=t, L=L, eta=eta, kern_cfg=kern_cfg, kp=kp,       tau_B=tau_B))
      ll_prp <- sum(sapply(Ynorm_list, loglik_icm, t=t, L=L, eta=eta, kern_cfg=kern_cfg, kp=kp_prop, tau_B=tau_B))
      lp_cur <- kern_cfg$prior(kp)
      lp_prp <- kern_cfg$prior(kp_prop)
      q_cgpr <- dlnorm(cur, log(prp), kern_cfg$prop_sd[[pn]], log=TRUE)
      q_prgc <- dlnorm(prp, log(cur), kern_cfg$prop_sd[[pn]], log=TRUE)
      a <- (ll_prp + lp_prp + q_cgpr) - (ll_cur + lp_cur + q_prgc)
      if (is.finite(a) && log(runif(1)) < a) { kp <- kp_prop; acc$kern <- acc$kern + 1 }
    }
  }
  list(kp=kp, acc=acc)
}

mh_update_L <- function(L, step, t, Ynorm_list, eta, kern_cfg, kp, tau_B, acc) {
  th  <- pack_L(L)
  thp <- th + rnorm(length(th), 0, step)
  Lp  <- unpack_L(thp, nrow(L))
  ll_cur <- sum(sapply(Ynorm_list, loglik_icm, t=t, L=L,  eta=eta, kern_cfg=kern_cfg, kp=kp, tau_B=tau_B))
  ll_prp <- sum(sapply(Ynorm_list, loglik_icm, t=t, L=Lp, eta=eta, kern_cfg=kern_cfg, kp=kp, tau_B=tau_B))
  lp_cur <- sum(dnorm(th,  0, 1, log=TRUE))
  lp_prp <- sum(dnorm(thp, 0, 1, log=TRUE))
  a <- (ll_prp + lp_prp) - (ll_cur + lp_cur)
  if (is.finite(a) && log(runif(1)) < a) { L <- Lp; acc$L <- acc$L + 1 }
  list(L=L, acc=acc)
}

mh_update_eta <- function(eta, step, t, Ynorm_list, L, kern_cfg, kp, tau_B, acc) {
  for(j in seq_along(eta)) {
    cur <- eta[j]
    prp <- rlnorm(1, log(cur), step)
    etap <- eta; etap[j] <- prp
    ll_cur <- sum(sapply(Ynorm_list, loglik_icm, t=t, L=L, eta=eta,  kern_cfg=kern_cfg, kp=kp, tau_B=tau_B))
    ll_prp <- sum(sapply(Ynorm_list, loglik_icm, t=t, L=L, eta=etap, kern_cfg=kern_cfg, kp=kp, tau_B=tau_B))
    lp_cur <- invgamma::dinvgamma(cur, 2, 0.1, log=TRUE)
    lp_prp <- invgamma::dinvgamma(prp, 2, 0.1, log=TRUE)
    q_cgpr <- dlnorm(cur, log(prp), step, log=TRUE)
    q_prgc <- dlnorm(prp, log(cur), step, log=TRUE)
    a <- (ll_prp + lp_prp + q_cgpr) - (ll_cur + lp_cur + q_prgc)
    if (is.finite(a) && log(runif(1)) < a) { eta <- etap; acc$eta <- acc$eta + 1 }
  }
  list(eta=eta, acc=acc)
}

# ------------------------ Carlin–Chib kernel switch --------------------------

cc_switch_kernel <- function(kern_list, cur_idx, thetas, t, Ynorm_list, L, eta, tau_B) {
  Mmod <- length(kern_list)
  p_m  <- rep(1/Mmod, Mmod)
  theta_draws <- vector("list", Mmod)
  for (m in 1:Mmod) theta_draws[[m]] <- if (m == cur_idx) thetas[[m]] else kern_list[[m]]$pstar()

  logw <- rep(NA_real_, Mmod)
  for (m in 1:Mmod) {
    kc <- kern_list[[m]]
    kp <- theta_draws[[m]]

    # NOTE: now evaluates with tau_B
    ll_m <- sum(sapply(
      Ynorm_list,
      loglik_icm,
      t = t, L = L, eta = eta, kern_cfg = kc, kp = kp, tau_B = tau_B
    ))

    # pseudo-priors for theaux kernels (unchanged)
    rest <- sum(sapply(setdiff(1:Mmod, m), function(j) {
      th <- theta_draws[[j]]
      if (kern_list[[j]]$name == "Periodic") {
        dgamma(th$l_scale, 3, 2, log = TRUE) + dbeta(th$period, 5, 5, log = TRUE)
      } else {
        dgamma(th$l_scale, 2, 2, log = TRUE)
      }
    }))

    logw[m] <- log(p_m[m]) + ll_m + rest
  }
  w <- exp(logw - max(logw)); w <- w / sum(w)
  new_idx <- sample.int(Mmod, 1, prob = w)
  list(idx = new_idx, thetas = theta_draws, weights = w)
}

# ------------------------ Walker slice sampler (retrospective) ---------------

stick_to_pi <- function(v) {
  K <- length(v)
  pi <- numeric(K); cum <- 1
  for(k in 1:K) { pi[k] <- v[k] * cum; cum <- cum * (1 - v[k]) }
  pi
}

extend_sticks_until <- function(v, alpha, threshold) {
  tail <- prod(1 - v)
  while (tail > threshold) {
    v_new <- rbeta(1, 1, alpha)
    v <- c(v, v_new)
    tail <- tail * (1 - v_new)
  }
  v
}

update_v_given_z <- function(v, z, alpha) {
  K <- length(v)
  n_k <- tabulate(z, nbins = K)
  n_tail <- rev(cumsum(rev(n_k)))
  for(k in 1:K) {
    a <- 1 + n_k[k]
    b <- alpha + if(k<K) n_tail[k+1] else 0
    v[k] <- rbeta(1, a, b)
  }
  v
}

# ------------------------ Initialization (adds tau_B) ------------------------

draw_new_cluster_params <- function(M, P, t, kern_list, wf="la8", J=NULL, boundary="periodic",
                                    use_besov_pi=FALSE, b=0.5, c2=0.8,
                                    use_besov_g =TRUE, a=0.5) {
  if (is.null(J)) J <- floor(log2(P))
  zeros <- matrix(0, nrow=P, ncol=M)
  tmp   <- wt_forward_mat(zeros, wf=wf, J=J, boundary=boundary)
  lev_names <- names(tmp[[1]]$map$idx)
  det_names <- lev_names[grepl("^d", lev_names)]
  ncoeff    <- length(tmp[[1]]$coeff)
  beta_ch <- lapply(1:M, function(m) numeric(ncoeff))  # <--- NEW: cluster mean wavelet coeffs

  pi_level <- if (use_besov_pi) besov_pi_schedule(J, b=b, c2=c2) else setNames(rep(0.5, length(det_names)), det_names)
  g_hyp    <- if (use_besov_g)  besov_g_hyper(J, a=a) else NULL
  g_level  <- setNames(rep(2.0, length(det_names)), det_names)
  gamma_ch <- lapply(1:M, function(m) rbinom(ncoeff, 1, 0.2))
  thetas   <- lapply(kern_list, function(kc) kc$pstar())

  list(
    wpar = list(lev_names=lev_names, pi_level=pi_level, g_level=g_level, gamma_ch=gamma_ch),
    g_hyp = g_hyp,
    kern_idx = sample.int(length(kern_list),1),
    thetas   = thetas,
    L = diag(M),
    eta = rep(0.05, M),
    tau_B = 1.0,                     # <--- NEW: per-cluster amplitude
    acc = list(kern=0, L=0, eta=0, tauB=0),
    sigma2_i = numeric(0),    beta_ch = beta_ch,           # <--- NEW
    sigma2  = 0.05              # <--- NEW: cluster-level noise in wavelet space
  )
}

# ------------ Hardened inverse: sanitize and validate before idwt ------------

wt_inverse_1d <- function(coeff_vec, map) {
  # 1) Rebuild the dwt list with expected names
  J <- map$J
  w <- vector("list", J + 1L)
  names(w) <- c(paste0("d", 1:J), paste0("s", J))
  for (lev in 1:J) {
    ids <- map$idx[[paste0("d", lev)]]
    stopifnot(!is.null(ids), length(ids) > 0)
    w[[paste0("d", lev)]] <- as.numeric(coeff_vec[ids])
  }
  ids_s <- map$idx[[paste0("s", J)]]
  stopifnot(!is.null(ids_s), length(ids_s) > 0)
  w[[paste0("s", J)]] <- as.numeric(coeff_vec[ids_s])

  # 2) Sanitize: replace non-finite with 0 (safe default for reconstruction)
  for (nm in names(w)) {
    bad <- !is.finite(w[[nm]])
    if (any(bad)) w[[nm]][bad] <- 0
  }

  # 3) Attach attributes waveslim expects
  attr(w, "wavelet")  <- map$wf
  attr(w, "boundary") <- map$boundary
  class(w) <- "dwt"

  # 4) Final structural checks
  stopifnot(
    is.character(attr(w, "wavelet")),
    is.character(attr(w, "boundary")),
    all(vapply(w, is.numeric, logical(1L)))
  )

  # 5) Reconstruct
  waveslim::idwt(w)
}

# --------------- Hardened forward: store explicit metadata -------------------

wt_forward_mat <- function(y_mat, wf = "la8", J = NULL, boundary = "periodic") {
  # Apply wavelet transform to each channel of a matrix
  M <- ncol(y_mat)
  out <- vector("list", M)
  for (m in 1:M) {
    out[[m]] <- wt_forward_1d(y_mat[, m], wf = wf, J = J, boundary = boundary)
  }
  out
}

wt_forward_1d <- function(y, wf = "la8", J = NULL, boundary = "periodic") {
  P <- length(y)
  J <- ensure_dyadic_J(P, J)

  w <- waveslim::dwt(y, wf = wf, n.levels = J, boundary = boundary)

  # Flatten in canonical order: d1, d2, ..., dJ, sJ
  vec <- c(w$d1)
  idx <- list(d1 = seq_along(w$d1))
  off <- length(w$d1)
  if (J >= 2) {
    for (lev in 2:J) {
      nm <- paste0("d", lev)
      v  <- w[[nm]]
      vec <- c(vec, v)
      idx[[nm]] <- (off + 1):(off + length(v))
      off <- off + length(v)
    }
  }
  s_nm <- paste0("s", J)
  vec <- c(vec, w[[s_nm]])
  idx[[s_nm]] <- (off + 1):(off + length(w[[s_nm]]))

  # Keep full metadata so inverse can reproduce the dwt object precisely
  list(
    coeff = as.numeric(vec),
    map   = list(J = J, wf = wf, boundary = boundary, P = P, idx = idx)
  )
}

#

# ------- Extra guard where you reconstruct curves from per-curve transforms ---

wt_inverse_mat_soft <- function(y_mat, wtf_list, wpar, shrink_scaling = FALSE) {
  stopifnot(length(wtf_list) == ncol(y_mat))
  P <- nrow(y_mat); M <- ncol(y_mat)
  y_norm <- matrix(0, P, M)
  `%OR%` <- function(a, b) if (is.null(a)) b else a

  for (m in 1:M) {
    cvec <- as.numeric(wtf_list[[m]]$coeff)
    map  <- wtf_list[[m]]$map
    gam  <- wpar$gamma_ch[[m]]

    # Realign if needed
    if (length(gam) != length(cvec)) {
      warning(sprintf("gamma length (%d) != coeff length (%d); truncating",
                      length(gam), length(cvec)))
      L <- min(length(gam), length(cvec))
      gam  <- gam[seq_len(L)]
      cvec <- cvec[seq_len(L)]
    }

    lev_names <- names(map$idx)
    det_names <- lev_names[grepl("^d", lev_names)]
    s_name    <- lev_names[grepl("^s", lev_names)]

    for (lev in det_names) {
      ids <- map$idx[[lev]]
      if (!length(ids)) next
      g_j   <- (wpar$g_level[[lev]] %OR% 1.0)
      if (!is.finite(g_j) || g_j < 0) g_j <- 1.0
      kappa <- g_j / (1 + g_j)
      is_on <- (gam[ids] == 1L)
      if (any(!is_on)) cvec[ids[!is_on]] <- 0.0
      if (any(is_on))  cvec[ids[ is_on]] <- kappa * cvec[ids[ is_on]]
    }

    if (length(s_name) == 1 && !is.na(s_name) && shrink_scaling) {
      ids_s <- map$idx[[s_name]]
      if (length(ids_s) > 0) {
        g_s <- (wpar$g_level[[s_name]] %OR% Inf)
        if (is.finite(g_s) && g_s >= 0) {
          kappa_s <- g_s / (1 + g_s)
          cvec[ids_s] <- kappa_s * cvec[ids_s]
        }
      }
    }

    # Sanitize
    bad <- !is.finite(cvec)
    if (any(bad)) cvec[bad] <- 0

    y_norm[, m] <- wt_inverse_1d(cvec, map)
  }
  y_norm
}

# Build a P x M mean matrix mu from cluster beta_ch
compute_mu_from_beta <- function(beta_ch, wf, J, boundary, P) {
  M <- length(beta_ch)
  zeros <- matrix(0, nrow=P, ncol=M)
  tmpl  <- wt_forward_mat(zeros, wf=wf, J=J, boundary=boundary)
  mu <- matrix(0, P, M)
  for (m in 1:M) {
    mu[,m] <- wt_inverse_1d(beta_ch[[m]], tmpl[[m]]$map)
  }
  mu
}

# Convenience wrapper: reconstruct a list of series using soft shrinkage
reconstruct_normals_for_cluster_soft <- function(Y_list, wf, J, boundary, wpar,
                                                 shrink_scaling = FALSE) {
  N <- length(Y_list)
  out <- vector("list", N)
  for (i in 1:N) {
    wtf <- wt_forward_mat(Y_list[[i]], wf = wf, J = J, boundary = boundary)
    out[[i]] <- wt_inverse_mat_soft(Y_list[[i]], wtf, wpar, shrink_scaling = shrink_scaling)
  }
  out
}

# ------------------------ The full driver (tau_B-enabled) --------------------

run_model <- function(
    Y, t,
    n_iter=6000, burn=3000, thin=5,
    alpha_prior = c(10,1),
    wf="la8", J=NULL, boundary="periodic",
    mh_step_L=0.03, mh_step_eta=0.10, mh_step_tauB=0.15,   # <--- added tauB step
    use_besov_pi=TRUE, use_besov_g=TRUE,
    revealed_idx = integer(0),
    emp_bayes_init_iter = 80,
    unpin_after_warmstart = FALSE,
    K_init = 5,
    besov_c2 = 0.5
) {
  N <- length(Y); P <- nrow(Y[[1]]); M <- ncol(Y[[1]]); J <- ensure_dyadic_J(P, J)

  alpha <- rgamma(1, shape=alpha_prior[1], rate=alpha_prior[2])

  v <- rbeta(K_init, 1, alpha)
  pi <- stick_to_pi(v)
  K  <- length(v)

  params <- vector("list", K)
  for (k in 1:K) {
    params[[k]] <- draw_new_cluster_params(M, P, t, kernels, wf, J, boundary,
                                           use_besov_pi=use_besov_pi,
                                           use_besov_g=use_besov_g)
  }

  z <- sample.int(K, N, replace=TRUE)
  if (length(revealed_idx)) z[revealed_idx] <- 1L

  # ---- warm-start on cluster 1 (revealed normals only)
  if (length(revealed_idx) > 0) {
    idx <- sort(unique(revealed_idx))
    Yk  <- Y[idx]
    if (length(params[[1]]$sigma2_i) != length(idx) || any(!is.finite(params[[1]]$sigma2_i))) {
      params[[1]]$sigma2_i <- rep(0.05, length(idx))
    }
    
    cat(sprintf("🔄 Starting Empirical Bayes warm-start (%d iterations) on %d revealed normals\n",
                emp_bayes_init_iter, length(idx))); flush.console()
    
    for (it in 1:emp_bayes_init_iter) {
  # Cluster-level Ray–Mallick updates: (wpar, beta_ch, sigma2)
  upd <- update_cluster_wavelet_params_rm(
    Y_list = Yk, wf = wf, J = J, boundary = boundary,
    wpar   = params[[1]]$wpar,
    sigma2 = params[[1]]$sigma2,     # cluster σ² (scalar)
    a_pi = 1, b_pi = 1,
    g_hyp = params[[1]]$g_hyp,
    a_sig = 2.5, b_sig = 0.02
  )
  params[[1]]$wpar    <- upd$wpar
  params[[1]]$beta_ch <- upd$beta_ch      # store cluster β (per channel, wavelet domain)
  params[[1]]$sigma2  <- upd$sigma2       # store cluster σ²

  # Build μ from β (not from per-curve σ_i²)
  mu_k <- compute_mu_from_beta(
    beta_ch  = params[[1]]$beta_ch,
    wf       = wf, J = J, boundary = boundary,
    P        = P
  )

  # Residuals for GP updates
  Yk_centered <- lapply(Yk, function(y) y - mu_k)

  # (kernel/η/L/τ_B updates unchanged)
  cc <- cc_switch_kernel(
    kernels,
    params[[1]]$kern_idx, params[[1]]$thetas,
    t, Yk_centered, params[[1]]$L, params[[1]]$eta,
    params[[1]]$tau_B
  )
  params[[1]]$kern_idx <- cc$idx
  params[[1]]$thetas   <- cc$thetas

  kc <- kernels[[ params[[1]]$kern_idx ]]
  kp <- params[[1]]$thetas[[ params[[1]]$kern_idx ]]

  tmp <- mh_update_kernel(kc, kp, t, Yk_centered, params[[1]]$L, params[[1]]$eta, params[[1]]$tau_B, params[[1]]$acc)
  params[[1]]$thetas[[ params[[1]]$kern_idx ]] <- tmp$kp; params[[1]]$acc <- tmp$acc

  tmp <- mh_update_L(params[[1]]$L, mh_step_L, t, Yk_centered, params[[1]]$eta, kc,
                     params[[1]]$thetas[[ params[[1]]$kern_idx ]], params[[1]]$tau_B, params[[1]]$acc)
  params[[1]]$L   <- tmp$L; params[[1]]$acc <- tmp$acc

  tmp <- mh_update_eta(params[[1]]$eta, mh_step_eta, t, Yk_centered, params[[1]]$L, kc,
                       params[[1]]$thetas[[ params[[1]]$kern_idx ]], params[[1]]$tau_B, params[[1]]$acc)
  params[[1]]$eta <- tmp$eta; params[[1]]$acc <- tmp$acc

  tmp <- mh_update_tauB(params[[1]]$tau_B, mh_step_tauB, t, Yk_centered,
                        params[[1]]$L, params[[1]]$eta, kc,
                        params[[1]]$thetas[[ params[[1]]$kern_idx ]], params[[1]]$acc)
  params[[1]]$tau_B <- tmp$tau_B; params[[1]]$acc <- tmp$acc

  if (it %% 100 == 0) {
    progress_pct <- round(100 * it / emp_bayes_init_iter, 1)
    cat(sprintf("\r📈 Warm-start iteration %d/%d (%.1f%%)",
                it, emp_bayes_init_iter, progress_pct)); flush.console()
  }
}

    
    cat("\n✅ Warm-start completed! Starting main MCMC...\n"); flush.console()
  }

  keep <- floor((n_iter - burn)/thin)
  Z_s <- matrix(NA_integer_, keep, N); K_s <- integer(keep); alpha_s <- numeric(keep); kern_s <- integer(keep)
  acc_hist <- data.frame(iter=integer(0), accL=double(0), accEta=double(0), accKer=double(0), accTauB=double(0))
  
  # NEW: Track kernel parameters for the largest cluster
  kernel_params_trace <- list()
  
  # NEW: Track wavelet coefficients for the largest cluster
  wavelet_coeff_trace <- list()
  
  # Track K values for running average
  K_history <- integer(0)

  cat(sprintf("🔄 Starting MCMC run (N=%d, P=%d, M=%d, iters=%d, burn=%d, thin=%d, revealed=%d, K_init=%d)\n",
              N,P,M,n_iter,burn,thin,length(revealed_idx), K)); flush.console()

  pin_revealed <- length(revealed_idx) > 0 && !unpin_after_warmstart
  sidx <- 0

  for(iter in 1:n_iter) {
    # Walker slice expansion
    pi <- stick_to_pi(v)
    u  <- sapply(1:N, function(i) runif(1, 0, pi[z[i]]))
    u_star <- min(u)
    v <- extend_sticks_until(v, alpha, u_star)
    pi <- stick_to_pi(v)
    K  <- length(v)
    while(length(params) < K) params[[length(params)+1]] <-
      draw_new_cluster_params(M,P,t,kernels,wf,J,boundary,use_besov_pi=use_besov_pi,use_besov_g=use_besov_g)

    # Assignments
    for(i in 1:N) {
      if (pin_revealed && (i %in% revealed_idx)) { z[i] <- 1L; next }
      S <- which(pi > u[i]); if(length(S)==0) S <- 1
      logw <- rep(-Inf, length(S))
      for(ss in seq_along(S)) {
        k <- S[ss]
        # Compute or reuse cluster mean μ_k built from cluster-level shrinkage
       # Compute or reuse μ_k from current β_ch (cluster-level)
if (is.null(params[[k]]$mu_cached) || is.null(params[[k]]$mu_cached_iter) || params[[k]]$mu_cached_iter != iter) {
  if (!is.null(params[[k]]$beta_ch) && length(params[[k]]$beta_ch)) {
    mu_k <- compute_mu_from_beta(
      beta_ch = params[[k]]$beta_ch,
      wf = wf, J = J, boundary = boundary,
      P = P
    )
  } else {
    mu_k <- matrix(0, nrow(Y[[1]]), ncol(Y[[1]]))  # safe fallback
  }
  params[[k]]$mu_cached <- mu_k
  params[[k]]$mu_cached_iter <- iter
} else {
  mu_k <- params[[k]]$mu_cached
}


        # Residual for candidate assignment = raw curve minus shrunk cluster mean
        y_resid_i <- Y[[i]] - mu_k

        kc <- kernels[[ params[[k]]$kern_idx ]]
        kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]
        lp <- log(pi[k])
        ll <- loglik_icm(y_resid_i, t, params[[k]]$L, params[[k]]$eta, kc, kp,
                         tau_B = params[[k]]$tau_B, mu = NULL)  # already centered
        logw[ss] <- lp + ll

      }
      w <- exp(logw - max(logw)); w <- w / sum(w)
      z[i] <- S[sample.int(length(S), 1, prob=w)]
    }

    # Update sticks
    v <- update_v_given_z(v, z, alpha)
    pi <- stick_to_pi(v)
    K  <- length(v)

    # Per-cluster updates
    for(k in 1:K) {
     idx <- which(z == k)
if (length(idx) == 0) next
Yk <- Y[idx]

# Cluster-level Ray–Mallick update (wpar, beta_ch, sigma2)
upd <- update_cluster_wavelet_params_rm(
  Y_list = Yk, wf = wf, J = J, boundary = boundary,
  wpar   = params[[k]]$wpar,
  sigma2 = params[[k]]$sigma2,
  a_pi = 1, b_pi = 1,
  g_hyp = params[[k]]$g_hyp,
  a_sig = 2.5, b_sig = 0.02
)
params[[k]]$wpar    <- upd$wpar
params[[k]]$beta_ch <- upd$beta_ch
params[[k]]$sigma2  <- upd$sigma2

# Build μ from β (cluster mean in the data space)
mu_k <- compute_mu_from_beta(
  beta_ch = params[[k]]$beta_ch,
  wf = wf, J = J, boundary = boundary,
  P = P
)

# Residuals for GP updates
Yk_centered <- lapply(Yk, function(y) y - mu_k)

# Cache μ for diagnostics/plots
params[[k]]$mu_cached <- mu_k
params[[k]]$mu_cached_iter <- iter

# Kernel selection + parameter MH updates (unchanged)
cc <- cc_switch_kernel(
  kernels,
  params[[k]]$kern_idx, params[[k]]$thetas,
  t, Yk_centered, params[[k]]$L, params[[k]]$eta,
  params[[k]]$tau_B
)
params[[k]]$kern_idx <- cc$idx
params[[k]]$thetas   <- cc$thetas

kc <- kernels[[ params[[k]]$kern_idx ]]
kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]

tmp <- mh_update_kernel(kc, kp, t, Yk_centered, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$acc)
params[[k]]$thetas[[ params[[k]]$kern_idx ]] <- tmp$kp; params[[k]]$acc <- tmp$acc

tmp <- mh_update_L(params[[k]]$L, mh_step_L, t, Yk_centered, params[[k]]$eta, kc,
                   params[[k]]$thetas[[ params[[k]]$kern_idx ]], params[[k]]$tau_B, params[[k]]$acc)
params[[k]]$L   <- tmp$L; params[[k]]$acc <- tmp$acc

tmp <- mh_update_eta(params[[k]]$eta, mh_step_eta, t, Yk_centered, params[[k]]$L, kc,
                     params[[k]]$thetas[[ params[[k]]$kern_idx ]], params[[k]]$tau_B, params[[k]]$acc)
params[[k]]$eta <- tmp$eta; params[[k]]$acc <- tmp$acc

tmp <- mh_update_tauB(params[[k]]$tau_B, mh_step_tauB, t, Yk_centered,
                      params[[k]]$L, params[[k]]$eta, kc,
                      params[[k]]$thetas[[ params[[k]]$kern_idx ]], params[[k]]$acc)
params[[k]]$tau_B <- tmp$tau_B; params[[k]]$acc <- tmp$acc


      cc <- cc_switch_kernel(
        kernels,
        params[[k]]$kern_idx, params[[k]]$thetas,
        t, Yk_centered, params[[k]]$L, params[[k]]$eta,
        params[[k]]$tau_B
      )
      params[[k]]$kern_idx <- cc$idx
      params[[k]]$thetas   <- cc$thetas

      kc <- kernels[[ params[[k]]$kern_idx ]]
      kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]

      tmp <- mh_update_kernel(kc, kp, t, Yk_centered, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$acc)
      params[[k]]$thetas[[ params[[k]]$kern_idx ]] <- tmp$kp; params[[k]]$acc <- tmp$acc

      tmp <- mh_update_L(params[[k]]$L, mh_step_L, t, Yk_centered, params[[k]]$eta, kc,
                   params[[k]]$thetas[[ params[[k]]$kern_idx ]], params[[k]]$tau_B, params[[k]]$acc)

      params[[k]]$L   <- tmp$L; params[[k]]$acc <- tmp$acc

      tmp <- mh_update_eta(params[[k]]$eta, mh_step_eta, t, Yk_centered, params[[k]]$L, kc,
                     params[[k]]$thetas[[ params[[k]]$kern_idx ]], params[[k]]$tau_B, params[[k]]$acc)

      params[[k]]$eta <- tmp$eta; params[[k]]$acc <- tmp$acc

      # NEW: MH for tau_B
      tmp <- mh_update_tauB(params[[k]]$tau_B, mh_step_tauB, t, Yk_centered,
                      params[[k]]$L, params[[k]]$eta, kc,
                      params[[k]]$thetas[[ params[[k]]$kern_idx ]], params[[k]]$acc)

      params[[k]]$tau_B <- tmp$tau_B; params[[k]]$acc <- tmp$acc

    }

    # Update alpha (Escobar–West)
    Kocc <- length(unique(z))
    eta_aux <- rbeta(1, alpha + 1, N)
    mix <- (alpha_prior[1] + Kocc - 1) /
      (N*(alpha_prior[2] - log(eta_aux)) + alpha_prior[1] + Kocc - 1)
    if(runif(1) < mix) {
      alpha <- rgamma(1, alpha_prior[1] + Kocc, alpha_prior[2] - log(eta_aux))
    } else {
      alpha <- rgamma(1, alpha_prior[1] + Kocc - 1, alpha_prior[2] - log(eta_aux))
    }

    # Progress / bookkeeping
    if(iter %% 10 == 0) {
      # Track current K value
      current_K <- length(unique(z))
      K_history <- c(K_history, current_K)
      
      occ <- sort(unique(z))
      A <- sapply(occ, function(k) c(params[[k]]$acc$L, params[[k]]$acc$eta, params[[k]]$acc$kern, params[[k]]$acc$tauB))
      if(is.matrix(A)) A <- rowMeans(A) else A <- rep(NA,4)
      acc_hist <- rbind(acc_hist, data.frame(iter=iter, accL=A[1], accEta=A[2], accKer=A[3], accTauB=A[4]))
      for(k in occ) params[[k]]$acc <- list(kern=0, L=0, eta=0, tauB=0)

      # Calculate average K
      avg_K <- round(mean(K_history), 1)
      
      progress_pct <- round(100 * iter / n_iter, 1)
      # Calculate expected number of clusters: alpha * log(1 + N/alpha)
      expected_K <- round(alpha * log(1 + N/alpha), 1)
      # Use \r to overwrite the same line instead of creating new lines
      cat(sprintf("\r📈 Iteration %d/%d (%.1f%%), K=%d (avg=%.1f, exp=%.1f), alpha=%.3f", 
                  iter, n_iter, progress_pct, current_K, avg_K, expected_K, alpha))
      flush.console()
    }

    if(iter > burn && ((iter - burn) %% thin == 0)) {
      sidx <- sidx + 1
      Z_s[sidx,] <- z
      K_s[sidx]  <- length(unique(z))
      alpha_s[sidx] <- alpha
      tab <- table(z); k_big <- as.integer(names(which.max(tab)))
      kern_s[sidx] <- params[[k_big]]$kern_idx
      
      # NEW: Store kernel parameters for the largest cluster
      if (!is.null(params[[k_big]]$thetas) && !is.null(params[[k_big]]$kern_idx)) {
        current_kern_idx <- params[[k_big]]$kern_idx
        current_kern_params <- params[[k_big]]$thetas[[current_kern_idx]]
        if (!is.null(current_kern_params)) {
          kernel_params_trace[[sidx]] <- list(
            iter = iter,
            cluster = k_big,
            kern_idx = current_kern_idx,
            kern_name = kernels[[current_kern_idx]]$name,
            params = current_kern_params
          )
        }
      }
      
      # NEW: Store wavelet coefficients for the largest cluster
      if (!is.null(params[[k_big]]$wpar) && !is.null(params[[k_big]]$wpar$gamma_ch)) {
        # Extract wavelet coefficient statistics for each channel
        coeff_stats <- list()
        for (m in 1:length(params[[k_big]]$wpar$gamma_ch)) {
          gamma_m <- params[[k_big]]$wpar$gamma_ch[[m]]
          coeff_stats[[paste0("channel_", m)]] <- list(
            active_coeffs = sum(gamma_m),
            total_coeffs = length(gamma_m),
            active_ratio = sum(gamma_m) / length(gamma_m),
            mean_active = mean(gamma_m),
            sd_active = sd(gamma_m)
          )
        }
        
        wavelet_coeff_trace[[sidx]] <- list(
          iter = iter,
          cluster = k_big,
          channel = paste0("cluster_", k_big),
          coeffs = coeff_stats
        )
      }
    }

    # No progress bar needed - using single-line updates
  }

  cat(sprintf("\n✅ MCMC completed! Final K=%d, kept %d samples\n",
            length(unique(z)), nrow(Z_s))); flush.console()

  list(Z=Z_s, K=K_s, alpha=alpha_s, kern=kern_s,
       acc=acc_hist, params=params, v=v, pi=stick_to_pi(v),
       revealed_idx=revealed_idx, kernel_params_trace=kernel_params_trace,
       wavelet_coeff_trace=wavelet_coeff_trace)
}

# ------------------------ Dahl (2006) least-squares consensus partition -------

dahl_partition <- function(Z) {
  if (is.null(dim(Z)) || length(dim(Z)) != 2) stop("Z must be an S x N matrix of labels.")
  S <- nrow(Z); N <- ncol(Z)
  PSM <- matrix(0, N, N)
  for (s in 1:S) {
    zs <- Z[s, ]
    A  <- outer(zs, zs, FUN = "==") * 1L
    PSM <- PSM + A
  }
  PSM <- PSM / S
  score <- numeric(S)
  for (s in 1:S) {
    zs <- Z[s, ]
    A  <- outer(zs, zs, FUN = "==") * 1L
    score[s] <- sum((A - PSM)^2)
  }
  s_hat <- which.min(score)
  z_hat <- as.integer(factor(Z[s_hat, ]))
  K_hat <- length(unique(z_hat))
  list(z_hat = z_hat, K_hat = K_hat, s_hat = s_hat, PSM = PSM, score = score)
}

# plot_psm function moved to plotting_functions.R

dahl_from_res <- function(res) {
  if (is.null(res$Z)) stop("res$Z not found; ensure run_model() stored label draws.")
  dahl_partition(res$Z)
}

# ------------------------ MAP (modal) partition over sampled labelings --------

# relabel by order-of-appearance (canonical form), not alphabetical factor order
.canon_labels <- function(z) {
  u <- unique(z)
  match(z, u)  # 1..K in first-appearance order
}

map_partition <- function(Z) {
  if (is.null(dim(Z)) || length(dim(Z)) != 2) stop("Z must be an S x N matrix of labels.")
  S <- nrow(Z); N <- ncol(Z)
  keys <- character(S)
  for (s in 1:S) keys[s] <- paste(.canon_labels(Z[s, ]), collapse = "-")
  tab <- table(keys)
  key_hat <- names(which.max(tab))
  s_hat <- which(keys == key_hat)[1]
  z_hat <- .canon_labels(Z[s_hat, ])
  K_hat <- length(unique(z_hat))
  list(z_hat = z_hat, K_hat = K_hat, s_hat = s_hat, key_hat = key_hat, freq = max(tab))
}

map_from_res <- function(res) {
  if (is.null(res$Z)) stop("res$Z not found; ensure run_model() stored label draws.")
  map_partition(res$Z)
}
