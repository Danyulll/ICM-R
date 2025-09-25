
###############################################################################
# SEMI-SUPERVISED DP–ICM–GP WITH WALKER SLICE, RAY–MALLICK WAVELET SHRINKAGE,
# AND CARLIN–CHIB KERNEL SELECTION (Empirical Bayes warm-start for 5% normals)
###############################################################################

utils::packageVersion("utils") # ensure utils is loaded

suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(invgamma)
  library(waveslim)
  library(foreach)
  library(doParallel)
})


set.seed(42)

ensure_dyadic_J <- function(P, J) {
  if (is.null(J)) J <- log2(P)
  J_int <- as.integer(round(J))
  if (abs(P - 2^J_int) > .Machine$double.eps * max(1, P)) {
    stop(sprintf("P must be 2^J (dyadic). Got P=%s, J≈%.6f (rounded J=%d gives 2^J=%s).",
                 P, J, J_int, 2^J_int))
  }
  J_int
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

kernels <- make_kernels()



# --------------------------- Utilities ---------------------------------------

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
  K <- kronecker(B, Kx)
  P <- nrow(Kx); M <- length(eta)
  K <- K + diag(rep(eta, each=P))
  K
}

pack_L <- function(L) L[lower.tri(L, diag=TRUE)]
unpack_L <- function(theta, m) {
  L <- matrix(0, m, m)
  L[lower.tri(L, diag=TRUE)] <- theta
  diag(L) <- abs(diag(L)) + 1e-8
  L
}

# ------------------------ Wavelet wrappers (waveslim) ------------------------

wt_forward_1d <- function(y, wf="la8", J=NULL, boundary="periodic") {
  P <- length(y)
  J <- ensure_dyadic_J(P, J)
  
  w <- waveslim::dwt(y, wf=wf, n.levels=J, boundary=boundary)
  
  vec <- c(w$d1)
  idx <- list(d1 = seq_along(w$d1))
  off <- length(w$d1)
  for (lev in 2:J) {
    v <- w[[paste0("d",lev)]]
    vec <- c(vec, v)
    idx[[paste0("d",lev)]] <- (off+1):(off+length(v))
    off <- off + length(v)
  }
  vec <- c(vec, w[[paste0("s",J)]])
  idx[[paste0("s",J)]] <- (off+1):(off+length(w[[paste0("s",J)]]))
  list(coeff=vec, map=list(J=J, wf=wf, boundary=boundary, P=P, idx=idx))
}

wt_inverse_1d <- function(coeff_vec, map) {
  J <- map$J
  w <- list()
  for (lev in 1:J) w[[paste0("d",lev)]] <- coeff_vec[ map$idx[[paste0("d",lev)]] ]
  w[[paste0("s",J)]] <- coeff_vec[ map$idx[[paste0("s",J)]] ]
  attr(w, "wavelet")  <- map$wf
  attr(w, "boundary") <- map$boundary
  class(w) <- "dwt"
  waveslim::idwt(w)
}

wt_forward_mat <- function(y_mat, wf="la8", J=NULL, boundary="periodic") {
  P <- nrow(y_mat); M <- ncol(y_mat)
  out <- vector("list", M)
  for (m in 1:M) out[[m]] <- wt_forward_1d(y_mat[,m], wf=wf, J=J, boundary=boundary)
  out
}

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
      v_slab  <- g_j * v_spike
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
          ss <- ss + sum(dvec[sel]^2 / g_j)
          n_eff <- n_eff + sum(sel)
        }
        if (any(!sel)) {
          ss <- ss + sum(dvec[!sel]^2 / 1.0)
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

# ------------------------ Likelihood & MH updates ----------------------------

# ------------------------ Likelihood with tau_B -------------------------------
loglik_icm <- function(y_norm, t, L, eta, kern_cfg, kp, tau_B = 1.0) {
  # Build B as: tau_B * (shape normalized to unit trace for stability)
  M <- length(eta)
  Bshape <- tcrossprod(L)
  trB <- sum(diag(Bshape))
  if (trB > 0) Bshape <- Bshape * (M / trB)
  B  <- tau_B * Bshape
  
  Kx <- kern_cfg$fun(t, kp)
  K  <- kronecker_icm(B, Kx, eta)
  yv <- as.numeric(y_norm)
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

cc_switch_kernel <- function(kern_list, cur_idx, thetas, t, Ynorm_list, L, eta) {
  Mmod <- length(kern_list)
  p_m  <- rep(1/Mmod, Mmod)
  theta_draws <- vector("list", Mmod)
  for(m in 1:Mmod) theta_draws[[m]] <- if(m==cur_idx) thetas[[m]] else kern_list[[m]]$pstar()
  logw <- rep(NA_real_, Mmod)
  for(m in 1:Mmod) {
    kc <- kern_list[[m]]; kp <- theta_draws[[m]]
    ll_m <- sum(sapply(Ynorm_list, loglik_icm, t=t, L=L, eta=eta, kern_cfg=kc, kp=kp))
    rest <- sum(sapply(setdiff(1:Mmod, m), function(j) {
      th <- theta_draws[[j]]
      if(kern_list[[j]]$name == "Periodic") {
        dgamma(th$l_scale, 3, 2, log=TRUE) + dbeta(th$period, 5, 5, log=TRUE)
      } else {
        dgamma(th$l_scale, 2, 2, log=TRUE)
      }
    }))
    logw[m] <- log(p_m[m]) + ll_m + rest
  }
  w <- exp(logw - max(logw)); w <- w / sum(w)
  new_idx <- sample.int(Mmod, 1, prob = w)
  list(idx=new_idx, thetas=theta_draws, weights=w)
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

# ------------------------ Data generation (ICM) ------------------------------

gen_icm_curves <- function(N=50, P=128, M=3, t=seq(0,1,len=P),
                           B = matrix(c(1,.7,.5,.7,1,.4,.5,.4,1),3,3),
                           eta = rep(0.02,3),
                           kern = kernels[[1]], par = list(l_scale=0.2)) {
  Kx <- kern$fun(t, par)
  K  <- kronecker_icm(B, Kx, eta)
  Y  <- vector("list", N)
  for(i in 1:N) {
    yvec <- as.numeric(MASS::mvrnorm(1, mu=rep(0, P*M), Sigma=K + diag(1e-6, P*M)))
    
    Y[[i]] <- matrix(yvec, nrow=P, ncol=M)
  }
  Y
}

# ------------------------ Initialization -------------------------------------

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
    sigma2_i = numeric(0)
  )
}


## ------------------------ The full driver (tau_B-enabled) --------------------
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
  N <- length(Y); P <- nrow(Y[[1]]); M <- ncol(Y[[1]])
  J <- ensure_dyadic_J(P, J)

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
    
    worker_id <- Sys.getpid()
    cat(sprintf("🔄 Worker %d: Starting Empirical Bayes warm-start (%d iterations) on %d revealed normals\n",
                worker_id, emp_bayes_init_iter, length(idx)))
    flush.console()
    
    for (it in 1:emp_bayes_init_iter) {
      upd <- update_wavelet_params_rm(
        Y_list = Yk, wf=wf, J=J, boundary=boundary,
        wpar = params[[1]]$wpar, sigma2_i = params[[1]]$sigma2_i,
        a_pi=1, b_pi=1,
        g_hyp = params[[1]]$g_hyp,
        a_sig=2.5, b_sig=0.02,
        use_besov_pi = use_besov_pi,
        pi_sched = if (use_besov_pi) besov_pi_schedule(J, b=0.5, c2=besov_c2) else NULL
      )
      params[[1]]$wpar     <- upd$wpar
      params[[1]]$sigma2_i <- upd$sigma2_i

      Yk_norm <- reconstruct_normals_for_cluster(Yk, wf=wf, J=J, boundary=boundary, wpar=params[[1]]$wpar)

      cc <- cc_switch_kernel(kernels, params[[1]]$kern_idx, params[[1]]$thetas, t, Yk_norm, params[[1]]$L, params[[1]]$eta)
      params[[1]]$kern_idx <- cc$idx
      params[[1]]$thetas   <- cc$thetas

      kc <- kernels[[ params[[1]]$kern_idx ]]
      kp <- params[[1]]$thetas[[ params[[1]]$kern_idx ]]

      tmp <- mh_update_kernel(kc, kp, t, Yk_norm, params[[1]]$L, params[[1]]$eta, params[[1]]$tau_B, params[[1]]$acc)
      params[[1]]$thetas[[ params[[1]]$kern_idx ]] <- tmp$kp; params[[1]]$acc <- tmp$acc

      tmp <- mh_update_L(params[[1]]$L, mh_step_L, t, Yk_norm, params[[1]]$eta, kc,
                         params[[1]]$thetas[[ params[[1]]$kern_idx ]], params[[1]]$tau_B, params[[1]]$acc)
      params[[1]]$L   <- tmp$L; params[[1]]$acc <- tmp$acc

      tmp <- mh_update_eta(params[[1]]$eta, mh_step_eta, t, Yk_norm, params[[1]]$L, kc,
                           params[[1]]$thetas[[ params[[1]]$kern_idx ]], params[[1]]$tau_B, params[[1]]$acc)
      params[[1]]$eta <- tmp$eta; params[[1]]$acc <- tmp$acc

      # NEW: MH for tau_B
      tmp <- mh_update_tauB(params[[1]]$tau_B, mh_step_tauB, t, Yk_norm,
                            params[[1]]$L, params[[1]]$eta, kc,
                            params[[1]]$thetas[[ params[[1]]$kern_idx ]], params[[1]]$acc)
      params[[1]]$tau_B <- tmp$tau_B; params[[1]]$acc <- tmp$acc

      # NOTE: no "normalize B" rescale here anymore
      
      # Progress reporting every 100 iterations
      if (it %% 100 == 0) {
        progress_pct <- round(100 * it / emp_bayes_init_iter, 1)
        cat(sprintf("📈 Worker %d: Warm-start iteration %d/%d (%.1f%%)\n",
                    worker_id, it, emp_bayes_init_iter, progress_pct))
        flush.console()
      }
    }
    
    cat(sprintf("✅ Worker %d: Warm-start completed! Starting main MCMC...\n", worker_id))
    flush.console()
  }

  keep <- floor((n_iter - burn)/thin)
  Z_s <- matrix(NA_integer_, keep, N); K_s <- integer(keep); alpha_s <- numeric(keep); kern_s <- integer(keep)
  acc_hist <- data.frame(iter=integer(0), accL=double(0), accEta=double(0), accKer=double(0), accTauB=double(0))
  
  # Track K values for running average
  K_history <- integer(0)

  worker_id <- Sys.getpid()
  cat(sprintf("🔄 Worker %d: Starting MCMC run (N=%d, P=%d, M=%d, iters=%d, burn=%d, thin=%d, revealed=%d, K_init=%d)\n",
              worker_id, N,P,M,n_iter,burn,thin,length(revealed_idx), K))
  flush.console()

  pin_revealed <- length(revealed_idx) > 0 && !unpin_after_warmstart
  pb <- utils::txtProgressBar(min = 1, max = n_iter, style = 3,
                              title = sprintf("Worker %d MCMC", worker_id))
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
        wtf <- wt_forward_mat(Y[[i]], wf=wf, J=J, boundary=boundary)
        gam <- params[[k]]$wpar$gamma_ch
        y_norm <- wt_inverse_mat_keep(Y[[i]], wtf, gam)
        kc <- kernels[[ params[[k]]$kern_idx ]]
        kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]
        lp <- log(pi[k])
        ll <- loglik_icm(y_norm, t, params[[k]]$L, params[[k]]$eta, kc, kp, tau_B=params[[k]]$tau_B)
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
      if(length(idx) == 0) next
      Yk <- Y[idx]
      if (length(params[[k]]$sigma2_i) != length(idx) || any(!is.finite(params[[k]]$sigma2_i))) {
        params[[k]]$sigma2_i <- rep(0.05, length(idx))
      }
      upd <- update_wavelet_params_rm(
        Y_list = Yk, wf=wf, J=J, boundary=boundary,
        wpar = params[[k]]$wpar, sigma2_i = params[[k]]$sigma2_i,
        a_pi=1, b_pi=1,
        g_hyp = params[[k]]$g_hyp,
        a_sig=2.5, b_sig=0.02,
        use_besov_pi = use_besov_pi,
        pi_sched = if (use_besov_pi) besov_pi_schedule(J, b=0.5, c2=besov_c2) else NULL
      )
      params[[k]]$wpar     <- upd$wpar
      params[[k]]$sigma2_i <- upd$sigma2_i

      Yk_norm <- reconstruct_normals_for_cluster(Yk, wf=wf, J=J, boundary=boundary, wpar=params[[k]]$wpar)

      cc <- cc_switch_kernel(kernels, params[[k]]$kern_idx, params[[k]]$thetas, t, Yk_norm, params[[k]]$L, params[[k]]$eta)
      params[[k]]$kern_idx <- cc$idx
      params[[k]]$thetas   <- cc$thetas

      kc <- kernels[[ params[[k]]$kern_idx ]]
      kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]

      tmp <- mh_update_kernel(kc, kp, t, Yk_norm, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$acc)
      params[[k]]$thetas[[ params[[k]]$kern_idx ]] <- tmp$kp; params[[k]]$acc <- tmp$acc

      tmp <- mh_update_L(params[[k]]$L, mh_step_L, t, Yk_norm, params[[k]]$eta, kc,
                         params[[k]]$thetas[[ params[[k]]$kern_idx ]], params[[k]]$tau_B, params[[k]]$acc)
      params[[k]]$L   <- tmp$L; params[[k]]$acc <- tmp$acc

      tmp <- mh_update_eta(params[[k]]$eta, mh_step_eta, t, Yk_norm, params[[k]]$L, kc,
                           params[[k]]$thetas[[ params[[k]]$kern_idx ]], params[[k]]$tau_B, params[[k]]$acc)
      params[[k]]$eta <- tmp$eta; params[[k]]$acc <- tmp$acc

      # NEW: MH for tau_B
      tmp <- mh_update_tauB(params[[k]]$tau_B, mh_step_tauB, t, Yk_norm,
                            params[[k]]$L, params[[k]]$eta, kc,
                            params[[k]]$thetas[[ params[[k]]$kern_idx ]], params[[k]]$acc)
      params[[k]]$tau_B <- tmp$tau_B; params[[k]]$acc <- tmp$acc

      # NOTE: removed the old "normalize B" block entirely
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
    if(iter %% 100 == 0) {
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
      cat(sprintf("📈 Worker %d: Iteration %d/%d (%.1f%%), K=%d (avg=%.1f), alpha=%.3f\n",
                  worker_id, iter, n_iter, progress_pct, current_K, avg_K, alpha))
      flush.console()
      progress_msg <- sprintf("[%s] Worker %d MCMC: Iter %d/%d (%.1f%%), K=%d (avg=%.1f), alpha=%.3f\n",
                              format(Sys.time(), "%H:%M:%S"), worker_id, iter, n_iter,
                              progress_pct, current_K, avg_K, alpha)
      cat(progress_msg, file = progress_file, append = TRUE)
    }


    if(iter > burn && ((iter - burn) %% thin == 0)) {
      sidx <- sidx + 1
      Z_s[sidx,] <- z
      K_s[sidx]  <- length(unique(z))
      alpha_s[sidx] <- alpha
      tab <- table(z); k_big <- as.integer(names(which.max(tab)))
      kern_s[sidx] <- params[[k_big]]$kern_idx
    }

    utils::setTxtProgressBar(pb, iter)
  }
  close(pb)

  cat(sprintf("✅ Worker %d: MCMC completed! Final K=%d, kept %d samples\n",
              worker_id, length(unique(z)), nrow(Z_s)))
  flush.console()

  list(Z=Z_s, K=K_s, alpha=alpha_s, kern=kern_s,
       acc=acc_hist, params=params, v=v, pi=stick_to_pi(v),
       revealed_idx=revealed_idx)
}



###############################################################################
# Dahl (2006) least-squares consensus partition
###############################################################################

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

plot_psm <- function(PSM, z_hat = NULL, main = "Posterior similarity matrix") {
  N <- nrow(PSM)
  if (!is.null(z_hat)) {
    ord <- order(z_hat)
    PSMo <- PSM[ord, ord]
  } else {
    PSMo <- PSM
  }
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar = c(4, 4, 3, 2))
  heatmap(PSMo, Rowv = NA, Colv = NA, scale = "none",
          col = gray.colors(256, start = 1, end = 0), margins = c(5,5),
          main = main)
}

dahl_from_res <- function(res) {
  if (is.null(res$Z)) stop("res$Z not found; ensure run_model() stored label draws.")
  dahl_partition(res$Z)
}

###############################################################################
# MAP (modal) partition over sampled labelings (up to relabeling)
###############################################################################

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
  key_map <- names(which.max(tab))
  # pick the first draw achieving the mode key
  s_map <- which(keys == key_map)[1]
  z_map <- .canon_labels(Z[s_map, ])
  K_map <- length(unique(z_map))
  list(z_map = z_map, K_map = K_map, s_map = s_map, counts = tab)
}

map_from_res <- function(res) {
  if (is.null(res$Z)) stop("res$Z not found; ensure run_model() stored label draws.")
  map_partition(res$Z)
}


# ======================= Anomaly generators & dataset ========================

.pick_win <- function(P, min_len = max(4, ceiling(P*0.05)), max_len = max(6, ceiling(P*0.20))) {
  L <- sample(min_len:max_len, 1)
  s <- sample(1:(P-L+1), 1)
  list(s=s, e=s+L-1, L=L)
}

anom_spikes <- function(y, k = 5, amp = NULL) {  # Increased k from 3 to 5
  P <- nrow(y); M <- ncol(y)
  if (is.null(amp)) amp <- 5 * sd(as.vector(y))  # Increased from 3x to 5x
  for (u in 1:k) {
    t0 <- sample(1:P, 1); m <- sample(1:M, 1)
    y[t0, m] <- y[t0, m] + rnorm(1, 0, amp)
  }
  y
}

anom_burst <- function(y, amp = NULL, f = NULL) {
  P <- nrow(y); M <- ncol(y); w <- .pick_win(P)
  if (is.null(amp)) amp <- 1.5 * sd(as.vector(y))
  if (is.null(f))   f   <- runif(1, 4, 12)
  tt <- seq(0, 1, length.out = w$L)
  msel <- sample(1:M, sample(1:M,1))
  for (m in msel) y[w$s:w$e, m] <- y[w$s:w$e, m] + amp * sin(2*pi*f*tt)
  y
}

anom_varburst <- function(y, factor = 5) {  # Increased from 3 to 5
  P <- nrow(y); M <- ncol(y); w <- .pick_win(P)
  add <- matrix(rnorm(w$L*M, 0, factor * sd(as.vector(y))), nrow = w$L, ncol = M)
  y[w$s:w$e, ] <- y[w$s:w$e, ] + add
  y
}

anom_step <- function(y, delta = NULL) {
  P <- nrow(y); M <- ncol(y); w <- .pick_win(P)
  if (is.null(delta)) delta <- rnorm(M, 0, 4*apply(y,2,sd))  # Increased from 2x to 4x
  y[w$s:w$e, ] <- sweep(y[w$s:w$e, ], 2, delta, "+")
  y
}

anom_ramp <- function(y, slope = NULL) {
  P <- nrow(y); M <- ncol(y); w <- .pick_win(P)
  if (is.null(slope)) slope <- rnorm(M, 0, 3*apply(y,2,sd)/w$L)
  ramp <- outer(seq(0, 1, length.out = w$L), slope, "*")
  y[w$s:w$e, ] <- y[w$s:w$e, ] + ramp
  y
}

anom_freqchange <- function(y, amp = NULL) {
  P <- nrow(y); M <- ncol(y); w <- .pick_win(P)
  if (is.null(amp)) amp <- apply(y,2,sd)
  f_new <- runif(1, 6, 14)
  tt <- seq(0, 1, length.out = w$L)
  msel <- sample(1:M, sample(1:M,1))
  for (m in msel) y[w$s:w$e, m] <- y[w$s:w$e, m] + amp[m] * sin(2*pi*f_new*tt)
  y
}

anom_phase <- function(y, shift = NULL) {
  P <- nrow(y); M <- ncol(y); w <- .pick_win(P)
  if (is.null(shift)) shift <- runif(1, 0.2, 0.8)
  Ls <- floor(shift * w$L)
  for (m in 1:M) {
    seg <- y[w$s:w$e, m]
    if (Ls>0 && Ls < length(seg)) y[w$s:w$e, m] <- c(seg[(Ls+1):length(seg)], seg[1:Ls])
  }
  y
}

anom_corrbreak <- function(y, theta = NULL) {
  P <- nrow(y); M <- ncol(y); if (M < 2) return(y)
  w <- .pick_win(P); if (is.null(theta)) theta <- runif(1, pi/6, pi/3)
  R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  idx <- 1:min(2, M)
  y[w$s:w$e, idx] <- as.matrix(y[w$s:w$e, idx]) %*% R
  y
}

anom_warp <- function(y, factor = NULL) {
  P <- nrow(y); M <- ncol(y); w <- .pick_win(P)
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
  P <- nrow(y); M <- ncol(y); w <- .pick_win(P)
  y[w$s:w$e, ] <- matrix(y[max(1, w$s-1), ], nrow = w$L, ncol = M, byrow = TRUE)
  y
}

anom_swap <- function(y) {
  P <- nrow(y); M <- ncol(y); if (M < 2) return(y)
  w <- .pick_win(P); ch <- sample(1:M, 2)
  tmp <- y[w$s:w$e, ch[1]]
  y[w$s:w$e, ch[1]] <- y[w$s:w$e, ch[2]]
  y[w$s:w$e, ch[2]] <- tmp
  y
}

anom_global_scale <- function(y, factor = NULL) {
  if (is.null(factor)) factor <- runif(1, 1.5, 2.2)
  y * factor
}

ANOM_REG <- list(
  spikes = anom_spikes,
  burst  = anom_burst,
  varb   = anom_varburst,
  step   = anom_step,
  ramp   = anom_ramp,
  fchg   = anom_freqchange,
  phase  = anom_phase,
  corr   = anom_corrbreak,
  warp   = anom_warp,
  drop   = anom_dropout,
  swap   = anom_swap,
  gsc    = anom_global_scale
)

make_anomaly_dataset <- function(N = 200, P = 128, M = 3, t = seq(0,1,len=P),
                                 base_B = diag(M), base_eta = rep(0.02, M),
                                 base_kern = kernels[[1]], base_par = list(l_scale=0.2),
                                 anom_rates = NULL, frac_anom = 0.05, seed = 1) {
  set.seed(seed)
  Y <- gen_icm_curves(N=N, P=P, M=M, t=t, B=base_B, eta=base_eta, kern=base_kern, par=base_par)
  y_all <- Y
  ylab  <- rep(0L, N)   # 0=normal, >0 anomalous type id
  n_anom <- floor(frac_anom * N)
  idx_anom <- if (n_anom > 0) sort(sample(1:N, n_anom)) else integer(0)
  
  types <- names(ANOM_REG)
  if (is.null(anom_rates)) {
    w <- rep(1/length(types), length(types)); names(w) <- types
  } else {
    stopifnot(all(names(anom_rates) %in% types))
    w <- rep(1e-6, length(types)); names(w) <- types
    w[names(anom_rates)] <- anom_rates
    w <- w / sum(w)
  }
  for (i in idx_anom) {
    ty <- sample(types, 1, prob = w)
    y_all[[i]] <- ANOM_REG[[ty]](y_all[[i]])
    ylab[i] <- match(ty, types)
  }
  list(Y = y_all, labels = ylab, types = types, idx_anom = idx_anom)
}

# ------------------------ Plotting helpers (multi-channel, robust) -----------

.palette_k <- function(K) {
  base <- c("#1f77b4","#d62728","#2ca02c","#9467bd","#ff7f0e",
            "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
  if (K <= length(base)) base[1:K] else grDevices::rainbow(K)
}

.channel_limits <- function(Y) {
  M <- ncol(Y[[1]])
  sapply(1:M, function(m) {
    y <- unlist(lapply(Y, function(y) y[,m]))
    y <- y[is.finite(y)]
    if (!length(y)) return(c(0,1))
    rng <- range(y)
    if (rng[1] == rng[2]) rng <- rng + c(-1,1) * (ifelse(abs(rng[1]) > 0, abs(rng[1])*0.05, 1))
    rng
  })
}

.save_dataset_png <- function(Y, t, idx_anom = integer(0),
                              reveal_idx = integer(0),
                              title = "Dataset",
                              outfile = "dataset.png") {
  M <- ncol(Y[[1]]); lims <- .channel_limits(Y)
  H <- max(400, 260*M)
  png(outfile, width = 1600, height = H, res = 150)
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(M,1), mar = c(3.5,4,2,1), mgp = c(2.1,0.6,0))
  idx_norm <- setdiff(seq_along(Y), idx_anom)
  idx_unrev_norm <- setdiff(idx_norm, reveal_idx)
  for (m in 1:M) {
    plot(NA, xlim = range(t), ylim = lims[,m], xlab = "t", ylab = paste0("channel ", m),
         main = if (m==1) title else NULL)
    for (i in idx_unrev_norm) lines(t, Y[[i]][,m], col = adjustcolor("gray30", 0.20), lwd = 1)
    if (length(reveal_idx)) for (i in reveal_idx) lines(t, Y[[i]][,m], col = adjustcolor("#1f77b4", 0.90), lwd = 2)
    if (length(idx_anom)) for (i in idx_anom) lines(t, Y[[i]][,m], col = adjustcolor("#d62728", 0.75), lwd = 2)
    if (m==1) {
      leg <- c("unrevealed normal","revealed normal","anomaly")
      col <- c("gray30","#1f77b4","#d62728")
      keep <- c(length(idx_unrev_norm)>0, length(reveal_idx)>0, length(idx_anom)>0)
      legend("topright", leg[keep], col = col[keep], lwd = 2, bty="n", cex=0.9)
    }
  }
  par(op); dev.off()
}

.save_clustered_png <- function(Y, t, z_hat,
                                title = "Clustered curves (Dahl)",
                                outfile = "clustered.png") {
  M <- ncol(Y[[1]])
  K <- length(unique(z_hat))
  cols <- .palette_k(K); col_map <- cols[ as.integer(factor(z_hat)) ]
  lims <- .channel_limits(Y)
  H <- max(400, 260*M)
  png(outfile, width = 1600, height = H, res = 150)
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(M,1), mar = c(3.5,4,2,1), mgp = c(2.1,0.6,0))
  for (m in 1:M) {
    plot(NA, xlim = range(t), ylim = lims[,m], xlab = "t", ylab = paste0("channel ", m),
         main = if (m==1) title else NULL)
    for (i in seq_along(Y)) {
      lines(t, Y[[i]][,m], col = adjustcolor(col_map[i], 0.5), lwd = 1.5)
    }
    if (m==1) {
      legend("topright",
             legend = paste("Cluster", sort(unique(as.integer(factor(z_hat))))),
             col = cols[sort(unique(as.integer(factor(z_hat))))], lwd = 2, bty = "n", cex = 0.9)
    }
  }
  par(op); dev.off()
}

.save_confusion_png <- function(true_anom, pred_anom,
                                outfile = "confusion.png",
                                title = "Confusion: anomaly detection") {
  stopifnot(length(true_anom) == length(pred_anom))
  tbl <- table(True = factor(true_anom, levels = c(0,1), labels = c("Normal","Anomaly")),
               Pred = factor(pred_anom, levels = c(0,1), labels = c("Normal","Anomaly")))
  Z <- as.matrix(tbl); Zimg <- Z[2:1, , drop=FALSE]
  acc <- if (sum(Z) > 0) sum(diag(Z)) / sum(Z) else NA_real_
  png(outfile, width = 1200, height = 800, res = 150)
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)
  par(mfrow = c(1,2), mar = c(4,4,3,1), mgp = c(2.1,0.7,0))
  image(x = 1:2, y = 1:2, z = Zimg, axes = FALSE, col = hcl.colors(50, "YlOrRd"))
  axis(1, at = 1:2, labels = colnames(Z))
  axis(2, at = 1:2, labels = rev(rownames(Z)))
  title(main = title)
  mtext(sprintf("Accuracy = %s", if (is.na(acc)) "n/a" else sprintf("%.3f", acc)),
        side = 3, line = 0.2, cex = 0.9)
  barplot(Z, beside = TRUE, legend.text = TRUE,
          args.legend = list(bty="n"), main = "Counts", ylab = "N")
}

.smooth <- function(x, k = 5) { if (length(x) < k) x else as.numeric(stats::filter(x, rep(1/k, k), sides = 2)) }

.save_diagnostics <- function(res, dahl, outfile_prefix = "run",
                              true_anom = NULL, pred_anom = NULL) {
  Ktr  <- res$K;  Ktr <- Ktr[is.finite(Ktr)]
  atr  <- res$alpha; atr <- atr[is.finite(atr)]
  ktr  <- res$kern;  ktr <- ktr[is.finite(ktr)]
  acc  <- res$acc
  PSM  <- dahl$PSM
  zhat <- dahl$z_hat
  
  png(sprintf("%s_diagnostics.png", outfile_prefix), width = 1800, height = 1200, res = 150)
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(2,3), mar = c(3.5,4,2.5,1), mgp = c(2.1,0.6,0))
  if (length(Ktr)) plot(seq_along(Ktr), Ktr, type="s", xlab="kept iter", ylab="K", main="Trace of K") else { plot.new(); title("Trace of K (no draws)") }
  if (length(atr)) plot(seq_along(atr), atr, type="l", xlab="kept iter", ylab=expression(alpha), main="Trace of alpha") else { plot.new(); title("Trace of alpha (no draws)") }
  if (length(ktr)) plot(seq_along(ktr), ktr, type="s", xlab="kept iter", ylab="kernel idx", main="Kernel index (largest cluster)") else { plot.new(); title("Kernel index (no draws)") }
  if (!is.null(acc) && nrow(acc) > 0) {
    matplot(acc$iter, cbind(.smooth(acc$accL,7), .smooth(acc$accEta,7), .smooth(acc$accKer,7)),
            type="l", lty=1, lwd=2, xlab="iter", ylab="acceptance", main="MH acceptance (smoothed)")
    legend("bottomright", c("L (B)","eta","kernel params"), col=1:3, lty=1, bty="n", cex=0.9)
  } else { plot.new(); title("MH acceptance (no samples)") }
  if (length(Ktr)) {
    if (length(unique(Ktr)) > 1) { d <- density(Ktr); plot(d, main="Density of K", xlab="K") }
    else { hist(Ktr, breaks=(min(Ktr)-0.5):(max(Ktr)+0.5), main="Density of K", xlab="K") }
  } else { plot.new(); title("Density of K (no draws)") }
  if (!is.null(PSM) && all(is.finite(PSM))) {
    ord <- order(zhat); PSMo <- PSM[ord, ord, drop=FALSE]
    image(seq_len(nrow(PSMo)), seq_len(ncol(PSMo)), PSMo[rev(seq_len(nrow(PSMo))), ], col = gray.colors(256, start=1, end=0),
          xlab = "", ylab = "", main = "PSM (ordered by Dahl)"); box()
  } else { plot.new(); title("PSM unavailable") }
  dev.off()
  
  if (length(Ktr) && length(atr)) {
    K_mode <- as.integer(names(which.max(table(Ktr))))
    cat(sprintf("K_mean=%.2f, K_mode=%d, alpha_mean=%.3f\n", mean(Ktr), K_mode, mean(atr)),
        file = sprintf("%s_diag_summary.txt", outfile_prefix))
  } else {
    cat("No kept draws; diagnostics incomplete\n", file = sprintf("%s_diag_summary.txt", outfile_prefix))
  }
}

# ------------------------ Convergence Diagnostics ---------------------------

.save_convergence_diagnostics <- function(res, outfile_prefix = "run") {
  # Extract parameter traces from the largest cluster
  if (is.null(res$params) || length(res$params) == 0) {
    cat("No parameter information available for convergence diagnostics\n")
    return(invisible(NULL))
  }
  
  # Find the largest cluster (most samples)
  if (is.null(res$Z) || nrow(res$Z) == 0) {
    cat("No cluster assignments available for convergence diagnostics\n")
    return(invisible(NULL))
  }
  
  # Get the most frequent cluster across all samples
  all_assignments <- as.vector(res$Z)
  cluster_counts <- table(all_assignments)
  largest_cluster <- as.numeric(names(which.max(cluster_counts)))
  
  # Extract parameters for the largest cluster
  cluster_params <- res$params[[largest_cluster]]
  if (is.null(cluster_params)) {
    cat("No parameters available for largest cluster\n")
    return(invisible(NULL))
  }
  
  # Create convergence diagnostics plot
  png(sprintf("%s_convergence.png", outfile_prefix), width = 2000, height = 1600, res = 150)
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(3,4), mar = c(3,3,2,1), mgp = c(1.5,0.5,0))
  
  # 1. ACF plots for key parameters
  if (!is.null(res$K) && length(res$K) > 10) {
    acf(res$K, main = "ACF: K", lag.max = min(50, length(res$K)/4))
  } else {
    plot.new(); title("ACF: K (insufficient data)")
  }
  
  if (!is.null(res$alpha) && length(res$alpha) > 10) {
    acf(res$alpha, main = "ACF: alpha", lag.max = min(50, length(res$alpha)/4))
  } else {
    plot.new(); title("ACF: alpha (insufficient data)")
  }
  
  if (!is.null(res$kern) && length(res$kern) > 10) {
    acf(res$kern, main = "ACF: kernel index", lag.max = min(50, length(res$kern)/4))
  } else {
    plot.new(); title("ACF: kernel index (insufficient data)")
  }
  
  # 2. Trace plots for key parameters
  if (!is.null(res$K) && length(res$K) > 0) {
    plot(seq_along(res$K), res$K, type="l", xlab="iteration", ylab="K", main="Trace: K")
    abline(h=mean(res$K, na.rm=TRUE), col="red", lty=2)
  } else {
    plot.new(); title("Trace: K (no data)")
  }
  
  if (!is.null(res$alpha) && length(res$alpha) > 0) {
    plot(seq_along(res$alpha), res$alpha, type="l", xlab="iteration", ylab="alpha", main="Trace: alpha")
    abline(h=mean(res$alpha, na.rm=TRUE), col="red", lty=2)
  } else {
    plot.new(); title("Trace: alpha (no data)")
  }
  
  if (!is.null(res$kern) && length(res$kern) > 0) {
    plot(seq_along(res$kern), res$kern, type="l", xlab="iteration", ylab="kernel", main="Trace: kernel index")
    abline(h=mean(res$kern, na.rm=TRUE), col="red", lty=2)
  } else {
    plot.new(); title("Trace: kernel index (no data)")
  }
  
  # 3. Wavelet parameters (if available)
  if (!is.null(cluster_params$wpar) && !is.null(cluster_params$wpar$pi_level)) {
    pi_vals <- unlist(cluster_params$wpar$pi_level)
    if (length(pi_vals) > 0) {
      barplot(pi_vals, main="Wavelet pi_level", ylab="probability", las=2)
    } else {
      plot.new(); title("Wavelet pi_level (no data)")
    }
  } else {
    plot.new(); title("Wavelet pi_level (no data)")
  }
  
  if (!is.null(cluster_params$wpar) && !is.null(cluster_params$wpar$g_level)) {
    g_vals <- unlist(cluster_params$wpar$g_level)
    if (length(g_vals) > 0) {
      barplot(g_vals, main="Wavelet g_level", ylab="scale factor", las=2)
    } else {
      plot.new(); title("Wavelet g_level (no data)")
    }
  } else {
    plot.new(); title("Wavelet g_level (no data)")
  }
  
  # 4. Kernel parameters
  if (!is.null(cluster_params$thetas) && length(cluster_params$thetas) > 0) {
    # Get the current kernel parameters
    kern_idx <- cluster_params$kern_idx
    if (!is.null(kern_idx) && kern_idx <= length(cluster_params$thetas)) {
      current_params <- cluster_params$thetas[[kern_idx]]
      if (!is.null(current_params)) {
        param_vals <- unlist(current_params)
        if (length(param_vals) > 0) {
          barplot(param_vals, main="Kernel parameters", ylab="value", las=2)
        } else {
          plot.new(); title("Kernel parameters (no data)")
        }
      } else {
        plot.new(); title("Kernel parameters (no data)")
      }
    } else {
      plot.new(); title("Kernel parameters (no data)")
    }
  } else {
    plot.new(); title("Kernel parameters (no data)")
  }
  
  # 5. Correlation matrix L (if available)
  if (!is.null(cluster_params$L)) {
    L_mat <- cluster_params$L
    if (nrow(L_mat) > 1 && ncol(L_mat) > 1) {
      image(1:nrow(L_mat), 1:ncol(L_mat), L_mat, main="Correlation matrix L", 
            xlab="", ylab="", col=heat.colors(20))
    } else {
      plot.new(); title("Correlation matrix L (too small)")
    }
  } else {
    plot.new(); title("Correlation matrix L (no data)")
  }
  
  # 6. Noise parameters eta
  if (!is.null(cluster_params$eta) && length(cluster_params$eta) > 0) {
    barplot(cluster_params$eta, main="Noise parameters eta", ylab="value", 
            names.arg=paste0("eta", seq_along(cluster_params$eta)))
  } else {
    plot.new(); title("Noise parameters eta (no data)")
  }
  
  # 7. Tau_B parameter
  if (!is.null(cluster_params$tau_B)) {
    barplot(cluster_params$tau_B, main="Tau_B parameter", ylab="value")
  } else {
    plot.new(); title("Tau_B parameter (no data)")
  }
  
  dev.off()
  
  # Save convergence summary
  summary_file <- sprintf("%s_convergence_summary.txt", outfile_prefix)
  cat("=== CONVERGENCE DIAGNOSTICS SUMMARY ===\n", file = summary_file)
  cat(sprintf("Largest cluster: %d\n", largest_cluster), file = summary_file, append = TRUE)
  cat(sprintf("Cluster size: %d samples\n", cluster_counts[largest_cluster]), file = summary_file, append = TRUE)
  
  if (!is.null(res$K) && length(res$K) > 0) {
    cat(sprintf("K - Mean: %.2f, SD: %.2f, Range: [%d, %d]\n", 
                mean(res$K, na.rm=TRUE), sd(res$K, na.rm=TRUE), 
                min(res$K, na.rm=TRUE), max(res$K, na.rm=TRUE)), 
        file = summary_file, append = TRUE)
  }
  
  if (!is.null(res$alpha) && length(res$alpha) > 0) {
    cat(sprintf("Alpha - Mean: %.3f, SD: %.3f, Range: [%.3f, %.3f]\n", 
                mean(res$alpha, na.rm=TRUE), sd(res$alpha, na.rm=TRUE), 
                min(res$alpha, na.rm=TRUE), max(res$alpha, na.rm=TRUE)), 
        file = summary_file, append = TRUE)
  }
  
  if (!is.null(cluster_params$tau_B)) {
    cat(sprintf("Tau_B: %.3f\n", cluster_params$tau_B), file = summary_file, append = TRUE)
  }
  
  if (!is.null(cluster_params$eta)) {
    cat("Eta values: ", paste(round(cluster_params$eta, 4), collapse=", "), "\n", 
        file = summary_file, append = TRUE)
  }
}

###############################################################################
# BATCH (each dataset in its own folder) + SEMI-SUPERVISED REVEAL = 5% normals
###############################################################################

dir.create("anomaly_outputs_s3", showWarnings = FALSE)

# ---- run one anomaly type and save everything into its own subfolder --------
run_one_anomaly_type <- function(
    type,
    # dataset size / structure
    N = 200, P = 32, M = 3, t = seq(0,1,length.out=P),
    frac_anom = 0.05,
    # dataset GP base (for normals)
    base_B   = diag(M),
    base_eta = rep(0.02, M),
    base_kern = kernels[[1]],
    base_par  = list(l_scale = 0.18),
    # IO
    outroot = "anomaly_outputs_s3",
    seed = 101,
    # reveal policy
    reveal_frac = 0.05,
    # model iterations
    n_iter = 2000, burn = 1000, thin = 2,
    # model (DP, wavelet, kernels, MH, etc.)
    alpha_prior = c(10, 1),
    wf = "la8", J = log2(P), boundary = "periodic",
    mh_step_L = 0.03, mh_step_eta = 0.10,
    use_besov_pi = TRUE, use_besov_g = TRUE,
    emp_bayes_init_iter = 150,
    unpin_after_warmstart = FALSE,
    K_init = 5,
    besov_c2 = 0.5
) {
  stopifnot(type %in% names(ANOM_REG))
  stopifnot(P == 2^round(log2(P)))  # dyadic length
  
  # make a per-dataset subfolder
  tag <- sprintf("%s_N%d_P%d_M%d_frac%02d_seed%d",
                 type, N, P, M, round(100*frac_anom), seed)
  subdir <- file.path(outroot, tag)
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  # only this anomaly active
  rates <- rep(0, length(ANOM_REG)); names(rates) <- names(ANOM_REG)
  rates[type] <- 1
  
  # dataset
  dat <- make_anomaly_dataset(
    N = N, P = P, M = M, t = t,
    base_B = base_B, base_eta = base_eta,
    base_kern = base_kern, base_par = base_par,
    anom_rates = rates, frac_anom = frac_anom, seed = seed
  )
  Y <- dat$Y
  idx_anom <- dat$idx_anom
  
  # choose reveal_frac of NORMALS to reveal (never anomalies)
  idx_normals <- setdiff(seq_len(N), idx_anom)
  n_reveal <- max(1, floor(reveal_frac * length(idx_normals)))
  set.seed(seed + 999)
  reveal_idx <- sort(sample(idx_normals, n_reveal))
  
  prefix <- file.path(subdir, "gen")
  
  # dataset plot (mark revealed normals)
  .save_dataset_png(
    Y, t, idx_anom, reveal_idx = reveal_idx,
    title = sprintf("Dataset (%.0f%% normal / %.0f%% %s)   Revealed normals = %d",
                    100*(1-frac_anom), 100*frac_anom, type, length(reveal_idx)),
    outfile = sprintf("%s_dataset.png", prefix)
  )
  
  # fit model (semi-supervised)
  set.seed(seed + which(names(ANOM_REG) == type))
  res <- run_model(
    Y, t,
    n_iter = n_iter, burn = burn, thin = thin,
    alpha_prior = alpha_prior,
    wf = wf, J = J, boundary = boundary,
    mh_step_L = mh_step_L, mh_step_eta = mh_step_eta,
    use_besov_pi = use_besov_pi, use_besov_g = use_besov_g,
    revealed_idx = reveal_idx,
    emp_bayes_init_iter = emp_bayes_init_iter,
    unpin_after_warmstart = unpin_after_warmstart,
    K_init = K_init,
    besov_c2 = besov_c2
  )
  
  # === Consensus (Dahl) and MAP (modal) partitions ==========================
  dahl <- dahl_from_res(res)
  z_hat <- dahl$z_hat
  
  map  <- map_from_res(res)
  z_map <- map$z_map
  
  # predicted anomaly = not in the cluster with most revealed normals (done for both)
  # Find cluster with most revealed normal samples for Dahl method
  if (length(reveal_idx) > 0) {
    revealed_clusters_hat <- z_hat[reveal_idx]
    tab_revealed_hat <- table(revealed_clusters_hat)
    normal_cluster_hat <- as.integer(names(which.max(tab_revealed_hat)))
  } else {
    # Fallback to largest cluster if no revealed samples
    tabz_hat <- table(z_hat)
    normal_cluster_hat <- as.integer(names(which.max(tabz_hat)))
  }
  pred_anom_hat <- as.integer(z_hat != normal_cluster_hat)
  
  # Find cluster with most revealed normal samples for MAP method
  if (length(reveal_idx) > 0) {
    revealed_clusters_map <- z_map[reveal_idx]
    tab_revealed_map <- table(revealed_clusters_map)
    normal_cluster_map <- as.integer(names(which.max(tab_revealed_map)))
  } else {
    # Fallback to largest cluster if no revealed samples
    tabz_map <- table(z_map)
    normal_cluster_map <- as.integer(names(which.max(tabz_map)))
  }
  pred_anom_map <- as.integer(z_map != normal_cluster_map)
  
  true_anom <- as.integer(seq_along(Y) %in% idx_anom)
  
  # save clustered plots (Dahl + MAP) and diagnostics
  .save_clustered_png(
    Y, t, z_hat,
    title = sprintf("Clustered curves (Dahl): %s anomalies", type),
    outfile = sprintf("%s_clustered_dahl.png", prefix)
  )
  .save_clustered_png(
    Y, t, z_map,
    title = sprintf("Clustered curves (MAP): %s anomalies", type),
    outfile = sprintf("%s_clustered_map.png", prefix)
  )
  
  .save_diagnostics(
    res, dahl, outfile_prefix = prefix,
    true_anom = true_anom, pred_anom = pred_anom_hat
  )
  
  # Save convergence diagnostics
  .save_convergence_diagnostics(res, outfile_prefix = prefix)
  
  # save confusion plots for both
  .save_confusion_png(
    true_anom, pred_anom_hat,
    outfile = sprintf("%s_confusion_dahl.png", prefix),
    title   = sprintf("Confusion (Dahl): Normal vs %s anomaly", type)
  )
  .save_confusion_png(
    true_anom, pred_anom_map,
    outfile = sprintf("%s_confusion_map.png", prefix),
    title   = sprintf("Confusion (MAP): Normal vs %s anomaly", type)
  )
  
  # metadata: add MAP/Dahl summaries
  meta_path <- file.path(subdir, "meta.txt")
  cat(
    sprintf(paste0(
      "type=%s\nN=%d\nP=%d\nM=%d\nfrac_anom=%.3f\nseed=%d\n",
      "n_iter=%d\nburn=%d\nthin=%d\nuse_besov_pi=%s\nuse_besov_g=%s\n",
      "revealed_normals=%s\n",
      "Dahl: K_hat=%d, s_hat=%d\n",
      "MAP:  K_map=%d, s_map=%d, mode_count=%d\n"
    ),
    type, N, P, M, frac_anom, seed,
    n_iter, burn, thin, use_besov_pi, use_besov_g,
    paste(reveal_idx, collapse=","),
    dahl$K_hat, dahl$s_hat,
    map$K_map, map$s_map, as.integer(map$counts[names(which.max(map$counts))])
    ), file = meta_path, append = TRUE)
  
  
  invisible(list(
    Y=Y, t=t, idx_anom=idx_anom, res=res,
    dahl=dahl, z_dahl=z_hat,
    map=map,  z_map=z_map,
    true_anom=true_anom,
    pred_anom_dahl=pred_anom_hat,
    pred_anom_map=pred_anom_map,
    reveal_idx=reveal_idx, outdir=subdir, prefix=prefix
  ))
  
}

# ---- run all anomaly types: each in its own subfolder -----------------------
run_all_anomaly_types <- function(
    # dataset size / structure
  N = 200, P = 32, M = 3, t = seq(0,1,length.out=P),
  frac_anom = 0.05,
  # dataset GP base (for normals)
  base_B   = diag(M),
  base_eta = rep(0.02, M),
  base_kern = kernels[[1]],
  base_par  = list(l_scale = 0.18),
  # IO
  outroot = "anomaly_outputs_s3",
  seed = 101,
  # reveal policy
  reveal_frac = 0.05,
  # model iterations
  n_iter = 2000, burn = 1000, thin = 2,
  # model (DP, wavelet, kernels, MH, etc.)
  alpha_prior = c(10, 1),
  wf = "la8", J = log2(P), boundary = "periodic",
  mh_step_L = 0.03, mh_step_eta = 0.10,
  use_besov_pi = TRUE, use_besov_g = TRUE,
  emp_bayes_init_iter = 150,
  unpin_after_warmstart = FALSE,
  K_init = 5,
  besov_c2 = 0.5,
  # parallel controls
  use_parallel = TRUE,
  n_workers = NULL  # Will be set to min(available_cores, num_tasks)
) {
  dir.create(outroot, showWarnings = FALSE, recursive = TRUE)
  outroot <- normalizePath(outroot, winslash = "/", mustWork = FALSE)
  types <- c("spikes", "step", "varb")  # Only test these three anomaly types
  
  # Set n_workers to be efficient: use min(available_cores, num_tasks)
  if (is.null(n_workers)) {
    available_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
    n_workers <- min(available_cores, length(types))
  }
  
  # -------- sequential path --------
  if (!use_parallel) {
    out <- setNames(vector("list", length(types)), types)
    for (i in seq_along(types)) {
      ty <- types[i]
      message(sprintf("=== Running anomaly type: %s ===", ty))
      out[[ty]] <- run_one_anomaly_type(
        type = ty,
        N = N, P = P, M = M, t = t,
        frac_anom = frac_anom,
        base_B = base_B, base_eta = base_eta,
        base_kern = base_kern, base_par = base_par,
        outroot = outroot,
        seed = seed + i*1000L,
        reveal_frac = reveal_frac,
        n_iter = n_iter, burn = burn, thin = thin,
        alpha_prior = alpha_prior,
        wf = wf, J = J, boundary = boundary,
        mh_step_L = mh_step_L, mh_step_eta = mh_step_eta,
        use_besov_pi = use_besov_pi, use_besov_g = use_besov_g,
        emp_bayes_init_iter = emp_bayes_init_iter,
        unpin_after_warmstart = unpin_after_warmstart,
        K_init = K_init,
        besov_c2 = besov_c2
      )
    }
    message(sprintf("All outputs saved under ./%s/<type_N_P_M_frac_seed>/", outroot))
    return(invisible(out))
  }
  
  # -------- parallel path: foreach + doParallel (Windows PSOCK) --------
  cl <- parallel::makeCluster(n_workers)
  on.exit({ 
    try(doParallel::stopImplicitCluster(), silent = TRUE)
    try(parallel::stopCluster(cl), silent = TRUE) 
  }, add = TRUE)
  doParallel::registerDoParallel(cl)
  
  # Load packages on workers
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(MASS); library(mvtnorm); library(invgamma); library(waveslim)
      library(foreach); library(doParallel)
    })
    NULL
  })
  
  # Export specific functions and variables needed by workers
  required_functions <- c(
    # Main functions
    "run_one_anomaly_type", "run_model", "make_anomaly_dataset", "gen_icm_curves",
    "update_wavelet_params_rm", "reconstruct_normals_for_cluster", "loglik_icm",
    "mh_update_kernel", "mh_update_L", "mh_update_eta", "cc_switch_kernel",
    "draw_new_cluster_params", "dahl_from_res", "map_from_res","mh_update_tauB",
    
    # Helper functions (dot-prefixed)
    ".save_dataset_png", ".save_clustered_png", ".save_confusion_png", 
    ".save_diagnostics", ".save_convergence_diagnostics", ".channel_limits", ".palette_k", ".smooth",
    ".canon_labels", ".pick_win",
    
    # Statistical functions
    "dahl_partition", "map_partition", "plot_psm",
    
    # Wavelet functions
    "wt_forward_mat", "wt_inverse_mat_keep", "wt_stack_channel", 
    "wt_forward_1d", "wt_inverse_1d",
    
    # Utility functions
    "ensure_dyadic_J", "make_kernels", "kernels", "ANOM_REG", 
    "besov_pi_schedule", "besov_g_hyper", "kronecker_icm", "chol_loglik", 
    "pack_L", "unpack_L", "logit", "invlogit", "stick_to_pi", 
    "extend_sticks_until", "update_v_given_z",
    
    # Anomaly generation functions
    "anom_spikes", "anom_burst", "anom_varburst", "anom_step", "anom_ramp",
    "anom_freqchange", "anom_phase", "anom_corrbreak", "anom_warp", 
    "anom_dropout", "anom_swap", "anom_global_scale"
  )
  
  # Progress tracking setup
  total_types <- length(types)
  cat(sprintf("🚀 Starting parallel execution with %d workers on %d anomaly types...\n", n_workers, total_types))
  cat(sprintf("📋 Processing types: %s\n", paste(types, collapse=", ")))
  flush.console()
  
  # Create a progress file to track worker status
  progress_file <- file.path(outroot, "parallel_progress.txt")
  if (file.exists(progress_file)) file.remove(progress_file)
  
  # Check which functions exist and export them
  existing_functions <- required_functions[sapply(required_functions, exists, envir = globalenv())]
  if (length(existing_functions) > 0) {
    parallel::clusterExport(cl, varlist = existing_functions, envir = globalenv())
  }
  
  # Export variables needed
  parallel::clusterExport(cl, varlist = c("N", "P", "M", "t", "frac_anom", 
                                         "base_B", "base_eta", "base_kern", "base_par",
                                         "outroot", "reveal_frac", "n_iter", "burn", "thin",
                                         "alpha_prior", "wf", "J", "boundary", 
                                         "mh_step_L", "mh_step_eta", "use_besov_pi", 
                                         "use_besov_g", "emp_bayes_init_iter", 
                                         "unpin_after_warmstart", "K_init", "besov_c2",
                                         "progress_file"),
                         envir = environment())
  
  # Run in parallel - iterate over types directly
  out_list <- foreach(
    i = seq_along(types),
    .packages = c("MASS","mvtnorm","invgamma","waveslim","foreach","doParallel"),
    .errorhandling = "stop",
    .combine = "list",
    .multicombine = TRUE
  ) %dopar% {
    ty <- types[i]
    current_seed <- seed + i*1000L
    worker_id <- Sys.getpid()
    
    # Thread-specific progress messages
    cat(sprintf("🚀 Worker %d starting anomaly type: %s (task %d/%d)\n", 
               worker_id, ty, i, total_types))
    flush.console()
    
    # Write to progress file
    progress_msg <- sprintf("[%s] Worker %d starting %s (task %d/%d)\n", 
                           format(Sys.time(), "%H:%M:%S"), worker_id, ty, i, total_types)
    cat(progress_msg, file = progress_file, append = TRUE)
    
    # Set seed for this worker
    set.seed(current_seed)
    
    # Start timing
    start_time <- Sys.time()
    
    # Run the analysis with progress updates
    result <- tryCatch({
      cat(sprintf("📊 Worker %d: Generating dataset for %s...\n", worker_id, ty))
      flush.console()
      
      # Update progress file
      progress_msg <- sprintf("[%s] Worker %d generating dataset for %s\n", 
                             format(Sys.time(), "%H:%M:%S"), worker_id, ty)
      cat(progress_msg, file = progress_file, append = TRUE)
      
      # Call the main function
      run_one_anomaly_type(
        type = ty,
        N = N, P = P, M = M, t = t,
        frac_anom = frac_anom,
        base_B = base_B, base_eta = base_eta,
        base_kern = base_kern, base_par = base_par,
        outroot = outroot,
        seed = current_seed,
        reveal_frac = reveal_frac,
        n_iter = n_iter, burn = burn, thin = thin,
        alpha_prior = alpha_prior,
        wf = wf, J = J, boundary = boundary,
        mh_step_L = mh_step_L, mh_step_eta = mh_step_eta,
        use_besov_pi = use_besov_pi, use_besov_g = use_besov_g,
        emp_bayes_init_iter = emp_bayes_init_iter,
        unpin_after_warmstart = unpin_after_warmstart,
        K_init = K_init,
        besov_c2 = besov_c2
      )
    }, error = function(e) {
      cat(sprintf("❌ Worker %d ERROR on %s: %s\n", worker_id, ty, e$message))
      flush.console()
      stop(e)
    })
    
    # Completion message with timing
    end_time <- Sys.time()
    elapsed <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
    cat(sprintf("✅ Worker %d completed %s in %.2f minutes (task %d/%d)\n", 
               worker_id, ty, elapsed, i, total_types))
    flush.console()
    
    # Extract confusion matrix information
    if (!is.null(result) && !is.null(result$true_anom) && !is.null(result$pred_anom_dahl)) {
      true_anom <- result$true_anom
      pred_anom_dahl <- result$pred_anom_dahl
      pred_anom_map <- result$pred_anom_map
      
      # Calculate confusion matrix for Dahl method
      cm_dahl <- table(True = factor(true_anom, levels = c(0,1), labels = c("Normal","Anomaly")),
                      Pred = factor(pred_anom_dahl, levels = c(0,1), labels = c("Normal","Anomaly")))
      
      # Calculate confusion matrix for MAP method
      cm_map <- table(True = factor(true_anom, levels = c(0,1), labels = c("Normal","Anomaly")),
                     Pred = factor(pred_anom_map, levels = c(0,1), labels = c("Normal","Anomaly")))
      
      # Calculate accuracy
      acc_dahl <- if (sum(cm_dahl) > 0) sum(diag(cm_dahl)) / sum(cm_dahl) else NA_real_
      acc_map <- if (sum(cm_map) > 0) sum(diag(cm_map)) / sum(cm_map) else NA_real_
      
      # Print confusion matrix to progress file
      confusion_msg <- paste0(
        sprintf("[%s] Worker %d CONFUSION MATRIX for %s:\n", 
                format(Sys.time(), "%H:%M:%S"), worker_id, ty),
        sprintf("Dahl Method - Accuracy: %.3f\n", acc_dahl),
        "  True\\Pred    Normal  Anomaly\n",
        sprintf("  Normal      %6d  %6d\n", cm_dahl[1,1], cm_dahl[1,2]),
        sprintf("  Anomaly     %6d  %6d\n", cm_dahl[2,1], cm_dahl[2,2]),
        sprintf("MAP Method - Accuracy: %.3f\n", acc_map),
        "  True\\Pred    Normal  Anomaly\n",
        sprintf("  Normal      %6d  %6d\n", cm_map[1,1], cm_map[1,2]),
        sprintf("  Anomaly     %6d  %6d\n", cm_map[2,1], cm_map[2,2]),
        "----------------------------------------\n"
      )
      cat(confusion_msg, file = progress_file, append = TRUE)
    }
    
    # Update progress file
    progress_msg <- sprintf("[%s] Worker %d completed %s in %.2f minutes (task %d/%d)\n", 
                           format(Sys.time(), "%H:%M:%S"), worker_id, ty, elapsed, i, total_types)
    cat(progress_msg, file = progress_file, append = TRUE)
    
    result
  }
  
  out <- setNames(out_list, types)
  
  # Final summary
  cat(sprintf("🎉 All parallel tasks completed! %d anomaly types processed.\n", total_types))
  cat(sprintf("📁 All outputs saved under ./%s/<type_N_P_M_frac_seed>/\n", outroot))
  
  # Show progress file contents
  if (file.exists(progress_file)) {
    cat("\n📋 Progress log:\n")
    progress_log <- readLines(progress_file)
    for (line in progress_log) {
      cat(line, "\n")
    }
    file.remove(progress_file)  # Clean up
  }
  
  flush.console()
  invisible(out)
}



# ------------------------------ Progress Monitoring --------------------------

# Simple function to monitor parallel progress in real-time
monitor_parallel_progress <- function(outroot = "anomaly_outputs_s3", refresh_seconds = 3) {
  progress_file <- file.path(outroot, "parallel_progress.txt")
  
  cat("🔍 Monitoring parallel progress... (Press Ctrl+C to stop)\n")
  cat("📁 Watching:", progress_file, "\n\n")
  
  last_line_count <- 0
  
  while (TRUE) {
    if (file.exists(progress_file)) {
      progress_log <- readLines(progress_file)
      current_line_count <- length(progress_log)
      
      # Only show new lines
      if (current_line_count > last_line_count) {
        cat("\n", rep("=", 60), "\n", sep="")
        cat("📊 Progress Update at", format(Sys.time(), "%H:%M:%S"), "\n")
        cat(rep("=", 60), "\n", sep="")
        
        # Show new lines
        new_lines <- progress_log[(last_line_count + 1):current_line_count]
        for (line in new_lines) {
          cat(line, "\n")
        }
        
        last_line_count <- current_line_count
      }
    } else {
      cat("⏳ Waiting for progress file to be created...\n")
    }
    
    Sys.sleep(refresh_seconds)
  }
}

# ------------------------------ Example --------------------------------------
# (uncomment to run a quick pass)
dim = 16
results_all <- run_all_anomaly_types(
  N = 25, P = dim, M = 3,
  base_B = diag(3),
  base_eta = c(0.01, 0.01, 0.01),          # DECREASED: Less noise for cleaner normals
  base_kern = kernels[[1]],                 # CHANGED: Use SE kernel for smoother normals
  base_par  = list(l_scale = 0.25),         # INCREASED: Longer correlation length
  outroot = "anomaly_outputs_s3",
  seed = 87460945,
  frac_anom = 0.15,
  reveal_frac = 0.50,  # INCREASED: Need more revealed normals for proper identification
  n_iter = 5000, burn = 1000, thin = 1,

  # --- TIER 1 CHANGES ---
  alpha_prior = c(15, 1),                   # MODERATE: Prefer more clusters but not too many
  besov_c2 = 0.20,                         # DECREASED: Make normal template even stiffer

  # --- TIER 2 CHANGES (Optional, try if Tier 1 is not enough) ---
  emp_bayes_init_iter = 300,               # INCREASED: Better warm-start
  K_init = 8,                              # MODERATE: Encourage exploration but not too much

  # --- Original parameters ---
  wf = "la8", J = log2(dim), boundary = "periodic",
  mh_step_L = 0.04, mh_step_eta = 0.12,
  use_besov_pi = TRUE, use_besov_g = TRUE,
  unpin_after_warmstart = FALSE
)
