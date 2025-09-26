###############################################################################
# PLOTTING FUNCTIONS FOR ICM-R ANOMALY DETECTION
# ---------------------------------------------------------------------------
# This file contains all plotting and visualization functions used in the
# ICM-R anomaly detection system, including:
# - Dataset visualization
# - Clustering results
# - Classification metrics
# - MCMC diagnostics
# - ROC curves and confusion matrices
###############################################################################

# ------------------------ Color and Utility Functions ----------------------

.palette_k <- function(K) {
  base <- c("#1f77b4","#d62728","#2ca02c","#9467bd","#ff7f0e",
            "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
  if (K <= length(base)) base[1:K] else rainbow(K)
}

.channel_limits <- function(Y) {
  M <- ncol(Y[[1]])
  sapply(1:M, function(m) {
    v <- unlist(lapply(Y, function(y) y[,m]))
    c(min(v, na.rm=TRUE), max(v, na.rm=TRUE))
  })
}

# ------------------------ Dataset Visualization ----------------------------

.save_dataset_png <- function(Y, t, idx_anom = integer(0),
                              reveal_idx = integer(0),
                              title = "Dataset",
                              outfile = "dataset.png") {
  M <- ncol(Y[[1]]); lims <- .channel_limits(Y)
  H <- max(400, 260*M)
  png(outfile, width = 1600, height = H, res = 150)
  op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)
  par(mfrow = c(M,1), mar = c(3,4,2,1), mgp = c(2.1,0.6,0))
  
  for (m in 1:M) {
    plot(NA, xlim = range(t), ylim = lims[,m], 
         xlab = if (m == M) "Time" else "", ylab = sprintf("Channel %d", m),
         main = if (m == 1) title else "")
    for (i in seq_along(Y)) {
      col <- if (i %in% idx_anom) "red" else if (i %in% reveal_idx) "blue" else "black"
      lwd <- if (i %in% reveal_idx) 2 else 1
      lines(t, Y[[i]][,m], col = col, lwd = lwd)
    }
    if (length(idx_anom) > 0) {
      legend("topright", c("Normal", "Anomaly", "Revealed"), 
             col = c("black", "red", "blue"), lwd = c(1,1,2), bty = "n")
    }
  }
  par(op); dev.off()
}

# ------------------------ Clustering Visualization -------------------------

.save_clustered_png <- function(Y, t, z_hat,
                                title = "Clustered curves (Dahl)",
                                outfile = "clustered.png") {
  M <- ncol(Y[[1]])
  K <- length(unique(z_hat))
  cols <- .palette_k(K); col_map <- cols[ as.integer(factor(z_hat)) ]
  lims <- .channel_limits(Y)
  H <- max(400, 260*M)
  png(outfile, width = 1600, height = H, res = 150)
  op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)
  par(mfrow = c(M,1), mar = c(3,4,2,1), mgp = c(2.1,0.6,0))
  
  for (m in 1:M) {
    plot(NA, xlim = range(t), ylim = lims[,m], 
         xlab = if (m == M) "Time" else "", ylab = sprintf("Channel %d", m),
         main = if (m == 1) title else "")
    for (i in seq_along(Y)) {
      lines(t, Y[[i]][,m], col = col_map[i], lwd = 1)
    }
    legend("topright", sprintf("Cluster %d", 1:K), col = cols, lwd = 2, bty = "n")
  }
  par(op); dev.off()
}

.save_binary_clustered_png <- function(Y, t, pred_anom, normal_cluster, anomaly_clusters,
                                      title = "Binary Classification: Normal vs Anomaly",
                                      outfile = "binary_clustered.png") {
  M <- ncol(Y[[1]])
  lims <- .channel_limits(Y)
  H <- max(400, 260*M)
  
  # Create color mapping
  col_map <- ifelse(pred_anom, "red", "black")
  
  png(outfile, width = 1600, height = H, res = 150)
  op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)
  par(mfrow = c(M,1), mar = c(3,4,2,1), mgp = c(2.1,0.6,0))
  
  for (m in 1:M) {
    plot(NA, xlim = range(t), ylim = lims[,m], 
         xlab = if (m == M) "Time" else "", ylab = sprintf("Channel %d", m),
         main = if (m == 1) title else "")
    for (i in seq_along(Y)) {
      lines(t, Y[[i]][,m], col = col_map[i], lwd = 1)
    }
    legend("topright", c("Normal", "Anomaly"), col = c("black", "red"), lwd = 2, bty = "n")
  }
  par(op); dev.off()
}

# ------------------------ Confusion Matrix Visualization ------------------

.save_confusion_png <- function(true_anom, pred_anom,
                                outfile = "confusion.png",
                                title = "Confusion: anomaly detection") {
  cm <- table(true_anom, pred_anom)
  png(outfile, width = 800, height = 600, res = 150)
  op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)
  
  # Create a more detailed confusion matrix plot
  par(mar = c(5,5,4,2))
  image(1:nrow(cm), 1:ncol(cm), cm, 
        col = heat.colors(max(cm)), 
        xlab = "Predicted", ylab = "True",
        main = title, axes = FALSE)
  axis(1, at = 1:ncol(cm), labels = colnames(cm))
  axis(2, at = 1:nrow(cm), labels = rownames(cm))
  
  # Add text labels
  for (i in 1:nrow(cm)) {
    for (j in 1:ncol(cm)) {
      text(j, i, cm[i,j], cex = 2, font = 2)
    }
  }
  
  par(op); dev.off()
}

# ------------------------ ROC Curve and AUC Visualization -----------------

plot_auc_roc <- function(true_labels, predicted_labels, title = "ROC Curve") {
  if (length(unique(true_labels)) <= 1 || length(unique(predicted_labels)) <= 1) {
    plot.new()
    title(main = paste(title, "(insufficient data)"))
    return(NA)
  }
  
  tryCatch({
    pred_obj <- ROCR::prediction(predicted_labels, true_labels)
    perf <- ROCR::performance(pred_obj, "tpr", "fpr")
    auc <- ROCR::performance(pred_obj, "auc")@y.values[[1]]
    
    plot(perf, main = paste(title, sprintf("(AUC = %.3f)", auc)),
         col = "blue", lwd = 2)
    abline(0, 1, lty = 2, col = "red")
    
    return(auc)
  }, error = function(e) {
    plot.new()
    title(main = paste(title, "(error computing ROC)"))
    return(NA)
  })
}

# ------------------------ MCMC Diagnostics Visualization ------------------

plot_psm <- function(PSM, z_hat = NULL, main = "Posterior similarity matrix") {
  N <- nrow(PSM)
  if (!is.null(z_hat)) {
    # Reorder by cluster assignment
    ord <- order(z_hat)
    PSM <- PSM[ord, ord]
    z_hat <- z_hat[ord]
  }
  
  image(1:N, 1:N, PSM, col = heat.colors(100), 
        xlab = "Observation", ylab = "Observation", main = main)
  
  if (!is.null(z_hat)) {
    # Add cluster boundaries
    cluster_boundaries <- which(diff(z_hat) != 0)
    abline(v = cluster_boundaries + 0.5, h = cluster_boundaries + 0.5, 
           col = "white", lwd = 2)
  }
}

# ------------------------ Diagnostic Plotting Functions -------------------

.smooth <- function(x, k = 5) { 
  if (length(x) < k) x else as.numeric(stats::filter(x, rep(1/k, k), sides = 2)) 
}

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
    plot(acc$iter, .smooth(acc$accL), type="l", xlab="iter", ylab="acc rate", main="Acceptance: L", ylim=c(0,1))
    plot(acc$iter, .smooth(acc$accEta), type="l", xlab="iter", ylab="acc rate", main="Acceptance: eta", ylim=c(0,1))
    plot(acc$iter, .smooth(acc$accKer), type="l", xlab="iter", ylab="acc rate", main="Acceptance: kernel", ylim=c(0,1))
  } else {
    plot.new(); title("Acceptance rates (no data)")
    plot.new(); title("Acceptance rates (no data)")
    plot.new(); title("Acceptance rates (no data)")
  }
  par(op); dev.off()
  
  # Save diagnostic summary
  diag_file <- sprintf("%s_diag_summary.txt", outfile_prefix)
  cat("=== MCMC DIAGNOSTICS SUMMARY ===\n", file = diag_file)
  cat(sprintf("Total iterations: %d\n", nrow(res$Z)), file = diag_file, append = TRUE)
  cat(sprintf("Final K: %d\n", dahl$K_hat), file = diag_file, append = TRUE)
  cat(sprintf("Mean alpha: %.4f\n", mean(atr)), file = diag_file, append = TRUE)
  if (!is.null(acc) && nrow(acc) > 0) {
    cat(sprintf("Mean acceptance L: %.3f\n", mean(acc$accL, na.rm=TRUE)), file = diag_file, append = TRUE)
    cat(sprintf("Mean acceptance eta: %.3f\n", mean(acc$accEta, na.rm=TRUE)), file = diag_file, append = TRUE)
    cat(sprintf("Mean acceptance kernel: %.3f\n", mean(acc$accKer, na.rm=TRUE)), file = diag_file, append = TRUE)
  }
  invisible(NULL)
}

# ------------------------ Wavelet Coefficient Diagnostics -----------------

.save_wavelet_coeff_diagnostics <- function(res, outfile_prefix = "run", metric = "active_ratio") {
  td <- extract_wavelet_coeff_traces(res, metric = metric)
  if (is.null(td) || is.null(td$channel_traces) || length(td$channel_traces) == 0) {
    cat("No wavelet coefficient traces available for diagnostics\n")
    return(invisible(NULL))
  }

  png(sprintf("%s_wavelet_coeff_diagnostics.png", outfile_prefix),
      width = 2000, height = 1600, res = 150)
  op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)

  n <- length(td$channel_traces)
  n_cols <- min(4, n); n_rows <- ceiling(n / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(3,3,3,1), mgp = c(1.5,0.5,0))

  for (nm in names(td$channel_traces)) {
    v <- as.numeric(td$channel_traces[[nm]])
    v <- v[is.finite(v)]
    if (length(v) < 5) { plot.new(); title(main = sprintf("%s (%s: n<5)", nm, metric)); next }

    # Trace
    plot(seq_along(v), v, type = "l",
         xlab = "iteration", ylab = metric,
         main = sprintf("Trace: %s", nm))
    abline(h = mean(v), lty = 2)

    # ACF (second panel per facet)
    if (length(v) > 20) {
      acf(v, main = sprintf("ACF: %s", nm), lag.max = min(50, length(v)/4))
    } else {
      plot.new(); title(main = sprintf("ACF: %s (n<20)", nm))
    }
  }

  # optional summary file
  summary_file <- sprintf("%s_wavelet_coeff_summary.txt", outfile_prefix)
  cat("=== WAVELET COEFFICIENT DIAGNOSTICS SUMMARY ===\n", file = summary_file)
  cat(sprintf("Metric: %s\nChannels: %s\n",
              td$metric, paste(names(td$channel_traces), collapse = ", ")),
      file = summary_file, append = TRUE)
  invisible(NULL)
}

# ------------------------ Kernel Parameter Diagnostics ---------------------

.save_kernel_param_diagnostics <- function(res, outfile_prefix = "run") {
  if (is.null(res$kernel_params_trace) || length(res$kernel_params_trace) == 0) {
    cat("No kernel parameter traces available for diagnostics\n")
    return(invisible(NULL))
  }
  td <- extract_kernel_param_traces(res)
  if (is.null(td$param_traces) || length(td$param_traces) == 0) {
    cat("No valid kernel parameter traces found\n")
    return(invisible(NULL))
  }

  png(sprintf("%s_kernel_param_diagnostics.png", outfile_prefix),
      width = 2000, height = 1600, res = 150)
  op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)

  n <- length(td$param_traces)
  n_cols <- min(4, n); n_rows <- ceiling(n / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(3,3,3,1), mgp = c(1.5,0.5,0))

  for (param_name in names(td$param_traces)) {
    x <- as.numeric(td$param_traces[[param_name]])
    x <- x[is.finite(x)]
    if (length(x) < 5) { plot.new(); title(main = sprintf("%s (n<5)", param_name)); next }

    # Trace
    plot(seq_along(x), x, type = "l",
         xlab = "iteration", ylab = "value",
         main = sprintf("Trace: %s", param_name))
    abline(h = mean(x), lty = 2)

    # ACF
    if (length(x) > 20) {
      acf(x, main = sprintf("ACF: %s", param_name), lag.max = min(50, length(x)/4))
    } else {
      plot.new(); title(main = sprintf("ACF: %s (n<20)", param_name))
    }
  }

  # optional summary file
  summary_file <- sprintf("%s_kernel_param_summary.txt", outfile_prefix)
  cat("=== KERNEL PARAMETER DIAGNOSTICS SUMMARY ===\n", file = summary_file)
  cat(sprintf("Parameters: %s\n", paste(names(td$param_traces), collapse = ", ")), file = summary_file, append = TRUE)
  cat(sprintf("Kernels: %s\n", paste(td$kernel_names, collapse = ", ")), file = summary_file, append = TRUE)
  invisible(NULL)
}

# ------------------------ Coda Diagnostics --------------------------------

.save_coda_diagnostics <- function(res, outfile_prefix = "run") {
  if (is.null(res) || length(res) == 0) {
    cat("No results available for coda diagnostics\n")
    return(invisible(NULL))
  }

  # Helper function to safely create mcmc objects
  safe_mcmc <- function(data, name) {
    if (is.null(data) || length(data) == 0) return(NULL)
    clean_data <- data[is.finite(data)]
    if (length(clean_data) < 2) {
      cat(sprintf("Warning: Insufficient data in %s (length=%d), skipping\n", name, length(clean_data)))
      return(NULL)
    }
    
    tryCatch({
      coda::mcmc(clean_data)
    }, error = function(e) {
      cat(sprintf("Warning: Failed to create mcmc object for %s: %s\n", name, e$message))
      NULL
    })
  }

  # Extract and clean main parameters
  mcmc_list <- list()
  if (!is.null(res$K)) {
    mcmc_obj <- safe_mcmc(res$K, "K")
    if (!is.null(mcmc_obj)) mcmc_list$K <- mcmc_obj
  }
  
  if (!is.null(res$alpha)) {
    mcmc_obj <- safe_mcmc(res$alpha, "alpha")
    if (!is.null(mcmc_obj)) mcmc_list$alpha <- mcmc_obj
  }
  
  if (!is.null(res$kern)) {
    mcmc_obj <- safe_mcmc(res$kern, "kernel_idx")
    if (!is.null(mcmc_obj)) mcmc_list$kernel_idx <- mcmc_obj
  }
  
  # Extract acceptance rates if available
  if (!is.null(res$acc) && is.data.frame(res$acc)) {
    acc_data <- res$acc
    if (nrow(acc_data) > 1) {
      if ("accL" %in% names(acc_data)) {
        mcmc_obj <- safe_mcmc(acc_data$accL, "acceptance_L")
        if (!is.null(mcmc_obj)) mcmc_list$acceptance_L <- mcmc_obj
      }
      if ("accEta" %in% names(acc_data)) {
        mcmc_obj <- safe_mcmc(acc_data$accEta, "acceptance_eta")
        if (!is.null(mcmc_obj)) mcmc_list$acceptance_eta <- mcmc_obj
      }
      if ("accKer" %in% names(acc_data)) {
        mcmc_obj <- safe_mcmc(acc_data$accKer, "acceptance_kernel")
        if (!is.null(mcmc_obj)) mcmc_list$acceptance_kernel <- mcmc_obj
      }
      if ("accTauB" %in% names(acc_data)) {
        mcmc_obj <- safe_mcmc(acc_data$accTauB, "acceptance_tauB")
        if (!is.null(mcmc_obj)) mcmc_list$acceptance_tauB <- mcmc_obj
      }
    }
  }

  if (length(mcmc_list) == 0) {
    cat("No valid mcmc objects created for diagnostics\n")
    return(invisible(NULL))
  }

  # Create individual diagnostic plots
  for (name in names(mcmc_list)) {
    mcmc_obj <- mcmc_list[[name]]
    if (is.null(mcmc_obj) || !coda::is.mcmc(mcmc_obj) || length(mcmc_obj) < 2) next

    fn <- sprintf("%s_coda_%s.png", outfile_prefix, name)
    png(fn, width = 1200, height = 900, res = 150)
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mfrow = c(2,2), mar = c(4,4,3,1))

    # Trace plot
    try(coda::traceplot(mcmc_obj, main = sprintf("Trace: %s", name)), silent = TRUE)

    # Density plot
    try(coda::densplot(mcmc_obj, main = sprintf("Density: %s", name)), silent = TRUE)

    # ACF
    try(coda::autocorr.plot(mcmc_obj, main = sprintf("ACF: %s", name)), silent = TRUE)

    # ESS (removed barplot since ESS is already in text files)

    dev.off()
  }

  # Optional: write a single summary file
  summary_file <- sprintf("%s_coda_summary.txt", outfile_prefix)
  cat("=== CODA DIAGNOSTICS SUMMARY ===\n", file = summary_file)
  for (name in names(mcmc_list)) {
    mcmc_obj <- mcmc_list[[name]]
    if (!is.null(mcmc_obj) && coda::is.mcmc(mcmc_obj)) {
      cat(sprintf("\n%s:\n", name), file = summary_file, append = TRUE)
      cat(sprintf("  Mean: %.4f\n", mean(mcmc_obj)), file = summary_file, append = TRUE)
      cat(sprintf("  SD: %.4f\n", sd(mcmc_obj)), file = summary_file, append = TRUE)
      cat(sprintf("  ESS: %.1f\n", coda::effectiveSize(mcmc_obj)), file = summary_file, append = TRUE)
      cat(sprintf("  Geweke p-value: %.4f\n", coda::geweke.diag(mcmc_obj)$p), file = summary_file, append = TRUE)
    }
  }
  invisible(NULL)
}

# ------------------------ Kernel Parameter Coda Diagnostics --------------

.save_coda_kernel_diagnostics <- function(res, outfile_prefix = "run") {
  if (is.null(res$kernel_params_trace) || length(res$kernel_params_trace) == 0) {
    cat("No kernel parameter traces available for diagnostics\n")
    return(invisible(NULL))
  }
  td <- extract_kernel_param_traces(res)
  if (is.null(td$param_traces) || length(td$param_traces) == 0) {
    cat("No valid kernel parameter traces found\n")
    return(invisible(NULL))
  }

  for (param_name in names(td$param_traces)) {
    x <- as.numeric(td$param_traces[[param_name]])
    x <- x[is.finite(x)]
    if (length(x) < 2) next

    fn <- sprintf("%s_coda_kernel_%s.png", outfile_prefix, param_name)
    png(fn, width = 1200, height = 900, res = 150)
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mfrow = c(2,2), mar = c(4,4,3,1))

    # Trace
    plot(seq_along(x), x, type = "l", main = sprintf("Trace: %s", param_name),
         xlab = "Iteration", ylab = "Value")
    abline(h = mean(x), lty = 2)

    # Density
    try(plot(stats::density(x), main = sprintf("Density: %s", param_name),
             xlab = "Value", ylab = "Density"), silent = TRUE)

    # ACF
    try(acf(x, main = sprintf("ACF: %s", param_name), lag.max = min(50, floor(length(x)/4))), silent = TRUE)

    # ESS (removed barplot since ESS is already in text files)

    dev.off()
  }

  invisible(NULL)
}

# ------------------------ Wavelet Coefficient Coda Diagnostics ------------

.save_coda_wavelet_diagnostics <- function(res, outfile_prefix = "run", metric = "active_ratio") {
  td <- extract_wavelet_coeff_traces(res, metric = metric)
  if (is.null(td) || is.null(td$channel_traces) || length(td$channel_traces) == 0) {
    cat("No wavelet coefficient traces available for diagnostics\n")
    return(invisible(NULL))
  }

  for (nm in names(td$channel_traces)) {
    v <- as.numeric(td$channel_traces[[nm]])
    v <- v[is.finite(v)]
    if (length(v) < 2) next

    fn <- sprintf("%s_coda_wavelet_%s_%s.png", outfile_prefix, metric, nm)
    png(fn, width = 1200, height = 900, res = 150)
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mfrow = c(2,2), mar = c(4,4,3,1))

    # Trace
    plot(seq_along(v), v, type = "l",
         xlab = "Iteration", ylab = metric,
         main = sprintf("Trace: %s (%s)", nm, metric))
    abline(h = mean(v), lty = 2)

    # Density
    try(plot(stats::density(v), main = sprintf("Density: %s (%s)", nm, metric),
             xlab = metric, ylab = "Density"), silent = TRUE)

    # ACF
    try(acf(v, main = sprintf("ACF: %s (%s)", nm, metric), lag.max = min(50, floor(length(v)/4))), silent = TRUE)

    # ESS (removed barplot since ESS is already in text files)

    dev.off()
  }

  invisible(NULL)
}
