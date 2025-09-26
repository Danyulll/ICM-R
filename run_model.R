###############################################################################
# MAIN EXECUTION SCRIPT AND DIAGNOSTIC FUNCTIONS
# ---------------------------------------------------------------------------
# This file contains the main execution script with different configuration
# options and all diagnostic functions for the semi-supervised DP-ICM-GP model.
###############################################################################

# Source the required files
source("simulated_data.R")
source("mcmc_functions.R")

# Initialize kernels
kernels <- make_kernels()

# ------------------------ Coda-based MCMC Diagnostics ------------------

# Convert MCMC results to coda mcmc objects for better diagnostics
convert_to_coda <- function(res) {
  mcmc_list <- list()

  # Helper function to safely create mcmc objects
  safe_mcmc <- function(data, name) {
    if (is.null(data) || length(data) == 0) {
      return(NULL)
    }

    # Clean the data
    clean_data <- data[is.finite(data)]
    if (length(clean_data) == 0) {
      cat(sprintf("Warning: No finite values in %s, skipping\n", name))
      return(NULL)
    }

    if (length(clean_data) < 2) {
      cat(sprintf("Warning: Insufficient data in %s (length=%d), skipping\n", name, length(clean_data)))
      return(NULL)
    }

    tryCatch(
      {
        coda::mcmc(clean_data)
      },
      error = function(e) {
        cat(sprintf("Warning: Failed to create mcmc object for %s: %s\n", name, e$message))
        NULL
      }
    )
  }

  # Extract and clean main parameters
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

  # Add acceptance rates if available
  if (!is.null(res$acc) && nrow(res$acc) > 0) {
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

  return(mcmc_list)
}

# Comprehensive coda-based diagnostics
.save_coda_diagnostics <- function(res, outfile_prefix = "run") {
  mcmc_list <- convert_to_coda(res)
  if (length(mcmc_list) == 0) {
    cat("No MCMC samples available for diagnostics\n")
    return(invisible(NULL))
  }

  for (name in names(mcmc_list)) {
    mcmc_obj <- mcmc_list[[name]]
    if (is.null(mcmc_obj) || !coda::is.mcmc(mcmc_obj) || length(mcmc_obj) < 2) next

    fn <- sprintf("%s_coda_%s.png", outfile_prefix, name)
    png(fn, width = 1200, height = 900, res = 150)
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

    # Trace + Density
    try(plot(mcmc_obj, main = sprintf("Trace: %s", name)), silent = TRUE)
    try(coda::densplot(mcmc_obj, main = sprintf("Density: %s", name)), silent = TRUE)

    # ACF
    try(coda::autocorr.plot(mcmc_obj, main = sprintf("ACF: %s", name)), silent = TRUE)

    # ESS (removed barplot since ESS is already in text files)

    dev.off()
  }

  # Optional: write a single summary file
  summary_file <- sprintf("%s_coda_summary.txt", outfile_prefix)
  cat("Coda MCMC Diagnostics Summary\n============================\n", file = summary_file)
  for (name in names(mcmc_list)) {
    mcmc_obj <- mcmc_list[[name]]
    if (!is.null(mcmc_obj) && coda::is.mcmc(mcmc_obj)) {
      cat(sprintf("\n%s:\n", name), file = summary_file, append = TRUE)
      cat(sprintf("  Mean: %.4f\n", mean(mcmc_obj)), file = summary_file, append = TRUE)
      cat(sprintf("  SD: %.4f\n", sd(mcmc_obj)), file = summary_file, append = TRUE)
      cat(sprintf("  ESS: %.1f\n", coda::effectiveSize(mcmc_obj)), file = summary_file, append = TRUE)
      # Geweke / Heidel-Welch (best-effort)
      suppressWarnings(try(
        {
          gw <- coda::geweke.diag(mcmc_obj)
          z <- if (is.list(gw) && !is.null(gw$z)) gw$z else gw
          cat(sprintf("  Geweke z (mean): %.4f\n", mean(as.numeric(z), na.rm = TRUE)), file = summary_file, append = TRUE)
        },
        silent = TRUE
      ))
      suppressWarnings(try(
        {
          hw <- coda::heidel.diag(mcmc_obj)
          pass <- tryCatch(hw[1, "stest"], error = function(e) NA)
          cat(sprintf("  Heidel-Welch passed: %s\n", ifelse(isTRUE(pass), "Yes", "No/NA")), file = summary_file, append = TRUE)
        },
        silent = TRUE
      ))
    }
  }

  invisible(NULL)
}

# ------------------------ Classification Metrics and AUC Plots ------------------

# Calculate classification metrics
calculate_classification_metrics <- function(true_labels, predicted_labels) {
  # Convert to binary (0 = normal, 1 = anomaly)
  true_binary <- as.numeric(true_labels)
  pred_binary <- as.numeric(predicted_labels)

  # Calculate confusion matrix components
  tp <- sum(true_binary == 1 & pred_binary == 1) # True Positives
  tn <- sum(true_binary == 0 & pred_binary == 0) # True Negatives
  fp <- sum(true_binary == 0 & pred_binary == 1) # False Positives
  fn <- sum(true_binary == 1 & pred_binary == 0) # False Negatives

  # Calculate metrics
  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  recall <- if (tp + fn > 0) tp / (tp + fn) else 0
  specificity <- if (tn + fp > 0) tn / (tn + fp) else 0
  f1_score <- if (precision + recall > 0) 2 * (precision * recall) / (precision + recall) else 0
  accuracy <- (tp + tn) / (tp + tn + fp + fn)

  # Calculate AUC using ROCR package
  auc_score <- NA
  if (length(unique(true_binary)) > 1 && length(unique(pred_binary)) > 1) {
    tryCatch(
      {
        if (requireNamespace("ROCR", quietly = TRUE)) {
          pred_obj <- ROCR::prediction(pred_binary, true_binary)
          auc_obj <- ROCR::performance(pred_obj, "auc")
          auc_score <- as.numeric(auc_obj@y.values[[1]])
        }
      },
      error = function(e) {
        cat("Warning: AUC calculation failed:", e$message, "\n")
      }
    )
  }

  return(list(
    tp = tp, tn = tn, fp = fp, fn = fn,
    precision = precision,
    recall = recall,
    specificity = specificity,
    f1_score = f1_score,
    accuracy = accuracy,
    auc = auc_score
  ))
}

# Plotting functions are now in plotting_functions.R

# Save classification metrics and AUC plots
.save_classification_metrics <- function(true_anom, pred_anom_dahl, pred_anom_map,
                                         outfile_prefix = "run") {
  if (is.null(true_anom) || is.null(pred_anom_dahl) || is.null(pred_anom_map)) {
    cat("Missing classification data for metrics\n")
    return(invisible(NULL))
  }

  # Calculate metrics for DAHL
  metrics_dahl <- calculate_classification_metrics(true_anom, pred_anom_dahl)
  metrics_map <- calculate_classification_metrics(true_anom, pred_anom_map)

  # Create AUC plots
  png(sprintf("%s_classification_metrics.png", outfile_prefix),
    width = 1600, height = 800, res = 150
  )

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  # DAHL ROC curve
  auc_dahl <- plot_auc_roc(true_anom, pred_anom_dahl, "DAHL Clustering ROC")

  # MAP ROC curve
  auc_map <- plot_auc_roc(true_anom, pred_anom_map, "MAP Clustering ROC")

  dev.off()

  # Save detailed metrics to file
  metrics_file <- sprintf("%s_classification_metrics.txt", outfile_prefix)
  cat("=== CLASSIFICATION METRICS ===\n", file = metrics_file)
  cat(sprintf("Dataset: %s\n", outfile_prefix), file = metrics_file, append = TRUE)
  cat(sprintf("Total samples: %d\n", length(true_anom)), file = metrics_file, append = TRUE)
  cat(sprintf("True anomalies: %d\n", sum(true_anom)), file = metrics_file, append = TRUE)
  cat(sprintf("True normals: %d\n", sum(1 - true_anom)), file = metrics_file, append = TRUE)

  # DAHL metrics
  cat("\n--- DAHL CLUSTERING ---\n", file = metrics_file, append = TRUE)
  cat(sprintf("Predicted anomalies: %d\n", sum(pred_anom_dahl)), file = metrics_file, append = TRUE)
  cat(sprintf("True Positives: %d\n", metrics_dahl$tp), file = metrics_file, append = TRUE)
  cat(sprintf("True Negatives: %d\n", metrics_dahl$tn), file = metrics_file, append = TRUE)
  cat(sprintf("False Positives: %d\n", metrics_dahl$fp), file = metrics_file, append = TRUE)
  cat(sprintf("False Negatives: %d\n", metrics_dahl$fn), file = metrics_file, append = TRUE)
  cat(sprintf("Accuracy: %.4f\n", metrics_dahl$accuracy), file = metrics_file, append = TRUE)
  cat(sprintf("Precision: %.4f\n", metrics_dahl$precision), file = metrics_file, append = TRUE)
  cat(sprintf("Recall (Sensitivity): %.4f\n", metrics_dahl$recall), file = metrics_file, append = TRUE)
  cat(sprintf("Specificity: %.4f\n", metrics_dahl$specificity), file = metrics_file, append = TRUE)
  cat(sprintf("F1-Score: %.4f\n", metrics_dahl$f1_score), file = metrics_file, append = TRUE)
  cat(sprintf("AUC: %s\n", ifelse(is.na(metrics_dahl$auc), "N/A", sprintf("%.4f", metrics_dahl$auc))),
    file = metrics_file, append = TRUE
  )

  # MAP metrics
  cat("\n--- MAP CLUSTERING ---\n", file = metrics_file, append = TRUE)
  cat(sprintf("Predicted anomalies: %d\n", sum(pred_anom_map)), file = metrics_file, append = TRUE)
  cat(sprintf("True Positives: %d\n", metrics_map$tp), file = metrics_file, append = TRUE)
  cat(sprintf("True Negatives: %d\n", metrics_map$tn), file = metrics_file, append = TRUE)
  cat(sprintf("False Positives: %d\n", metrics_map$fp), file = metrics_file, append = TRUE)
  cat(sprintf("False Negatives: %d\n", metrics_map$fn), file = metrics_file, append = TRUE)
  cat(sprintf("Accuracy: %.4f\n", metrics_map$accuracy), file = metrics_file, append = TRUE)
  cat(sprintf("Precision: %.4f\n", metrics_map$precision), file = metrics_file, append = TRUE)
  cat(sprintf("Recall (Sensitivity): %.4f\n", metrics_map$recall), file = metrics_file, append = TRUE)
  cat(sprintf("Specificity: %.4f\n", metrics_map$specificity), file = metrics_file, append = TRUE)
  cat(sprintf("F1-Score: %.4f\n", metrics_map$f1_score), file = metrics_file, append = TRUE)
  cat(sprintf("AUC: %s\n", ifelse(is.na(metrics_map$auc), "N/A", sprintf("%.4f", metrics_map$auc))),
    file = metrics_file, append = TRUE
  )

  # Comparison
  cat("\n--- COMPARISON ---\n", file = metrics_file, append = TRUE)
  cat(sprintf("DAHL vs MAP Accuracy: %.4f vs %.4f\n", metrics_dahl$accuracy, metrics_map$accuracy),
    file = metrics_file, append = TRUE
  )
  cat(sprintf("DAHL vs MAP F1-Score: %.4f vs %.4f\n", metrics_dahl$f1_score, metrics_map$f1_score),
    file = metrics_file, append = TRUE
  )
  if (!is.na(metrics_dahl$auc) && !is.na(metrics_map$auc)) {
    cat(sprintf("DAHL vs MAP AUC: %.4f vs %.4f\n", metrics_dahl$auc, metrics_map$auc),
      file = metrics_file, append = TRUE
    )
  }

  # Print summary to console
  cat(sprintf("Classification Metrics Summary:\n"))
  cat(sprintf(
    "DAHL - Accuracy: %.3f, F1: %.3f, AUC: %s\n",
    metrics_dahl$accuracy, metrics_dahl$f1_score,
    ifelse(is.na(metrics_dahl$auc), "N/A", sprintf("%.3f", metrics_dahl$auc))
  ))
  cat(sprintf(
    "MAP  - Accuracy: %.3f, F1: %.3f, AUC: %s\n",
    metrics_map$accuracy, metrics_map$f1_score,
    ifelse(is.na(metrics_map$auc), "N/A", sprintf("%.3f", metrics_map$auc))
  ))

  invisible(list(dahl = metrics_dahl, map = metrics_map))
}

# ------------------------ Wavelet Coefficient Diagnostics ------------------

extract_wavelet_coeff_traces <- function(res, metric = c("active_ratio", "active_coeffs", "total_coeffs", "mean_active", "sd_active")) {
  metric <- match.arg(metric)

  if (is.null(res$wavelet_coeff_trace) || length(res$wavelet_coeff_trace) == 0) {
    return(NULL)
  }
  traces <- res$wavelet_coeff_trace
  traces <- traces[!sapply(traces, is.null)]
  if (!length(traces)) {
    return(NULL)
  }

  # channels are the names inside each $coeffs list (e.g., "channel_1", "channel_2", ...)
  channels <- unique(unlist(lapply(traces, function(x) names(x$coeffs))))

  channel_traces <- lapply(channels, function(ch) {
    v <- sapply(traces, function(x) {
      cf <- x$coeffs[[ch]]
      if (is.null(cf)) NA_real_ else as.numeric(cf[[metric]])
    })
    as.numeric(v)
  })
  names(channel_traces) <- channels

  list(
    metric = metric,
    iterations = sapply(traces, function(x) x$iter),
    channel_traces = channel_traces
  )
}

# Function to create diagnostic plots for wavelet coefficients
.save_wavelet_coeff_diagnostics <- function(res, outfile_prefix = "run", metric = "active_ratio") {
  td <- extract_wavelet_coeff_traces(res, metric = metric)
  if (is.null(td) || is.null(td$channel_traces) || length(td$channel_traces) == 0) {
    cat("No wavelet coefficient traces available for diagnostics\n")
    return(invisible(NULL))
  }

  png(sprintf("%s_wavelet_coeff_diagnostics.png", outfile_prefix),
    width = 2000, height = 1600, res = 150
  )
  op <- par(no.readonly = TRUE)
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
  )

  n <- length(td$channel_traces)
  n_cols <- min(4, n)
  n_rows <- ceiling(n / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))

  for (nm in names(td$channel_traces)) {
    v <- as.numeric(td$channel_traces[[nm]])
    v <- v[is.finite(v)]
    if (length(v) < 5) {
      plot.new()
      title(main = sprintf("%s (%s: n<5)", nm, metric))
      next
    }

    # Trace
    plot(seq_along(v), v,
      type = "l",
      xlab = "iteration", ylab = metric,
      main = sprintf("Trace: %s", nm)
    )
    abline(h = mean(v), lty = 2)

    # ACF (second panel per facet)
    if (length(v) > 20) {
      acf(v, main = sprintf("ACF: %s", nm), lag.max = min(50, length(v) / 4))
    } else {
      plot.new()
      title(main = sprintf("ACF: %s (n<20)", nm))
    }
  }

  # optional summary file
  summary_file <- sprintf("%s_wavelet_coeff_summary.txt", outfile_prefix)
  cat("=== WAVELET COEFFICIENT DIAGNOSTICS SUMMARY ===\n", file = summary_file)
  cat(
    sprintf(
      "Metric: %s\nChannels: %s\n",
      td$metric, paste(names(td$channel_traces), collapse = ", ")
    ),
    file = summary_file, append = TRUE
  )
  invisible(NULL)
}

# ------------------------ Kernel Parameter Diagnostics ---------------------

extract_kernel_param_traces <- function(res) {
  if (is.null(res$kernel_params_trace) || length(res$kernel_params_trace) == 0) {
    return(list(traces = NULL, param_names = NULL, kernel_names = NULL))
  }
  traces <- res$kernel_params_trace
  valid_traces <- traces[!sapply(traces, is.null)]
  if (length(valid_traces) == 0) {
    return(list(traces = NULL, param_names = NULL, kernel_names = NULL))
  }
  all_param_names <- unique(unlist(lapply(valid_traces, function(x) names(x$params))))
  kernel_names <- unique(sapply(valid_traces, function(x) x$kern_name))
  param_traces <- list()
  for (param_name in all_param_names) {
    param_values <- sapply(valid_traces, function(x) {
      if (param_name %in% names(x$params)) x$params[[param_name]] else NA_real_
    })
    param_traces[[param_name]] <- param_values[!is.na(param_values)]
  }
  list(
    traces = valid_traces,
    param_traces = param_traces,
    param_names = all_param_names,
    kernel_names = kernel_names,
    iterations = sapply(valid_traces, function(x) x$iter)
  )
}

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
    width = 2000, height = 1600, res = 150
  )
  op <- par(no.readonly = TRUE)
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
  )

  n <- length(td$param_traces)
  n_cols <- min(4, n)
  n_rows <- ceiling(n / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))

  for (param_name in names(td$param_traces)) {
    x <- as.numeric(td$param_traces[[param_name]])
    x <- x[is.finite(x)]
    if (length(x) < 5) {
      plot.new()
      title(main = sprintf("%s (n<5)", param_name))
      next
    }

    # Trace
    plot(seq_along(x), x,
      type = "l",
      xlab = "iteration", ylab = "value",
      main = sprintf("Trace: %s", param_name)
    )
    abline(h = mean(x), lty = 2)

    # ACF
    if (length(x) > 20) {
      acf(x, main = sprintf("ACF: %s", param_name), lag.max = min(50, length(x) / 4))
    } else {
      plot.new()
      title(main = sprintf("ACF: %s (n<20)", param_name))
    }
  }

  # optional summary file
  summary_file <- sprintf("%s_kernel_param_summary.txt", outfile_prefix)
  cat("=== KERNEL PARAMETER DIAGNOSTICS SUMMARY ===\n", file = summary_file)
  cat(sprintf("Parameters: %s\n", paste(names(td$param_traces), collapse = ", ")), file = summary_file, append = TRUE)
  cat(sprintf("Kernels: %s\n", paste(td$kernel_names, collapse = ", ")), file = summary_file, append = TRUE)
  invisible(NULL)
}

# Kernel parameter diagnostics using coda
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
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

    # Trace
    plot(seq_along(x), x,
      type = "l", main = sprintf("Trace: %s", param_name),
      xlab = "Iteration", ylab = "Value"
    )
    abline(h = mean(x), lty = 2)

    # Density
    try(plot(stats::density(x),
      main = sprintf("Density: %s", param_name),
      xlab = "Value", ylab = "Density"
    ), silent = TRUE)

    # ACF
    try(acf(x, main = sprintf("ACF: %s", param_name), lag.max = min(50, floor(length(x) / 4))), silent = TRUE)

    # ESS (removed barplot since ESS is already in text files)

    dev.off()
  }

  invisible(NULL)
}

# Wavelet coefficient diagnostics using coda
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
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

    # Trace
    plot(seq_along(v), v,
      type = "l",
      xlab = "Iteration", ylab = metric,
      main = sprintf("Trace: %s (%s)", nm, metric)
    )
    abline(h = mean(v), lty = 2)

    # Density
    try(plot(stats::density(v),
      main = sprintf("Density: %s (%s)", nm, metric),
      xlab = metric, ylab = "Density"
    ), silent = TRUE)

    # ACF
    try(acf(v, main = sprintf("ACF: %s (%s)", nm, metric), lag.max = min(50, floor(length(v) / 4))), silent = TRUE)

    # ESS (removed barplot since ESS is already in text files)

    dev.off()
  }

  invisible(NULL)
}

# ------------------------ Convergence Diagnostics ---------------------------

# .smooth function moved to plotting_functions.R

.save_diagnostics <- function(res, dahl, outfile_prefix = "run",
                              true_anom = NULL, pred_anom = NULL) {
  Ktr <- res$K
  Ktr <- Ktr[is.finite(Ktr)]
  atr <- res$alpha
  atr <- atr[is.finite(atr)]
  ktr <- res$kern
  ktr <- ktr[is.finite(ktr)]
  acc <- res$acc
  PSM <- dahl$PSM
  zhat <- dahl$z_hat

  png(sprintf("%s_diagnostics.png", outfile_prefix), width = 1800, height = 1200, res = 150)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(2, 3), mar = c(3.5, 4, 2.5, 1), mgp = c(2.1, 0.6, 0))
  if (length(Ktr)) {
    plot(seq_along(Ktr), Ktr, type = "s", xlab = "kept iter", ylab = "K", main = "Trace of K")
  } else {
    plot.new()
    title("Trace of K (no draws)")
  }
  if (length(atr)) {
    plot(seq_along(atr), atr, type = "l", xlab = "kept iter", ylab = expression(alpha), main = "Trace of alpha")
  } else {
    plot.new()
    title("Trace of alpha (no draws)")
  }
  if (length(ktr)) {
    plot(seq_along(ktr), ktr, type = "s", xlab = "kept iter", ylab = "kernel idx", main = "Kernel index (largest cluster)")
  } else {
    plot.new()
    title("Kernel index (no draws)")
  }
  if (!is.null(acc) && nrow(acc) > 0) {
    plot(acc$iter, .smooth(acc$accL), type = "l", xlab = "iter", ylab = "acc rate", main = "Acceptance: L", ylim = c(0, 1))
    plot(acc$iter, .smooth(acc$accEta), type = "l", xlab = "iter", ylab = "acc rate", main = "Acceptance: eta", ylim = c(0, 1))
    plot(acc$iter, .smooth(acc$accKer), type = "l", xlab = "iter", ylab = "acc rate", main = "Acceptance: kernel", ylim = c(0, 1))
  } else {
    plot.new()
    title("Acceptance rates (no data)")
    plot.new()
    title("Acceptance rates (no data)")
    plot.new()
    title("Acceptance rates (no data)")
  }
  par(op)
  dev.off()

  # Save diagnostic summary
  diag_file <- sprintf("%s_diag_summary.txt", outfile_prefix)
  cat("=== MCMC DIAGNOSTICS SUMMARY ===\n", file = diag_file)
  cat(sprintf("Total iterations: %d\n", nrow(res$Z)), file = diag_file, append = TRUE)
  cat(sprintf("Final K: %d\n", dahl$K_hat), file = diag_file, append = TRUE)
  cat(sprintf("Mean alpha: %.4f\n", mean(atr)), file = diag_file, append = TRUE)
  if (!is.null(acc) && nrow(acc) > 0) {
    cat(sprintf("Mean acceptance L: %.3f\n", mean(acc$accL, na.rm = TRUE)), file = diag_file, append = TRUE)
    cat(sprintf("Mean acceptance eta: %.3f\n", mean(acc$accEta, na.rm = TRUE)), file = diag_file, append = TRUE)
    cat(sprintf("Mean acceptance kernel: %.3f\n", mean(acc$accKer, na.rm = TRUE)), file = diag_file, append = TRUE)
  }
  invisible(NULL)
}

# ------------------------ Main Execution Functions --------------------------

run_one_anomaly_type <- function(anom_type, N, P, M, frac_anom, seed,
                                 n_iter = 2000, n_burn = 1000, n_thin = 1,
                                 outdir = "anomaly_outputs_s3",
                                 alpha_prior = c(10, 1),
                                 wf = "la8", boundary = "periodic",
                                 mh_step_L = 0.03, mh_step_eta = 0.10, mh_step_tauB = 0.15,
                                 use_besov_pi = TRUE, use_besov_g = TRUE,
                                 revealed_idx = integer(0),
                                 emp_bayes_init_iter = 80,
                                 unpin_after_warmstart = FALSE,
                                 K_init = 5,
                                 besov_c2 = 0.5) {
  cat(sprintf("Running anomaly type: %s\n", anom_type))
  cat(sprintf(
    "Parameters: N=%d, P=%d, M=%d, frac_anom=%.2f, seed=%d\n",
    N, P, M, frac_anom, seed
  ))

  # Set seed for reproducibility
  set.seed(seed)

  # Create output directory
  run_name <- sprintf(
    "%s_N%d_P%d_M%d_frac%02d_seed%d",
    anom_type, N, P, M, as.integer(frac_anom * 100), seed
  )
  outdir_run <- file.path(outdir, run_name)
  if (!dir.exists(outdir_run)) {
    dir.create(outdir_run, recursive = TRUE)
  }

  # Generate data
  # Generate data
  cat("Generating data...\n")
  base_B <- diag(M) # ensure MxM matrix
  stopifnot(is.matrix(base_B), nrow(base_B) == M, ncol(base_B) == M)
data <- make_anomaly_dataset(
  N = N,
  P = P,
  M = M,
  frac_anom = frac_anom,
  seed = seed
)


  # After data <- make_anomaly_dataset(...):
  # After data <- make_anomaly_dataset(...):
  P <- nrow(data$Y[[1]])
  t <- if (!is.null(data$t)) data$t else seq(0, 1, length.out = P) # length(t) == P

  # Save dataset plot
  .save_dataset_png(data$Y, t,
    idx_anom = data$idx_anom,
    outfile = file.path(outdir_run, "gen_dataset.png")
  )

  # Run model (second arg is t, not true_anom)
  cat("Running MCMC...\n")
  res <- run_model(
    data$Y, t,
    n_iter = n_iter, burn = n_burn, thin = n_thin,
    alpha_prior = alpha_prior,
    wf = wf, J = log2(P), boundary = boundary,
    mh_step_L = mh_step_L, mh_step_eta = mh_step_eta, mh_step_tauB = mh_step_tauB,
    use_besov_pi = use_besov_pi, use_besov_g = use_besov_g,
    revealed_idx = revealed_idx,
    emp_bayes_init_iter = emp_bayes_init_iter,
    unpin_after_warmstart = unpin_after_warmstart,
    K_init = K_init,
    besov_c2 = besov_c2
  )
  # Clusterings
  dahl <- dahl <- dahl_from_res(res)
  map <- map_from_res(res)

  # Binary mapping via helper in this file (good for >2 clusters)
  bin_dahl <- determine_normal_anomaly_clusters(dahl$z_hat, reveal_idx = revealed_idx)
  bin_map <- determine_normal_anomaly_clusters(map$z_hat, reveal_idx = revealed_idx)

  # Plots that need t
  .save_clustered_png(data$Y, t, dahl$z_hat,
    outfile = file.path(outdir_run, "gen_clustered_dahl.png")
  )
  .save_clustered_png(data$Y, t, map$z_hat,
    outfile = file.path(outdir_run, "gen_clustered_map.png")
  )

  .save_binary_clustered_png(
    data$Y, t,
    pred_anom = bin_dahl$pred_anom,
    normal_cluster = bin_dahl$normal_cluster,
    anomaly_clusters = bin_dahl$anomaly_clusters,
    outfile = file.path(outdir_run, "gen_binary_dahl.png")
  )
  .save_binary_clustered_png(
    data$Y, t,
    pred_anom = bin_map$pred_anom,
    normal_cluster = bin_map$normal_cluster,
    anomaly_clusters = bin_map$anomaly_clusters,
    outfile = file.path(outdir_run, "gen_binary_map.png")
  )

  # Confusion plots (true, pred)
  .save_confusion_png(
    true_anom = data$true_anom,
    pred_anom = bin_dahl$pred_anom,
    outfile = file.path(outdir_run, "gen_confusion_dahl.png")
  )
  .save_confusion_png(
    true_anom = data$true_anom,
    pred_anom = bin_map$pred_anom,
    outfile = file.path(outdir_run, "gen_confusion_map.png")
  )


  # Save diagnostics
  cat("Saving diagnostics...\n")
  .save_diagnostics(res, dahl, file.path(outdir_run, "gen"))
  .save_coda_diagnostics(res, file.path(outdir_run, "gen_coda"))
  .save_coda_kernel_diagnostics(res, file.path(outdir_run, "gen_coda"))
  .save_coda_wavelet_diagnostics(res, file.path(outdir_run, "gen_coda"))
  .save_kernel_param_diagnostics(res, file.path(outdir_run, "gen"))
  .save_wavelet_coeff_diagnostics(res, file.path(outdir_run, "gen"))

  # Calculate and save classification metrics
  cat("Calculating classification metrics...\n")
  metrics <- .save_classification_metrics(
    true_anom        = data$true_anom,
    pred_anom_dahl   = bin_dahl$pred_anom,
    pred_anom_map    = bin_map$pred_anom,
    outfile_prefix   = file.path(outdir_run, "gen")
  )


  # Save metadata
  meta_file <- file.path(outdir_run, "meta.txt")
  cat("=== RUN METADATA ===\n", file = meta_file)
  cat(sprintf("Anomaly type: %s\n", anom_type), file = meta_file, append = TRUE)
  cat(sprintf("N: %d\n", N), file = meta_file, append = TRUE)
  cat(sprintf("P: %d\n", P), file = meta_file, append = TRUE)
  cat(sprintf("M: %d\n", M), file = meta_file, append = TRUE)
  cat(sprintf("Fraction anomalies: %.2f\n", frac_anom), file = meta_file, append = TRUE)
  cat(sprintf("Seed: %d\n", seed), file = meta_file, append = TRUE)
  cat(sprintf("Iterations: %d\n", n_iter), file = meta_file, append = TRUE)
  cat(sprintf("Burn-in: %d\n", n_burn), file = meta_file, append = TRUE)
  cat(sprintf("Thinning: %d\n", n_thin), file = meta_file, append = TRUE)
  cat(sprintf("Final K (Dahl): %d\n", dahl$K_hat), file = meta_file, append = TRUE)
  cat(sprintf("Final K (MAP): %d\n", map$K_hat), file = meta_file, append = TRUE)

  cat(sprintf("Completed run: %s\n", run_name))
  return(list(
    run_name = run_name,
    data = data,
    res = res,
    dahl = dahl,
    map = map,
    metrics = metrics
  ))
}

run_all_anomaly_types <- function(N = 100, P = 16, M = 3, frac_anom = 0.5,
                                  seed = 8746, n_iter = 2000, n_burn = 1000, n_thin = 1,
                                  outdir = "anomaly_outputs_s3",
                                  alpha_prior = c(10, 1),
                                  wf = "la8", boundary = "periodic",
                                  mh_step_L = 0.03, mh_step_eta = 0.10, mh_step_tauB = 0.15,
                                  use_besov_pi = TRUE, use_besov_g = TRUE,
                                  revealed_pct = 0.05,  # Percentage of normals to reveal
                                  emp_bayes_init_iter = 80,
                                  unpin_after_warmstart = FALSE,
                                  K_init = 5,
                                  besov_c2 = 0.5) {
  cat("Running all anomaly types...\n")
  cat(sprintf(
    "Parameters: N=%d, P=%d, M=%d, frac_anom=%.2f, seed=%d\n",
    N, P, M, frac_anom, seed
  ))

  # Available anomaly types - run all of them
  anom_types <- c("spikes", "burst", "varb", "step", "ramp", "fchg", "phase", 
                  "corr", "warp", "drop", "swap", "gsc", "vglob", "vchan", 
                  "vgroup", "vstep", "vbump", "vdir")

  results <- list()
  for (anom_type in anom_types) {
    # Use different seeds for each anomaly type
    current_seed <- seed + which(anom_types == anom_type)

    tryCatch(
      {
      # Calculate revealed indices from percentage for this anomaly type
      set.seed(current_seed)
      temp_data <- make_anomaly_dataset(anom_type = anom_type, N = N, P = P, M = M, frac_anom = frac_anom, seed = current_seed)
      normal_indices <- which(temp_data$true_anom == 0)
      n_revealed <- max(1, round(revealed_pct * length(normal_indices)))
      revealed_idx <- if (n_revealed > 0 && length(normal_indices) > 0) {
        sample(normal_indices, min(n_revealed, length(normal_indices)))
      } else {
        integer(0)
      }
      
      result <- run_one_anomaly_type(
        anom_type = anom_type,
        N = N, P = P, M = M, 
        frac_anom = frac_anom,
        seed = current_seed,
        n_iter = n_iter, n_burn = n_burn, n_thin = n_thin,
        outdir = outdir,
        alpha_prior = alpha_prior,
        wf = wf, boundary = boundary,
        mh_step_L = mh_step_L, mh_step_eta = mh_step_eta, mh_step_tauB = mh_step_tauB,
        use_besov_pi = use_besov_pi, use_besov_g = use_besov_g,
        revealed_idx = revealed_idx,
        emp_bayes_init_iter = emp_bayes_init_iter,
        unpin_after_warmstart = unpin_after_warmstart,
        K_init = K_init,
        besov_c2 = besov_c2
      )
        results[[anom_type]] <- result
      },
      error = function(e) {
        cat(sprintf("Error running %s: %s\n", anom_type, e$message))
      }
    )
  }

  # Save summary
  summary_file <- file.path(outdir, "output.txt")
  cat("=== SUMMARY OF ALL RUNS ===\n", file = summary_file)
  for (anom_type in names(results)) {
    result <- results[[anom_type]]
    cat(sprintf("\n%s:\n", anom_type), file = summary_file, append = TRUE)
    cat(sprintf("  Run name: %s\n", result$run_name), file = summary_file, append = TRUE)
    cat(sprintf("  Final K (Dahl): %d\n", result$dahl$K_hat), file = summary_file, append = TRUE)
    cat(sprintf("  Final K (MAP): %d\n", result$map$K_hat), file = summary_file, append = TRUE)
    if (!is.null(result$metrics)) {
      cat(
        sprintf(
          "  Dahl - Accuracy: %.3f, F1: %.3f\n",
          result$metrics$dahl$accuracy, result$metrics$dahl$f1_score
        ),
        file = summary_file, append = TRUE
      )
      cat(
        sprintf(
          "  MAP  - Accuracy: %.3f, F1: %.3f\n",
          result$metrics$map$accuracy, result$metrics$map$f1_score
        ),
        file = summary_file, append = TRUE
      )
    }
  }

  cat("Completed all runs\n")
  return(results)
}

# ------------------------ Configuration and Main Script ---------------------

# Default configuration
DEFAULT_CONFIG <- list(
  N = 100,
  P = 16,
  M = 3,
  frac_anom = 0.5,
  seed = 8746,
  n_iter = 2000,
  n_burn = 1000,
  n_thin = 1,
  outdir = "anomaly_outputs_s3"
)

# Main execution function
main <- function(config = DEFAULT_CONFIG) {
  cat("Starting ICM-R anomaly detection analysis...\n")
  cat(sprintf("Configuration: %s\n", paste(names(config), config, sep = "=", collapse = ", ")))

  # Run all anomaly types
  results <- do.call(run_all_anomaly_types, config)

  cat("Analysis completed successfully!\n")
  return(results)
}

# Example usage (commented out to avoid automatic execution)
# results <- main()
