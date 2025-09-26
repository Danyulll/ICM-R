###############################################################################
# SEMI-SUPERVISED DP–ICM–GP WITH WALKER SLICE, RAY–MALLICK WAVELET SHRINKAGE,
# AND CARLIN–CHIB KERNEL SELECTION (Empirical Bayes warm-start for 5% normals)
# ---------------------------------------------------------------------------
# This file has been refactored into modular components:
#   - simulated_data.R: Data generation and anomaly types
#   - mcmc_functions.R: MCMC sampling and model functions  
#   - run_model.R: Main execution, diagnostics, and configuration
###############################################################################

# Load required packages
suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(invgamma)
  library(waveslim)
  library(coda)
  library(ROCR)
})

# Source the organized modules
source("simulated_data.R")
source("mcmc_functions.R") 
source("run_model.R")
source("plotting_functions.R")

# Set default seed
set.seed(42)

# ------------------------ Quick Start Functions -----------------------------

# Function to run a single anomaly type with default settings
quick_run <- function(anom_type = "spikes", N = 100, P = 16, M = 3, 
                      frac_anom = 0.5, seed = 42, n_iter = 1000, n_burn = 500,
                      alpha_prior = c(10,1),
                      wf = "la8", boundary = "periodic",
                      mh_step_L = 0.03, mh_step_eta = 0.10, mh_step_tauB = 0.15,
                      use_besov_pi = TRUE, use_besov_g = TRUE,
                      revealed_pct = 0.05,  # Percentage of normals to reveal (default 5%)
                      emp_bayes_init_iter = 80,
                      unpin_after_warmstart = FALSE,
                      K_init = 5,
                      besov_c2 = 0.5) {
  cat("Quick run with default settings...\n")
  
  # Calculate revealed indices from percentage
  # First generate data to determine which observations are normal
  set.seed(seed)
  data <- make_anomaly_dataset(anom_type = anom_type, N = N, P = P, M = M, frac_anom = frac_anom, seed = seed)
  normal_indices <- which(data$true_anom == 0)
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
    seed = seed,
    n_iter = n_iter,
    n_burn = n_burn,
    n_thin = 1,
    outdir = "quick_outputs",
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
  return(result)
}

# Function to run all 18 anomaly types with default settings
quick_run_all <- function(N = 100, P = 16, M = 3, frac_anom = 0.5, 
                          seed = 42, n_iter = 1000, n_burn = 500,
                          alpha_prior = c(10,1),
                          wf = "la8", boundary = "periodic",
                          mh_step_L = 0.03, mh_step_eta = 0.10, mh_step_tauB = 0.15,
                          use_besov_pi = TRUE, use_besov_g = TRUE,
                          revealed_pct = 0.05,  # Percentage of normals to reveal (default 5%)
                          emp_bayes_init_iter = 80,
                          unpin_after_warmstart = FALSE,
                          K_init = 5,
                          besov_c2 = 0.5) {
  cat("Quick run all 18 anomaly types with default settings...\n")
  results <- run_all_anomaly_types(
    N = N, P = P, M = M,
    frac_anom = frac_anom,
    seed = seed,
    n_iter = n_iter,
    n_burn = n_burn,
    n_thin = 1,
    outdir = "quick_outputs",
    alpha_prior = alpha_prior,
    wf = wf, boundary = boundary,
    mh_step_L = mh_step_L, mh_step_eta = mh_step_eta, mh_step_tauB = mh_step_tauB,
    use_besov_pi = use_besov_pi, use_besov_g = use_besov_g,
    revealed_pct = revealed_pct,
    emp_bayes_init_iter = emp_bayes_init_iter,
    unpin_after_warmstart = unpin_after_warmstart,
    K_init = K_init,
    besov_c2 = besov_c2
  )
  return(results)
}

# ------------------------ Example Usage -------------------------------------

# Uncomment the lines below to run examples:
options(error = NULL)
# Example 1: Run a single anomaly type
# result <- quick_run(anom_type = "spikes", N = 25, P = 16, M = 3, 
#                       frac_anom = 0.15, seed = 42, n_iter = 100, n_burn = 1,
#                       alpha_prior = c(1,1),
#                       wf = "la8", boundary = "periodic",
#                       mh_step_L = 0.03, mh_step_eta = 0.10, mh_step_tauB = 0.15,
#                       use_besov_pi = TRUE, use_besov_g = TRUE,
#                       revealed_pct = 0.05,
#                       emp_bayes_init_iter = 1,
#                       unpin_after_warmstart = FALSE,
#                       K_init = 2,
#                       besov_c2 = 0.5)


# Example 2: Run all anomaly types
results <- quick_run_all( N = 50, P = 16, M = 3, 
                      frac_anom = 0.15, seed = 42, n_iter = 5000, n_burn = 1000,
                      alpha_prior = c(2,1),
                      wf = "la8", boundary = "periodic",
                      mh_step_L = 0.03, mh_step_eta = 0.10, mh_step_tauB = 0.15,
                      use_besov_pi = TRUE, use_besov_g = TRUE,
                      revealed_pct = 0.15,
                      emp_bayes_init_iter = 250,
                      unpin_after_warmstart = FALSE,
                      K_init = 2,
                      besov_c2 = 0.5)

# Example 3: Run with custom configuration
# config <- list(N = 100, P = 16, M = 3, frac_anom = 0.3, seed = 123, 
#                n_iter = 2000, n_burn = 1000, n_thin = 1, outdir = "custom_outputs")
# results <- main(config)

# Example 4: Generate data only (for exploration)
# data <- make_anomaly_dataset("spikes", N = 100, P = 16, M = 3, frac_anom = 0.5, seed = 42)
# .save_dataset_png(data, "example_dataset.png")

cat("ICM-R refactored code loaded successfully!\n")
cat("Available functions:\n")
cat("  - quick_run(): Run single anomaly type with defaults\n")
cat("  - quick_run_all(): Run ALL 18 anomaly types with defaults\n") 
cat("  - main(): Run with custom configuration\n")
cat("  - make_anomaly_dataset(): Generate synthetic data\n")
cat("  - run_model(): Run MCMC inference\n")
cat("  - dahl_partition(), map_partition(): Get consensus partitions\n")
cat("  - Plotting functions: .save_dataset_png(), .save_clustered_png(), etc.\n")
cat("Key parameters:\n")
cat("  - revealed_pct: Percentage of normal observations to reveal (default 5%)\n")
cat("  - frac_anom: Fraction of observations that are anomalies\n")
cat("  - alpha_prior: Prior parameters for DP concentration parameter\n")
cat("See the source files for more detailed functions.\n")
cat("Files organized as:\n")
cat("  - simulated_data.R: Data generation and anomaly types\n")
cat("  - mcmc_functions.R: MCMC sampling and model functions\n")
cat("  - run_model.R: Main execution and diagnostics\n")
cat("  - plotting_functions.R: All plotting and visualization functions\n")
