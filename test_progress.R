# Test script to verify progress messages work
source("sampler3_semi_supervised.R")

cat("Testing progress messages with a small example...\n")
flush.console()

# Run a very small test to see progress messages
test_results <- run_all_anomaly_types(
  N = 5, P = 8, M = 2,
  base_B = diag(2),
  base_eta = c(0.02, 0.02),
  base_kern = kernels[[1]],
  base_par = list(l_scale = 0.15),
  outroot = "test_progress_output",
  seed = 12345,
  reveal_frac = 0.1,
  n_iter = 3, burn = 1, thin = 1,
  alpha_prior = c(5, 1),
  wf = "la8", J = 3, boundary = "periodic",
  mh_step_L = 0.03, mh_step_eta = 0.10,
  use_besov_pi = TRUE, use_besov_g = TRUE,
  emp_bayes_init_iter = 5,
  unpin_after_warmstart = FALSE,
  K_init = 2,
  besov_c2 = 0.5,
  use_parallel = TRUE,
  n_workers = 2
)

cat("Test completed!\n")
