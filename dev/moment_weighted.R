## ===== method_moment_weighted.R =====
compute_lambda_moment_weighted <- function(k, n, N, m, ...){
  ## ora_est の weighted と同じ： λ_i = U_i = N_i - n_i
  pmax(N - n, 0)
}
