## ===== method_weighted_calib.R =====
## build_lambda_from_weighted() を使って U を“マゲる”
compute_lambda_weighted_calib <- function(k, n, N, m,
                                          target_mean = NA_real_, target_var_frac = 0.10,
                                          weighted_calib_model = c("affine","power","scale"),
                                          power_alpha_range = c(0.5,3.5), power_s = NA_real_, ...
){
  weighted_calib_model <- match.arg(weighted_calib_model)
  U <- pmax(N-n,0)
  tm <- if (is.na(target_mean)) mean(U) else target_mean
  ps <- if (is.na(power_s)) NULL else power_s
  out <- build_lambda_from_weighted(U, target_mean=tm,
                                    target_var_frac=target_var_frac,
                                    model=weighted_calib_model,
                                    power_alpha_range=power_alpha_range,
                                    power_s=ps)
  out$lambda
}
