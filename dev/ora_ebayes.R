## ===== ora_ebayes.R =====
source("C:/Users/hyama/Documents/R/bayes_ora/utils.R")
source("C:/Users/hyama/Documents/R/bayes_ora/transforms.R")
source("C:/Users/hyama/Documents/R/bayes_ora/lambda_marglik.R")
source("C:/Users/hyama/Documents/R/bayes_ora/lambda_moment_log.R")
source("C:/Users/hyama/Documents/R/bayes_ora/lambda_moment_power.R")
source("C:/Users/hyama/Documents/R/bayes_ora/lambda_moment_true.R")
source("C:/Users/hyama/Documents/R/bayes_ora/moment_weighted.R")

ora_ebayes <- function(SIG, DET, M,
                       m_type=c("global","loo"),
                       method=c("moment_weighted","moment_true","moment_log","moment_power",
                                "weighted_calib","marglik","marglik_constrained","marglik_reg"),
                       eps=1e-9,
                       ## weighted_calib
                       weighted_calib_model=c("affine","power","scale"),
                       target_mean=NA_real_, target_var_frac=0.10,
                       power_alpha_range=c(0.5,3.5), power_s=NA_real_,
                       ## marglik
                       lower_alpha=0, reg_alpha=0, alpha0=1, init_c=10, init_alpha=1,
                       ## スケーリング（任意）
                       lambda_transform=c("none","ratio","minmax","logistic","lambda_over_c"),
                       lambda_rescale_policy=c("none","keep_mean","target_mean"),
                       logistic_alpha=0.2, lambda_over_c=10, target_lambda_mean=NA_real_
){
  m_type <- match.arg(m_type)
  method <- match.arg(method)
  lambda_transform      <- match.arg(lambda_transform)
  lambda_rescale_policy <- match.arg(lambda_rescale_policy)
  
  cnt <- .get_counts(SIG, DET, M)
  N_vec <- cnt$N; n_vec <- cnt$n; k_vec <- cnt$k
  
  ## ora_est の weighted と同じ比較をしやすいよう、global の既定は p_global
  if (m_type=="global"){
    p_global <- length(unique(SIG))/max(length(unique(DET)),1L)
    m_vec <- rep(pmin(pmax(p_global, eps), 1-eps), length(k_vec))
  } else {
    m_vec <- .compute_m(k_vec, n_vec, "loo", eps)
  }
  
  ## --- λ をメソッド別関数で生成 ---
  lambda_raw <- switch(method,
                       moment_weighted = compute_lambda_moment_weighted(k_vec, n_vec, N_vec, m_vec),
                       moment_true     = compute_lambda_moment_true    (k_vec, n_vec, N_vec, m_vec),
                       moment_log      = compute_lambda_moment_log     (k_vec, n_vec, N_vec, m_vec),
                       moment_power    = compute_lambda_moment_power   (k_vec, n_vec, N_vec, m_vec),
                       weighted_calib  = compute_lambda_weighted_calib (k_vec, n_vec, N_vec, m_vec,
                                                                        target_mean=target_mean, target_var_frac=target_var_frac,
                                                                        weighted_calib_model=match.arg(weighted_calib_model),
                                                                        power_alpha_range=power_alpha_range, power_s=power_s),
                       marglik         = compute_lambda_marglik        (k_vec, n_vec, N_vec, m_vec,
                                                                        lower_c=1e-9, lower_alpha=lower_alpha, init_c=init_c,
                                                                        init_alpha=init_alpha, reg_alpha=0, alpha0=alpha0),
                       marglik_constrained = compute_lambda_marglik    (k_vec, n_vec, N_vec, m_vec,
                                                                        lower_c=1e-9, lower_alpha=max(0.3, lower_alpha),
                                                                        init_c=init_c, init_alpha=init_alpha, reg_alpha=0, alpha0=alpha0),
                       marglik_reg     = compute_lambda_marglik        (k_vec, n_vec, N_vec, m_vec,
                                                                        lower_c=1e-9, lower_alpha=lower_alpha, init_c=init_c,
                                                                        init_alpha=init_alpha, reg_alpha=0.2, alpha0=alpha0)
  )
  
  ## --- 任意: λ の 0–1 変換と再スケール ---
  lambda01 <- lambda_transform01(lambda_raw, N_vec, n_vec,
                                 transform=lambda_transform,
                                 logistic_alpha=logistic_alpha,
                                 lambda_over_c=lambda_over_c)
  if (lambda_rescale_policy=="keep_mean"){
    s <- mean(lambda_raw, na.rm=TRUE) / pmax(mean(lambda01, na.rm=TRUE), .Machine$double.eps)
    lambda_use <- pmax(lambda01 * s, 0)
  } else if (lambda_rescale_policy=="target_mean"){
    tgt <- if (is.na(target_lambda_mean)) mean(lambda_raw, na.rm=TRUE) else target_lambda_mean
    s <- tgt / pmax(mean(lambda01, na.rm=TRUE), .Machine$double.eps)
    lambda_use <- pmax(lambda01 * s, 0)
  } else {
    lambda_use <- lambda_raw
  }
  
  ## --- 事後と Fisher ---
  a_post <- lambda_use*m_vec + k_vec
  b_post <- lambda_use*(1-m_vec) + (n_vec - k_vec)
  p_hat  <- a_post / (a_post + b_post)
  p_fisher <- fisher_from_phat(p_hat, k_vec, n_vec, m_vec, M, DET, N_vec)
  post_enrich <- 1 - pbeta(m_vec, a_post, b_post)
  
  data.frame(
    pathway = names(M),
    N = N_vec, n = n_vec, k = k_vec,
    m = m_vec,
    lambda_raw = lambda_raw,
    lambda01   = lambda01,
    lambda     = lambda_use,
    p_hat = p_hat, post_enrich = post_enrich,
    p_value = p_fisher,
    stringsAsFactors = FALSE
  )
}
