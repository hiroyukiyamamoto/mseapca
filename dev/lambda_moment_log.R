## ===== method_moment_log_fixalpha.R =====
## α=1 固定の傾き拘束ログ線形（MoM）
##   log λ_i = a + log(U_i + s)  <=>  λ_i = exp(a) * (U_i + s)
## - 目的関数： (観測二乗偏差 − βBin理論分散)^2 / v を最小化（MoM）
## - reg_lambda： log(λ/lam_ref)^2 の弱正則化（既定 lam_ref = U = N−n）
## - s_bounds：s の範囲（既定 [1e-6, 1e6]）※ s > 0 必須
## - lam_floor：数値安定の下駄
##
## 戻り値：list(lambda, c_hat, s_hat, alpha_hat=1, conv, value, message)

compute_lambda_moment_log_fixalpha <- function(
    k, n, N, m,
    ... # eps, maxit, reg_lambda, lam_ref, s_bounds, a_bounds, lam_floor, init_a, init_s
){
  dots         <- list(...)
  eps          <- if (!is.null(dots$eps))          dots$eps          else 1e-9
  maxit        <- if (!is.null(dots$maxit))        dots$maxit        else 2000
  reg_lambda   <- if (!is.null(dots$reg_lambda))   dots$reg_lambda   else 0
  lam_floor    <- if (!is.null(dots$lam_floor))    dots$lam_floor    else 1e-9
  
  s_bounds     <- if (!is.null(dots$s_bounds))     dots$s_bounds     else c(1e-6, 1e6) # s > 0
  a_bounds     <- if (!is.null(dots$a_bounds))     dots$a_bounds     else c(log(1e-6), log(1e6))
  
  init_a       <- if (!is.null(dots$init_a))       dots$init_a       else NA_real_
  init_s       <- if (!is.null(dots$init_s))       dots$init_s       else NA_real_
  
  ## s_bounds 検証（必須）
  if (length(s_bounds) != 2L || any(!is.finite(s_bounds)) ||
      s_bounds[1] <= 0 || s_bounds[2] <= 0 || s_bounds[1] >= s_bounds[2]) {
    stop("s_bounds must be c(lower, upper) with 0 < lower < upper and both finite.")
  }
  
  ## 整形・チェック
  k <- as.numeric(k); n <- as.numeric(n); N <- as.numeric(N)
  stopifnot(length(k)==length(n), length(n)==length(N))
  if (length(m)==1L) m <- rep(as.numeric(m), length(k)) else m <- as.numeric(m)
  stopifnot(length(m)==length(k))
  stopifnot(all(is.finite(k+n+N+m)))
  stopifnot(all(k >= 0 & n >= 0 & N >= 0 & k <= n & n <= N))
  
  ## m は (0,1) へ軽くクリップ
  m <- pmin(pmax(m, 1e-8), 1-1e-8)
  
  ## 基本量
  U  <- pmax(N - n, 0)
  Uj <- pmax(U, eps)
  
  lam_ref <- if (!is.null(dots$lam_ref)) {
    lr <- if (length(dots$lam_ref)==1L) rep(as.numeric(dots$lam_ref), length(k)) else as.numeric(dots$lam_ref)
    stopifnot(length(lr)==length(k)); pmax(lr, eps)
  } else {
    Uj
  }
  
  ## Beta–Binomial の理論分散
  v_model <- function(lam, n, m){
    den <- pmax(lam + 1, eps)
    n * m * (1 - m) * ((lam + n) / den)
  }
  
  ## 初期値
  if (is.na(init_s)) {
    pos <- Uj[Uj > 0]
    s0 <- if (length(pos)) max(stats::quantile(pos, probs = 0.1, na.rm=TRUE), 1e-6) else 1e-6
  } else s0 <- init_s
  
  a0 <- if (!is.na(init_a)) init_a else {
    mu <- mean(Uj + s0)
    log(pmax(mu, 1e-6))
  }
  par0 <- c(a0, log(pmax(s0, 1e-6)))  # (a, log s)
  
  ## 目的関数（L-BFGS-B）
  obj <- function(par){
    a     <- par[1]
    logs  <- par[2]
    
    ## 境界（a はそのまま、s は log から復元してチェック）
    if (!is.finite(a) || a < a_bounds[1] || a > a_bounds[2]) return(1e12)
    s     <- exp(logs)
    if (!is.finite(s) || s < s_bounds[1] || s > s_bounds[2]) return(1e12)
    
    lam   <- pmax(exp(a) * (Uj + s), lam_floor)
    
    v_th  <- v_model(lam, n, m)
    v_th[!is.finite(v_th)] <- 1e12
    
    resid2 <- (k - n * m)^2 - v_th
    
    ## 標準の MoM 重み（逆分散）。過大重みのクリップ。
    w <- 1 / pmax(v_th, eps)
    w <- pmin(w, 1e6)
    
    val <- sum((sqrt(w) * resid2)^2, na.rm=TRUE) / max(length(k), 1)
    
    ## λ の弱正則化（weighted への寄せ）
    if (reg_lambda > 0) {
      pen <- sum((log(lam) - log(lam_ref))^2, na.rm=TRUE) / max(length(k), 1)
      val <- val + 0.5 * reg_lambda * pen
    }
    
    if (!is.finite(val)) 1e12 else val
  }
  
  ## ===== 修正：境界の上限・下限（s の上限ログに 1+ を入れない）=====
  lower <- c(max(a_bounds[1], log(1e-12)), log(max(s_bounds[1], 1e-12)))
  upper <- c(max(a_bounds[2], log(1e-12)), log(max(s_bounds[2], 1e-12)))
  
  fit <- optim(
    par     = par0,
    fn      = obj,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(maxit = maxit)
  )
  
  a_hat  <- fit$par[1]
  s_hat  <- exp(fit$par[2])
  c_hat  <- exp(a_hat)               # λ = c * (U + s)
  lambda <- pmax(c_hat * (Uj + s_hat), lam_floor)
  
  list(
    lambda     = as.numeric(lambda),
    c_hat      = as.numeric(c_hat),
    s_hat      = as.numeric(s_hat),
    alpha_hat  = 1,
    conv       = fit$convergence,       # 0 なら収束
    value      = fit$value,             # 最終目的値（小さいほど良い）
    message    = fit$message
  )
}

## ===== 使い方の安全な比率確認（警告なし） =====
## fit_log_fix <- compute_lambda_moment_log_fixalpha(
##   k, n, N, m,
##   s_bounds   = c(1e-6, 1e3),
##   reg_lambda = 0.1,
##   lam_floor  = 1e-9
## )
## U <- pmax(N - n, 0)
## summary( fit_log_fix$lambda[U > 0] / U[U > 0] )
