## =========================================================
## Beta–Binomial: λ_i = s + c * U_i^α （U_i = N_i - n_i）
## - MoM 版（α固定）: compute_lambda_fixalpha()
## - 周辺尤度版（α固定/自由）: compute_lambda_marglik()
## 共通仕様:
##   * lam_ref の既定は weighted (= N - n)
##   * reg_lambda は log(λ) と log(lam_ref) の二乗距離ペナルティ
##   * s/c が上下同値（実質固定）の場合は、その変数を最適化から除外
##   * 数値安定: m を (0,1) にクリップ、lam_floor で λ に下駄
## =========================================================

## ---------- MoM: alpha 固定 ----------
compute_lambda_fixalpha <- function(
    k, n, N, m, alpha_fixed,
    ... # eps, maxit, reg_lambda, lam_ref, c_bounds, s_bounds
){
  dots       <- list(...)
  eps        <- if (!is.null(dots$eps))        dots$eps        else 1e-9
  maxit      <- if (!is.null(dots$maxit))      dots$maxit      else 2000
  reg_lambda <- if (!is.null(dots$reg_lambda)) dots$reg_lambda else 0
  c_bounds   <- if (!is.null(dots$c_bounds))   dots$c_bounds   else c(1e-9, 1e6)
  s_bounds   <- if (!is.null(dots$s_bounds))   dots$s_bounds   else c(0,    1e6)
  lam_ref_in <- if (!is.null(dots$lam_ref))    dots$lam_ref    else NULL
  
  ## 整形
  k <- as.numeric(k); n <- as.numeric(n); N <- as.numeric(N)
  stopifnot(length(k)==length(n), length(n)==length(N))
  if (length(m)==1L) m <- rep(as.numeric(m), length(k)) else m <- as.numeric(m)
  stopifnot(length(m)==length(k))
  
  ## U と lam_ref
  U  <- pmax(N - n, 0); epsp <- pmax(eps, 1e-12)
  Uj <- pmax(U, epsp)
  lam_ref <- if (is.null(lam_ref_in)) Uj else {
    lr <- if (length(lam_ref_in)==1L) rep(lam_ref_in, length(k)) else as.numeric(lam_ref_in)
    stopifnot(length(lr)==length(k)); pmax(lr, epsp)
  }
  
  a <- as.numeric(alpha_fixed)
  
  ## Beta–Binomial の理論分散（平均 n*m）
  v_model <- function(lam, n, m){
    den <- pmax(lam + 1, epsp)
    n * m * (1 - m) * ((lam + n) / den)
  }
  
  ## 初期値
  c0   <- max(stats::median(Uj, na.rm=TRUE), 1)
  s0   <- 0
  par0 <- c(log(c0), log(s0 + 1e-6))  # (log c, log s)
  
  ## 目的関数
  obj <- function(par){
    logc <- par[1]; logs <- par[2]
    cval <- exp(logc); sval <- exp(logs)
    lam  <- pmax(sval + cval * (Uj^a), epsp)
    v_th <- v_model(lam, n, m); v_th[!is.finite(v_th)] <- 1e12
    
    resid2 <- (k - n * m)^2 - v_th
    w <- 1 / pmax(v_th, epsp); w <- pmin(w, 1e6)
    
    val <- sum((sqrt(w) * resid2)^2, na.rm=TRUE) / max(length(k), 1)
    
    if (reg_lambda > 0) {
      pen <- sum((log(lam) - log(lam_ref))^2, na.rm=TRUE) / max(length(k), 1)
      val <- val + 0.5 * reg_lambda * pen
    }
    if (!is.finite(val)) 1e12 else val
  }
  
  lower <- c(log(max(c_bounds[1], 1e-12)), log(max(s_bounds[1], 1e-12)))
  upper <- c(log(c_bounds[2]),              log(max(s_bounds[2], 1+1e-12)))
  
  fit <- optim(par0, obj, method="L-BFGS-B", lower=lower, upper=upper,
               control=list(maxit=maxit))
  
  c_hat  <- exp(fit$par[1]); s_hat <- exp(fit$par[2])
  lambda <- pmax(s_hat + c_hat * (Uj^a), epsp)
  
  list(
    lambda = lambda,
    c_hat  = c_hat,
    s_hat  = s_hat,
    alpha  = a,
    conv   = fit$convergence,
    value  = fit$value,
    message= fit$message
  )
}

## ---------- Marginal Likelihood: alpha 固定/自由 ----------
compute_lambda_marglik <- function(
    k, n, N, m,
    ... # eps, maxit,
    # reg_alpha, alpha0,
    # reg_lambda, lam_ref,
    # lam_floor,
    # c_bounds, s_bounds, alpha_bounds,
    # init_c, init_s, init_alpha,
    # alpha_fixed
){
  dots          <- list(...)
  eps           <- if (!is.null(dots$eps))           dots$eps           else 1e-9
  maxit         <- if (!is.null(dots$maxit))         dots$maxit         else 2000
  reg_alpha     <- if (!is.null(dots$reg_alpha))     dots$reg_alpha     else 0
  alpha0        <- if (!is.null(dots$alpha0))        dots$alpha0        else 1
  reg_lambda    <- if (!is.null(dots$reg_lambda))    dots$reg_lambda    else 0
  lam_floor     <- if (!is.null(dots$lam_floor))     dots$lam_floor     else 1e-6
  c_bounds      <- if (!is.null(dots$c_bounds))      dots$c_bounds      else c(1e-9, 1e6)
  s_bounds      <- if (!is.null(dots$s_bounds))      dots$s_bounds      else c(0,    1e6)
  alpha_bounds  <- if (!is.null(dots$alpha_bounds))  dots$alpha_bounds  else c(0,    4)
  init_c        <- if (!is.null(dots$init_c))        dots$init_c        else 10
  init_s        <- if (!is.null(dots$init_s))        dots$init_s        else 0
  init_alpha    <- if (!is.null(dots$init_alpha))    dots$init_alpha    else 1
  alpha_fixed   <- if (!is.null(dots$alpha_fixed))   dots$alpha_fixed   else NULL
  lam_ref_in    <- if (!is.null(dots$lam_ref))       dots$lam_ref       else NULL
  
  ## 整形
  k <- as.numeric(k); n <- as.numeric(n); N <- as.numeric(N)
  stopifnot(length(k)==length(n), length(n)==length(N))
  if (length(m)==1L) m <- rep(as.numeric(m), length(k)) else m <- as.numeric(m)
  stopifnot(length(m)==length(k))
  stopifnot(all(k >= 0 & n >= 0 & N >= 0 & k <= n & n <= N))
  
  ## m は (0,1) へ軽くクリップ
  m <- pmin(pmax(m, 1e-8), 1-1e-8)
  
  ## U と lam_ref
  U  <- pmax(N - n, 0); epsp <- pmax(eps, 1e-12)
  Uj <- pmax(U, epsp)
  lam_ref <- if (is.null(lam_ref_in)) pmax(Uj, lam_floor) else {
    lr <- if (length(lam_ref_in)==1L) rep(lam_ref_in, length(k)) else as.numeric(lam_ref_in)
    stopifnot(length(lr)==length(k)); pmax(lr, lam_floor)
  }
  
  ## s/c が固定かどうかの検出（上下同値なら固定）
  s_fixed <- isTRUE(abs(s_bounds[2] - s_bounds[1]) < .Machine$double.eps^0.5)
  c_fixed <- isTRUE(abs(c_bounds[2] - c_bounds[1]) < .Machine$double.eps^0.5)
  s_const <- if (s_fixed) s_bounds[1] else NA_real_
  c_const <- if (c_fixed) c_bounds[1] else NA_real_
  
  ## 周辺対数尤度
  loglik <- function(c_, s_, a_){
    if (!is.finite(c_) || !is.finite(s_) || !is.finite(a_)) return(-Inf)
    if (!c_fixed && (c_ < c_bounds[1] || c_ > c_bounds[2])) return(-Inf)
    if (!s_fixed && (s_ < s_bounds[1] || s_ > s_bounds[2])) return(-Inf)
    lam <- pmax(s_ + c_ * (Uj^a_), lam_floor)
    a0  <- lam * m; b0 <- lam * (1 - m)
    val <- sum(lchoose(n, k) + lbeta(a0 + k, b0 + n - k) - lbeta(a0, b0))
    if (!is.finite(val)) -Inf else val
  }
  
  add_penalty <- function(val, lam, a_, alpha_pen=TRUE){
    if (!is.finite(val)) return(1e12)
    if (alpha_pen && reg_alpha > 0 && is.finite(a_)) {
      val <- val + 0.5 * reg_alpha * (a_ - alpha0)^2
    }
    if (reg_lambda > 0) {
      log_ratio <- log(lam) - log(lam_ref)
      val <- val + 0.5 * reg_lambda * sum(log_ratio * log_ratio) / max(length(lam), 1)
    }
    if (!is.finite(val)) 1e12 else val
  }
  
  ## ========== α 固定 ==========
  if (!is.null(alpha_fixed)) {
    a_const <- as.numeric(alpha_fixed)
    
    if (s_fixed && !c_fixed) {
      ## ---- s 固定: 1変数 (log c) ----
      nll1 <- function(par){
        logc <- par[1]
        c_ <- exp(logc); s_ <- s_const
        ll <- loglik(c_, s_, a_const)
        lam <- pmax(s_ + c_ * (Uj^a_const), lam_floor)
        add_penalty(-ll, lam, a_const, alpha_pen=FALSE)
      }
      lower <- log(pmax(c_bounds[1], 1e-12))
      upper <- log(c_bounds[2])
      par0  <- log(pmax(init_c, 1e-12))
      
      fit <- optim(par0, nll1, method="L-BFGS-B",
                   lower=lower, upper=upper, control=list(maxit=maxit))
      
      c_hat <- exp(fit$par[1]); s_hat <- s_const; a_hat <- a_const
      
    } else if (c_fixed && !s_fixed) {
      ## ---- c 固定: 1変数 (log s) ----
      nll1 <- function(par){
        logs <- par[1]
        c_ <- c_const; s_ <- exp(logs)
        ll <- loglik(c_, s_, a_const)
        lam <- pmax(s_ + c_ * (Uj^a_const), lam_floor)
        add_penalty(-ll, lam, a_const, alpha_pen=FALSE)
      }
      lower <- log(pmax(s_bounds[1], 1e-12))
      upper <- log(max(s_bounds[2], 1+1e-12))
      par0  <- log(pmax(max(init_s, 1e-9), 1e-12))
      
      fit <- optim(par0, nll1, method="L-BFGS-B",
                   lower=lower, upper=upper, control=list(maxit=maxit))
      
      c_hat <- c_const; s_hat <- exp(fit$par[1]); a_hat <- a_const
      
    } else if (s_fixed && c_fixed) {
      ## ---- 両方固定: 評価のみ ----
      c_hat <- c_const; s_hat <- s_const; a_hat <- a_const
      lam   <- pmax(s_hat + c_hat * (Uj^a_hat), lam_floor)
      val   <- add_penalty(-loglik(c_hat, s_hat, a_hat), lam, a_hat, alpha_pen=FALSE)
      return(list(lambda=lam, c_hat=c_hat, s_hat=s_hat, alpha_hat=a_hat,
                  conv=0, value=val, message="s and c fixed; no optimization"))
    } else {
      ## ---- 通常: 2変数 (log c, log s) ----
      nll2 <- function(par){
        logc <- par[1]; logs <- par[2]
        c_ <- exp(logc); s_ <- exp(logs)
        ll <- loglik(c_, s_, a_const)
        lam <- pmax(s_ + c_ * (Uj^a_const), lam_floor)
        add_penalty(-ll, lam, a_const, alpha_pen=FALSE)
      }
      lower <- log(pmax(c(c_bounds[1], s_bounds[1]), 1e-12))
      upper <- log(c(c_bounds[2], s_bounds[2]))
      par0  <- log(pmax(c(init_c, max(init_s, 1e-9)), 1e-12))
      
      fit <- optim(par0, nll2, method="L-BFGS-B",
                   lower=lower, upper=upper, control=list(maxit=maxit))
      
      c_hat <- exp(fit$par[1]); s_hat <- exp(fit$par[2]); a_hat <- a_const
    }
    
    lambda <- pmax(s_hat + c_hat * (Uj^a_hat), lam_floor)
    return(list(
      lambda=lambda, c_hat=c_hat, s_hat=s_hat, alpha_hat=a_hat,
      conv=fit$convergence, value=fit$value, message=fit$message
    ))
  }
  
  ## ========== α 自由: 3変数 (log c, log s, log α) ==========
  nll3 <- function(par){
    logc <- par[1]; logs <- par[2]; loga <- par[3]
    c_ <- if (c_fixed) c_const else exp(logc)
    s_ <- if (s_fixed) s_const else exp(logs)
    a_ <- exp(loga)
    if (a_ < alpha_bounds[1] || a_ > alpha_bounds[2]) return(1e12)
    ll <- loglik(c_, s_, a_)
    lam <- pmax(s_ + c_ * (Uj^a_), lam_floor)
    add_penalty(-ll, lam, a_, alpha_pen=TRUE)
  }
  
  ## 変数の次元を c/s 固定に応じて減らしたい場合はここを拡張
  lower <- log(pmax(c(c_bounds[1], s_bounds[1], alpha_bounds[1]), 1e-12))
  upper <- log(c(c_bounds[2], s_bounds[2], alpha_bounds[2]))
  par0  <- log(pmax(c(
    if (c_fixed) init_c else init_c,
    if (s_fixed) max(init_s, 1e-9) else max(init_s, 1e-9),
    init_alpha), 1e-12))
  
  fit <- optim(par0, nll3, method="L-BFGS-B",
               lower=lower, upper=upper, control=list(maxit=maxit))
  
  c_hat <- if (c_fixed) c_const else exp(fit$par[1])
  s_hat <- if (s_fixed) s_const else exp(fit$par[2])
  a_hat <- exp(fit$par[3])
  lambda <- pmax(s_hat + c_hat * (Uj^a_hat), lam_floor)
  
  list(
    lambda=lambda, c_hat=c_hat, s_hat=s_hat, alpha_hat=a_hat,
    conv=fit$convergence, value=fit$value, message=fit$message
  )
}

## ---------- 簡単な使い方例 ----------
## α=1 で weighted に寄せたい／U=0 の安定化：
# fit_ml <- compute_lambda_marglik(
#   k, n, N, m,
#   alpha_fixed = 1,
#   reg_lambda  = 0.01,        # 安全弁
#   # lam_ref   = N - n,       # 省略で自動
#   s_bounds    = c(1e-6, 1e-5), # ゼロ幅は避けて極小正
#   lam_floor   = 1e-6
# )
