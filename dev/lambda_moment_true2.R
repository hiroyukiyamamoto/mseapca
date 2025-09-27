## ===== method_moment_fixalpha.R (fixed-parameter safe) =====
## α を固定し、λ_i = s + c * U_i^α の c と s を MoM 目的で推定
## - 戻り値: lambda, c_hat, s_hat, alpha, conv, value, message
## - reg_lambda: lam_ref（既定 U=N-n）からの log 比 L2 正則化
## - c/s を上下同値にすると固定として扱い、最適化から除外

compute_lambda_fixalpha2 <- function(
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
  
  ## 固定パラメータ検出（上下同値なら固定）
  c_fixed <- isTRUE(abs(c_bounds[2] - c_bounds[1]) < .Machine$double.eps^0.5)
  s_fixed <- isTRUE(abs(s_bounds[2] - s_bounds[1]) < .Machine$double.eps^0.5)
  c_const <- if (c_fixed) c_bounds[1] else NA_real_
  s_const <- if (s_fixed) s_bounds[1] else NA_real_
  
  ## 目的関数（c,s は実数域）
  obj_cs <- function(cval, sval){
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
  
  ## --- ケース分岐 ---
  if (c_fixed && s_fixed) {
    ## 両方固定：評価のみ
    c_hat <- c_const; s_hat <- s_const
    lambda <- pmax(s_hat + c_hat * (Uj^a), epsp)
    val <- obj_cs(c_hat, s_hat)
    return(list(lambda=lambda, c_hat=c_hat, s_hat=s_hat, alpha=a,
                conv=0, value=val, message="c and s fixed; no optimization"))
  }
  
  if (c_fixed && !s_fixed) {
    ## 1変数最適化：log s
    obj_1 <- function(par){
      sval <- exp(par[1])
      if (sval < s_bounds[1] || sval > s_bounds[2]) return(1e12)
      obj_cs(c_const, sval)
    }
    par0  <- log(max(s_bounds[1], 1e-6))
    lower <- log(max(s_bounds[1], 1e-12))
    upper <- log(max(s_bounds[2], 1+1e-12))
    fit <- optim(par0, obj_1, method="L-BFGS-B", lower=lower, upper=upper,
                 control=list(maxit=maxit))
    s_hat <- exp(fit$par[1]); c_hat <- c_const
    lambda <- pmax(s_hat + c_hat * (Uj^a), epsp)
    return(list(lambda=lambda, c_hat=c_hat, s_hat=s_hat, alpha=a,
                conv=fit$convergence, value=fit$value, message=fit$message))
  }
  
  if (!c_fixed && s_fixed) {
    ## 1変数最適化：log c
    obj_1 <- function(par){
      cval <- exp(par[1])
      if (cval < c_bounds[1] || cval > c_bounds[2]) return(1e12)
      obj_cs(cval, s_const)
    }
    ## 初期値は median(U) を参考（元実装の流儀）
    c0   <- max(stats::median(Uj, na.rm=TRUE), 1)
    par0 <- log(max(c0, 1e-12))
    lower <- log(max(c_bounds[1], 1e-12))
    upper <- log(max(c_bounds[2], 1+1e-12))
    fit <- optim(par0, obj_1, method="L-BFGS-B", lower=lower, upper=upper,
                 control=list(maxit=maxit))
    c_hat <- exp(fit$par[1]); s_hat <- s_const
    lambda <- pmax(s_hat + c_hat * (Uj^a), epsp)
    return(list(lambda=lambda, c_hat=c_hat, s_hat=s_hat, alpha=a,
                conv=fit$convergence, value=fit$value, message=fit$message))
  }
  
  ## 2変数最適化：log c, log s
  c0   <- max(stats::median(Uj, na.rm=TRUE), 1)
  s0   <- max(min(s_bounds[2], 1e-6), s_bounds[1])
  par0 <- c(log(c0), log(s0))
  obj_2 <- function(par){
    cval <- exp(par[1]); sval <- exp(par[2])
    if (cval < c_bounds[1] || cval > c_bounds[2]) return(1e12)
    if (sval < s_bounds[1] || sval > s_bounds[2]) return(1e12)
    obj_cs(cval, sval)
  }
  lower <- c(log(max(c_bounds[1], 1e-12)), log(max(s_bounds[1], 1e-12)))
  upper <- c(log(max(c_bounds[2], 1+1e-12)), log(max(s_bounds[2], 1+1e-12)))
  fit <- optim(par0, obj_2, method="L-BFGS-B", lower=lower, upper=upper,
               control=list(maxit=maxit))
  
  c_hat  <- exp(fit$par[1]); s_hat <- exp(fit$par[2])
  lambda <- pmax(s_hat + c_hat * (Uj^a), epsp)
  
  list(lambda=lambda, c_hat=c_hat, s_hat=s_hat, alpha=a,
       conv=fit$convergence, value=fit$value, message=fit$message)
}
