## ===== method_moment_fixalpha.R =====
## α を固定し、 λ_i = s + c * U_i^α の c と s を推定（MoM 目的関数）
## - 戻り値: lambda, c_hat, s_hat, conv, value
## - reg_lambda: weighted(=N-n) から離れすぎを抑える弱正則化（log 比の2乗）
## - lam_ref を省略すると自動で N-n を使います

compute_lambda_fixalpha <- function(
    k, n, N, m, alpha_fixed,
    ... # eps, maxit, reg_lambda, lam_ref, c_bounds, s_bounds を渡せます
){
  dots       <- list(...)
  eps        <- if (!is.null(dots$eps))        dots$eps        else 1e-9
  maxit      <- if (!is.null(dots$maxit))      dots$maxit      else 2000
  reg_lambda <- if (!is.null(dots$reg_lambda)) dots$reg_lambda else 0
  lam_ref    <- if (!is.null(dots$lam_ref))    dots$lam_ref    else pmax(N - n, 0)
  c_bounds   <- if (!is.null(dots$c_bounds))   dots$c_bounds   else c(1e-9, 1e6)   # c >= 0
  s_bounds   <- if (!is.null(dots$s_bounds))   dots$s_bounds   else c(0,    1e6)   # s >= 0 （s=0固定は c(0,0)）
  
  ## ベクトル整形
  k <- as.numeric(k); n <- as.numeric(n); N <- as.numeric(N)
  stopifnot(length(k)==length(n), length(n)==length(N))
  if (length(m)==1L) m <- rep(as.numeric(m), length(k)) else m <- as.numeric(m)
  stopifnot(length(m)==length(k))
  stopifnot(length(lam_ref)==length(k))
  
  ## 基本量
  U  <- pmax(N - n, 0)
  Uj <- pmax(U, eps)
  a  <- as.numeric(alpha_fixed)
  
  ## MoM 用の理論分散（必要に応じて調整可）
  v_model <- function(lam, n, m){
    den <- pmax(lam + 1, eps)
    n * m * (1 - m) * ((lam + n) / den)
  }
  
  ## 初期値（c は median(U), s は 0 を起点に）
  c0   <- max(stats::median(Uj, na.rm=TRUE), 1)
  s0   <- 0
  par0 <- c(log(c0), log(s0 + 1e-6))  # (log c, log s) で最適化
  
  ## 目的関数（log 変数で L-BFGS-B）
  obj <- function(par){
    logc <- par[1]; logs <- par[2]
    cval <- exp(logc)
    sval <- exp(logs)
    
    lam  <- pmax(sval + cval * (Uj^a), eps)
    
    v_th <- v_model(lam, n, m)
    v_th[!is.finite(v_th)] <- 1e12
    
    resid2 <- (k - n * m)^2 - v_th
    w <- 1 / pmax(v_th, eps)
    w <- pmin(w, 1e6)
    
    val <- sum((sqrt(w) * resid2)^2, na.rm=TRUE)
    val <- val / max(length(k), 1)  # スケール正規化
    
    ## λ の離れすぎ抑制（weighted からの log 比の2乗）
    if (reg_lambda > 0) {
      pen <- sum((log((lam + eps)/(lam_ref + eps)))^2, na.rm=TRUE)
      val <- val + 0.5 * reg_lambda * (pen / max(length(k), 1))
    }
    
    if (!is.finite(val)) 1e12 else val
  }
  
  ## 境界（log で与える）
  lower <- c(log(max(c_bounds[1], 1e-12)), log(max(s_bounds[1], 1e-12)))
  upper <- c(log(max(c_bounds[2], 1+1e-12)), log(max(s_bounds[2], 1+1e-12)))
  
  fit <- optim(
    par     = par0,
    fn      = obj,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(maxit = maxit)
  )
  
  c_hat  <- exp(fit$par[1])
  s_hat  <- exp(fit$par[2])
  lambda <- pmax(s_hat + c_hat * (Uj^a), eps)
  
  list(
    lambda = lambda,
    c_hat  = c_hat,
    s_hat  = s_hat,
    alpha  = a,
    conv   = fit$convergence,
    value  = fit$value
  )
}
