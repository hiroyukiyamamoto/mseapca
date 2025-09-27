## =========================================
## ORA: λ~Gamma を周辺化（ラプラス）→ (α,β) 推定 → 事後予測で k を比較
## 必要: ベクトル k, n / スカラー m（固定, 0< m <1）
## =========================================

## ---- 0) 入力例（mが未計算なら参考：本来は全物質で算出）
# m <- sum(k) / sum(n)

## ---- 1) ユーティリティ
clip01 <- function(p) pmin(pmax(p, 1e-12), 1-1e-12)

## Beta–Binomial のログ尤度（常にスカラー numeric）
bb_log <- function(k, n, m, lam){
  m <- clip01(m)
  as.numeric( lchoose(n, k) +
                lbeta(k + m*lam, n - k + (1-m)*lam) -
                lbeta(m*lam, (1-m)*lam) )
}

## Gamma( shape = τ, rate = τ/μ ) のログ密度（μ=平均, τ=“強さ”）
ga_log_mu_tau <- function(lam, mu, tau){
  a <- tau; b <- tau/mu
  as.numeric( a*log(b) - lgamma(a) + (a-1)*log(lam) - b*lam )
}

## ---- 2) ラプラス近似（t = log λ 空間）
# 事後 log 密度: ℓ_i(t) = log f_BB(k|e^t) + log g_Gamma(e^t) + t
ell_i <- function(t, k, n, m, mu, tau){
  lam <- exp(t)
  as.numeric( bb_log(k,n,m,lam) + ga_log_mu_tau(lam, mu, tau) + t )
}

# t*（モード）と曲率 h（2階微分）を数値微分で取得（スカラー返し）
mode_curv_i <- function(k, n, m, mu, tau, t0 = log(mu + 1e-6), maxit = 25){
  t <- t0
  for (it in 1:maxit){
    eps <- 1e-4
    f0  <- ell_i(t,      k,n,m,mu,tau)
    f1  <- ell_i(t+eps,  k,n,m,mu,tau)
    f_1 <- ell_i(t-eps,  k,n,m,mu,tau)
    g   <- (f1 - f_1) / (2*eps)             # 一階
    h   <- (f1 - 2*f0 + f_1) / (eps^2)      # 二階（負が理想）
    if (any(!is.finite(c(g,h)))) break
    h_use <- ifelse(h < -1e-8, h, -1e-8)    # 凹性を保証
    step  <- - g / h_use
    step  <- max(min(step, 1.0), -1.0)      # ステップ幅制限
    t_new <- t + step
    if (!is.finite(t_new) || abs(t_new - t) < 1e-6){ t <- t_new; break }
    t <- t_new
  }
  list(tstar = as.numeric(t), curv = as.numeric(h)[1])  # σ^2 = -1/curv
}

# ラプラス近似の周辺対数尤度（必ずスカラーを返す）
log_marg_i_laplace <- function(k, n, m, mu, tau){
  mc <- mode_curv_i(k,n,m,mu,tau)
  if (!is.finite(mc$curv) || mc$curv >= 0) return(-1e12)
  val <- ell_i(mc$tstar, k,n,m,mu,tau) + 0.5*log(2*pi/(-mc$curv))
  as.numeric(val)[1]
}

## ---- 3) (μ,τ) の推定（＝ Γの α,β に対応）
# 目的関数：負の周辺対数尤度（kmin 未満は除外）
neg_ll <- function(p, k, n, m, kmin = 2){
  mu  <- exp(p[1])                 # μ>0
  tau <- exp(p[2])                 # τ>0（shape）
  tau <- pmin(tau, 50)             # 借力の強さを抑制（必要なら調整）
  idx <- which(is.finite(n) & n >= kmin)
  if (!length(idx)) return(1e12)
  s <- 0
  for (i in idx){
    li <- log_marg_i_laplace(k[i], n[i], m, mu, tau)
    if (any(!is.finite(li))) return(1e12)   # ★ ベクトル安全
    s <- s - as.numeric(li)[1]              # ★ スカラー化
  }
  penalty <- 1 * tau ## ここの係数を変える影響が大きい
  s + penalty
  #s
}

fit_mu_tau <- function(k, n, m, kmin = 2, start = c(log(5), log(5))){
  m <- clip01(m)
  fit <- optim(start, neg_ll, k = k, n = n, m = m, kmin = kmin,
               method = "L-BFGS-B",
               control = list(maxit = 300, factr = 1e-7))
  mu_hat  <- exp(fit$par[1])
  tau_hat <- pmin(exp(fit$par[2]), 50)
  list(mu = mu_hat, tau = tau_hat,
       alpha_hat = tau_hat, beta_hat = tau_hat/mu_hat,  # 変換：α=τ, β=τ/μ
       value = fit$value, conv = fit$convergence)
}

## ---- 4) 事後予測：各パスウェイの k_new をサンプル → 平均/95%区間と観測を比較
# Beta–Binomial 乱数（混合用）
rbbinom_mix <- function(n_draw, size, m, lambda){
  a <- m*lambda; b <- (1-m)*lambda
  p <- rbeta(n_draw, shape1 = a, shape2 = b)
  rbinom(n_draw, size = size, prob = p)
}

ppc_per_pathway <- function(k, n, m, alpha_hat, beta_hat, kmin = 2, S = 2000){
  m <- clip01(m)
  out <- lapply(seq_along(k), function(i){
    if (!is.finite(n[i]) || n[i] < kmin){
      return(data.frame(idx=i, n=n[i], obs=k[i], mean_pred=NA, q05=NA, q95=NA))
    }
    mu  <- alpha_hat / beta_hat
    tau <- alpha_hat
    mc  <- mode_curv_i(k[i], n[i], m, mu, tau)
    if (any(!is.finite(mc$curv)) || any(mc$curv >= 0)) {
      # フォールバック：事前 λ からの予測
      lam_s <- rgamma(S, shape = alpha_hat, rate = beta_hat)
      k_sim <- rbbinom_mix(S, size = n[i], m = m, lambda = lam_s)
      return(data.frame(idx=i, n=n[i], obs=k[i],
                        mean_pred = mean(k_sim),
                        q05 = unname(quantile(k_sim, 0.05)),
                        q95 = unname(quantile(k_sim, 0.95))))
    }
    sigma2 <- -1/mc$curv
    t_s   <- rnorm(S, mean = mc$tstar, sd = sqrt(sigma2))  # t ~ N(t*, σ^2)
    lam_s <- exp(t_s)                                      # λ = exp(t)
    k_sim <- rbbinom_mix(S, size = n[i], m = m, lambda = lam_s)
    data.frame(idx = i, n = n[i], obs = k[i],
               mean_pred = mean(k_sim),
               q05 = unname(quantile(k_sim, 0.05)),
               q95 = unname(quantile(k_sim, 0.95)))
  })
  do.call(rbind, out)
}

## ---- 5) 実行（推定 → 予測 → 比較）
set.seed(42)

m0 <- as.numeric(m[1])   # ← これを以後ずっと使う

# (a) (μ,τ) を推定 → α,β に変換
fit <- fit_mu_tau(k, n, m0, kmin = 0, start = c(log(5), log(5)))
alpha_hat <- fit$alpha_hat
beta_hat  <- fit$beta_hat
print(list(alpha_hat = alpha_hat, beta_hat = beta_hat, conv = fit$conv, nll = fit$value))

# (b) 事後予測で各パスウェイの k を要約
df_ppc <- ppc_per_pathway(k, n, m, alpha_hat, beta_hat, kmin = 5, S = 4000)
df_ppc <- subset(df_ppc, is.finite(obs) & is.finite(mean_pred))

# (c) 観測 vs 予測（表）と指標
comp <- transform(df_ppc,
                  k_obs  = obs,
                  k_pred = mean_pred,
                  resid  = obs - mean_pred)[, c("idx","n","k_obs","k_pred","resid","q05","q95")]

coverage95 <- mean(with(comp, q05 <= k_obs & k_obs <= q95), na.rm = TRUE)
MAE  <- mean(abs(comp$resid), na.rm = TRUE)
RMSE <- sqrt(mean(comp$resid^2, na.rm = TRUE))
COR  <- suppressWarnings(cor(comp$k_obs, comp$k_pred, use = "complete.obs"))

print(list(coverage95 = coverage95, MAE = MAE, RMSE = RMSE, cor = COR))

## ---- 6) 可視化（任意）
# install.packages("ggplot2") が未導入なら先に
# library(ggplot2)
# ggplot(comp, aes(x = k_pred, y = k_obs)) +
#   geom_abline(slope = 1, intercept = 0, linetype = 2) +
#   geom_point(alpha = 0.6) +
#   labs(x = "Predicted k (posterior predictive mean)", y = "Observed k")
#
# ggplot(comp, aes(x = factor(n), y = resid)) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_boxplot(outlier.alpha = 0.2) +
#   labs(x = "Pathway size n", y = "Observed - Predicted")

# -----------------------------------------------------
lambda_postmean <- function(k, n, m, alpha, beta, rel.tol = 1e-6){
  m <- clip01(as.numeric(m)[1])
  alpha <- as.numeric(alpha)[1]; beta <- as.numeric(beta)[1]
  
  num_fun <- function(t){
    lam <- exp(t)
    bad <- !is.finite(lam) | lam <= 0
    val <- log(lam) +
      bb_log(k, n, m, lam) +
      dgamma(lam, shape=alpha, rate=beta, log=TRUE) +
      t
    val[bad] <- -Inf
    exp(val)
  }
  
  den_fun <- function(t){
    lam <- exp(t)
    bad <- !is.finite(lam) | lam <= 0
    val <- bb_log(k, n, m, lam) +
      dgamma(lam, shape=alpha, rate=beta, log=TRUE) +
      t
    val[bad] <- -Inf
    exp(val)
  }
  
  num <- integrate(num_fun, -Inf, Inf, rel.tol=rel.tol, stop.on.error=FALSE)$value
  den <- integrate(den_fun, -Inf, Inf, rel.tol=rel.tol, stop.on.error=FALSE)$value
  
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_)
  num/den
}


lam_hat <- mapply(lambda_postmean, k, n,
                  MoreArgs=list(m=m, alpha=alpha_hat, beta=beta_hat))
head(lam_hat)
