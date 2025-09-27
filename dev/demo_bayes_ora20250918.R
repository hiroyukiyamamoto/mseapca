# ----------------
#   λが個別の場合
# ----------------

# 軽量MoM近似 with kmin フィルタ
light_mom <- function(k, n, N, m,
                      kmin = 2,
                      eps = 1e-8, lam_floor = 1e-6,
                      use = c("mean","var")) {
  use <- match.arg(use)
  
  # 有効観測（n>=2 ないと y_i が定義できない）
  i <- is.finite(k) & is.finite(n) & is.finite(N) & (n >= 2)
  # kmin フィルタ：有意物質がある程度あるものだけ
  i <- i & (k >= kmin)
  
  ki <- k[i]; ni <- n[i]; Ni <- N[i];mi <- m[i]
  Ui <- pmax(Ni - ni, eps)
  
  # 1) 階乗モーメント観測値
  yi <- (ki * pmax(ki - 1, 0)) / pmax(ni * (ni - 1), 1)
  
  # 2) MoM 粗推定: lambda_tilde = (m - y)/(y - m^2)
  denom <- yi - mi^2
  lam_tilde <- (mi - yi) / denom
  lam_tilde[!is.finite(lam_tilde)] <- NA_real_
  lam_tilde <- pmax(lam_tilde, lam_floor)
  
  # 3) z_i = log(lam_tilde) - log(U_i)
  zi <- log(lam_tilde) - log(Ui)
  zi <- zi[is.finite(zi)]
  
  if (length(zi) < 3) stop("有効な z_i が不足しています。")
  
  zbar <- mean(zi, na.rm = TRUE)
  zvar <- var(zi,  na.rm = TRUE)
  
  # 4) kappa を 1次元で解く
  if (use == "mean") {
    f <- function(kappa) digamma(kappa) - log(kappa) - zbar
  } else {
    f <- function(kappa) trigamma(kappa) - zvar
  }
  
  # 探索レンジ
  lo <- 1e-3; hi <- 1e4
  g_lo <- f(lo); g_hi <- f(hi)
  if (!is.finite(g_lo) || !is.finite(g_hi) || sign(g_lo) == sign(g_hi)) {
    lo <- 1e-6; hi <- 1e6
  }
  
  kappa_hat <- try(uniroot(function(x) f(x), lower = lo, upper = hi)$root, silent = TRUE)
  if (inherits(kappa_hat, "try-error")) {
    kappa_hat <- NA_real_
    warning("kappa の解が見つかりませんでした。")
  }
  
  list(
    kappa_hat   = kappa_hat,
    z_mean      = zbar,
    z_var       = zvar,
    lambda_tilde= lam_tilde,
    used        = use,
    kmin        = kmin,
    idx_used    = which(i)
  )
}

# k>=3 のパスウェイだけ使う
res <- light_mom(k, n, N, m, kmin=5, use="mean")
res$kappa_hat


# 観測と理論の階乗2次モーメントを比較する
# 入力:
#   k, n, N, m: ベクトル/スカラー
#   lam: 各経路iの λ 推定値（例: light_mom() の res$lambda_tilde）
#   idx_use: 解析に使ったインデックス（例: light_mom() の res$idx_used）
# 返り値:
#   data.frame と、簡単な誤差サマリ
moment_check <- function(k, n, N, m, lam, idx_use = NULL,
                         lam_floor = 1e-8, m_clip = 1e-6) {
  
  stopifnot(length(k) == length(n), length(n) == length(N), length(lam) == length(k))
  
  # 安定化
  m <- min(max(m, m_clip), 1 - m_clip)
  lam <- pmax(lam, lam_floor)
  
  # 使用インデックス（未指定なら全体。ただし n>=2 を必須に）
  if (is.null(idx_use)) {
    idx_use <- which(is.finite(k) & is.finite(n) & is.finite(N) & is.finite(lam) & n >= 2)
  }
  
  ki <- k[idx_use]; ni <- n[idx_use]; Ni <- N[idx_use]; lami <- lam[idx_use]
  Ui <- pmax(Ni - ni, 0)
  
  # 観測の階乗2次モーメント
  y_obs <- (ki * pmax(ki - 1, 0)) / (ni * pmax(ni - 1, 1))
  
  # 理論の階乗2次モーメント（Beta-Binomial, m, λ）
  y_the <- (m * (m * lami + 1)) / (lami + 1)
  
  # 参考：一次モーメント確認用
  p_hat <- ki / ni               # 観測一次
  p_the <- rep(m, length(ni))    # 理論一次（一定）
  
  # 重み（“mへの寄せ”の強さの目安）
  w <- lami / (ni + lami)
  
  df <- data.frame(
    idx   = idx_use,
    n     = ni,
    k     = ki,
    U     = Ui,
    lambda= lami,
    y_obs = y_obs,
    y_the = y_the,
    y_diff= y_obs - y_the,
    p_obs = p_hat,
    p_the = p_the,
    w     = w
  )
  
  # サマリ指標
  valid <- is.finite(df$y_obs) & is.finite(df$y_the)
  rmse  <- sqrt(mean((df$y_obs[valid] - df$y_the[valid])^2))
  mae   <- mean(abs(df$y_obs[valid] - df$y_the[valid]))
  corry <- suppressWarnings(cor(df$y_obs[valid], df$y_the[valid]))
  
  summary <- list(
    rmse_y   = rmse,
    mae_y    = mae,
    cor_yobs_ythe = corry,
    mean_p_obs = mean(df$p_obs, na.rm=TRUE),
    m_theory    = m,
    n_used      = sum(valid),
    idx_used    = idx_use
  )
  
  list(table = df, summary = summary)
}


# 1) 先に軽量MoMで lambda_tilde と idx_used を取得
res <- light_mom(k, n, N, m, kmin = 2, use = "mean")

# 2) 観測 vs 理論を比較
chk <- moment_check(k, n, N, m, lam = res$lambda_tilde, idx_use = res$idx_used)

# 表をざっと確認
head(chk$table[, c("idx","n","k","U","lambda","y_obs","y_the","y_diff","p_obs","w")], 10)

# 誤差サマリ
chk$summary

