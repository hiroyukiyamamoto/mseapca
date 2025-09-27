rm(list=ls(all=TRUE))

library(mseapca)

source("C:/Users/hyama/Documents/R/bayes_ora/utils.R")
source("C:/Users/hyama/Documents/R/bayes_ora/transforms.R")

## データ読み込み
data(fasting_mseapca); SIG <- unique(fasting_mseapca$SIG)
DET <- unique(fasting_mseapca$DET); M <- fasting_mseapca$pathway


cnt <- .get_counts(SIG, DET, M)
N <- cnt$N; n <- cnt$n; k <- cnt$k

eps <- 1e-9

p_global <- length(unique(SIG))/max(length(unique(DET)),1L)
m_vec <- rep(pmin(pmax(p_global, eps), 1-eps), length(k))
m <- m_vec

lam_ref <- N - n

# -------------
#   周辺尤度
# -------------

source("C:/Users/hyama/Documents/R/bayes_ora/lambda_marglik.R")

# c=1固定、sのみ推定（α=1固定）
fit_off <- compute_lambda_marglik(
  k, n, N, m,
  alpha_fixed = 1,
  c_bounds    = c(1, 1),          # ★ cを固定
  s_bounds    = c(1e-6, 1e6),     # sを推定（極小正〜広め）
  reg_lambda  = 0.1                 # モデル比較なら 0 推奨（後述）
)




# ----------------------------------------------------------

# fit_off$c_hat    # ≈ 1（固定）
# fit_off$s_hat    # 推定された下駄
# range(fit_off$lambda)  # λ = U + s の動き
# 
# 
# # lambdaを入力にして、フィッシャーの正確確率検定のp値を計算する
# lambda_use <- fit_off$lambda
# m_vec <- if (length(m)==1L) rep(m, length(M)) else as.numeric(m)
# a_post <- lambda_use * m_vec + k
# b_post <- lambda_use * (1 - m_vec) + (n - k)
# p_hat <- a_post / (a_post + b_post) # 事後平均（未検出割合を含めた推定割合）
# 
# EB_fisher_p <- fisher_from_phat(p_hat, k, n, m_vec, M, DET, N) ## 2) EB-Fisher を計算
# post_enrich <- 1 - pbeta(m_vec, a_post, b_post) ## 3) ベイズ事後確率も計算
# 
# N_vec <- N
# ## 4) 結果をまとめる
# res2 <- data.frame(
#   pathway = names(M),
#   N = N_vec,
#   n = n,
#   k = k,
#   lambda = lambda_use,
#   p_hat = p_hat,
#   EB_fisher_p = EB_fisher_p,
#   post_prob_enriched = post_enrich
# )
# 
# head(res2)
