rm(list=ls(all=TRUE))

library(mseapca)
source("C:/Users/hyama/Documents/R/bayes_ora/ora_ebayes.R")  # これが他の8ファイルを source します

## データ読み込み
data(fasting_mseapca); SIG <- unique(fasting_mseapca$SIG)
DET <- unique(fasting_mseapca$DET); M <- fasting_mseapca$pathway

# -------------------
#   weightedの確認
# ------------------

## weighted（= λ_i = N_i - n_i）
res_w <- ora_ebayes(SIG, DET, M, method="moment_weighted")

## 比率にスケーリング（(N-n)/N）して使う
res_ratio <- ora_ebayes(SIG, DET, M,
                        method="moment_weighted",
                        lambda_transform="ratio", lambda_rescale_policy="none")

### オリジナルのweighted
B <- ora_est(SIG, DET, M, method="weighted")
p <- B$`Result of MSEA (ORA with adjustment)`[,1]

res_w$p_value/p

# ---------------------------------------------
#   経験ベイズ(モーメント法)による割合ベイズ
# ---------------------------------------------

source("C:/Users/hyama/Documents/R/bayes_ora/lambda_moment_true.R")


cnt <- .get_counts(SIG, DET, M)
N <- cnt$N; n <- cnt$n; k <- cnt$k

eps <- 1e-9

p_global <- length(unique(SIG))/max(length(unique(DET)),1L)
m_vec <- rep(pmin(pmax(p_global, eps), 1-eps), length(k))
m <- m_vec

lam_ref <- N - n

## 2) α=1 固定、c と s を同時推定（s>0 も許容）
fit2 <- compute_lambda_fixalpha(k, n, N, m,
                                alpha_fixed = 1,
                                reg_lambda  = 0.01,       # 少しだけ安全弁
                                lam_ref     = lam_ref,
                                s_bounds    = c(0, 1e6),  # s を推定
                                c_bounds    = c(1e-6, 1e6)
)
fit2$c_hat; fit2$s_hat
summary(fit2$lambda / lam_ref)

# lambdaを入力にして、フィッシャーの正確確率検定のp値を計算する
lambda_use <- fit2$lambda
m_vec <- if (length(m)==1L) rep(m, length(M)) else as.numeric(m)
a_post <- lambda_use * m_vec + k
b_post <- lambda_use * (1 - m_vec) + (n - k)
p_hat <- a_post / (a_post + b_post) # 事後平均（未検出割合を含めた推定割合）

EB_fisher_p <- fisher_from_phat(p_hat, k, n, m_vec, M, DET, N) ## 2) EB-Fisher を計算
post_enrich <- 1 - pbeta(m_vec, a_post, b_post) ## 3) ベイズ事後確率も計算

N_vec <- N
## 4) 結果をまとめる
res <- data.frame(
  pathway = names(M),
  N = N_vec,
  n = n,
  k = k,
  lambda = lambda_use,
  p_hat = p_hat,
  EB_fisher_p = EB_fisher_p,
  post_prob_enriched = post_enrich
)

head(res)

### 周辺尤度
source("C:/Users/hyama/Documents/R/bayes_ora/lambda_marglik.R")

fit3 <- compute_lambda_marglik(
  k, n, N, m,
  alpha_fixed = 1,
  reg_lambda  = 0.5,
  lam_ref   = lam_ref,               # 省略で自動
  s_bounds    = c(-1e-6, 1e-6)       # U=0の安定化に有効（小さめ正）
)

# lambdaを入力にして、フィッシャーの正確確率検定のp値を計算する
lambda_use <- fit3$lambda
m_vec <- if (length(m)==1L) rep(m, length(M)) else as.numeric(m)
a_post <- lambda_use * m_vec + k
b_post <- lambda_use * (1 - m_vec) + (n - k)
p_hat <- a_post / (a_post + b_post) # 事後平均（未検出割合を含めた推定割合）

EB_fisher_p <- fisher_from_phat(p_hat, k, n, m_vec, M, DET, N) ## 2) EB-Fisher を計算
post_enrich <- 1 - pbeta(m_vec, a_post, b_post) ## 3) ベイズ事後確率も計算

N_vec <- N
## 4) 結果をまとめる
res2 <- data.frame(
  pathway = names(M),
  N = N_vec,
  n = n,
  k = k,
  lambda = lambda_use,
  p_hat = p_hat,
  EB_fisher_p = EB_fisher_p,
  post_prob_enriched = post_enrich
)

head(res2)





