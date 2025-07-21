rm(list=ls(all=TRUE))

library(mseapca)

data(fasting_mseapca)
SIG <- fasting_mseapca$SIG 		# SIG: 有意な物質
DET <- fasting_mseapca$DET 		# DET:検出された物質
M <- fasting_mseapca$pathway 	# M :パスウェイデータ

B1 <- msea_ora(SIG, DET, M, option = "det") # 検出された物質のみ
B2 <- msea_ora(SIG, DET, M, option = "est_naive")
B3 <- msea_ora_range(SIG, DET, M, option = "bino_naive", nsim=1000)

B4 <- msea_ora(SIG, DET, M, option = "est_weighted") # 検出された物質のみ
B5 <- msea_ora(SIG, DET, M, option = "est_shrink", lambda=5)

b1 <- B1[[1]][,1]
b2 <- B2[[1]][,1]
b3 <- B3[[1]][,2]
b4 <- B4[[1]][,1]
b5 <- B5[[1]][,1]

B <- cbind(b1,b2,b3,b4,b5)
colnames(B) <- c("det","est_naive","bino_naive","est_weighted","est_shrink")

# -------------------------------------

# ora_estを実行して結果を取得
B <- msea_ora(SIG, DET, M, option = "est_naive")

# クロス集計表を取り出す
tables <- B2[[3]]

### 経験ベイズに基づいたlambdaの推定

# ベクトルを初期化
r_vec <- c()
n_vec <- c()

# 各パスウェイごとに r = a / (a + b), n = a + b を計算
for (tab in tables) {
  a <- tab[1, 1]  # 有意かつ検出
  b <- tab[1, 2]  # 非有意かつ検出
  n <- a + b
  if (n > 0) {
    r <- a / n
    r_vec <- c(r_vec, r)
    n_vec <- c(n_vec, n)
  }
}

# method of moments による λ の推定
p_bar <- mean(r_vec)
s2 <- var(r_vec)
n_bar <- mean(n_vec)
lambda_hat <- max(0, (s2 - p_bar * (1 - p_bar) / n_bar) / (p_bar * (1 - p_bar)))

B6 <- msea_ora(SIG, DET, M, option = "est_shrink", lambda=lambda_hat)
b6 <- B6[[1]][,1]

B <- cbind(b1,b2,b3,b4,b5,b6)

colnames(B) <- c("det","est_naive","bino_naive","est_weighted","est_shrink","est_ebayes")