rm(list=ls(all=TRUE))

source("C:/Users/hyama/Documents/R/mseapca/ora_est.R") # 修正版

library(mseapca)

data(fasting_mseapca)
SIG <- fasting_mseapca$SIG 		# SIG: 有意な物質
DET <- fasting_mseapca$DET 		# DET:検出された物質
M <- fasting_mseapca$pathway 	# M :パスウェイデータ

B1 <- msea_ora(SIG, DET, M, option = "det") # 検出された物質のみ
B2 <- msea_ora(SIG, DET, M, option = "est_naive")
B2 <- ora_est(SIG, DET, M) # 修正版

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

# --------------------------------------------------------------

ora_est_adaptive_fisher <- function(SIG, DET, M, alpha = 5) {
  # 全物質リスト：Mに登場する全物質とDETの和集合
  ALL0 <- unique(unlist(M))
  ALL <- unique(c(ALL0, DET))
  
  # ラベリング関数（物質が各セットに含まれるか）
  setlabel <- function(x, sets) {
    sapply(sets, function(s) as.integer(ALL %in% intersect(x, s)))
  }
  
  L1 <- setlabel(ALL, M)     # in pathway
  L2 <- setlabel(DET, M)     # in DET
  Lsig <- setlabel(SIG, M)   # in SIG
  
  # 全検出物質中の有意割合（shrinkingの基本パラメータ）
  p <- length(SIG) / length(DET)
  TAB <- NULL
  
  for (i in seq_along(M)) {
    l1 <- colSums(L1)[i]      # パスウェイ全体の物質数（検出 + 未検出）
    l2 <- colSums(L2)[i]      # パスウェイ内で検出された数
    l3 <- colSums(Lsig)[i]    # パスウェイ内で有意だった数（検出のみ）
    
    # lambdaの計算：adaptiveに設定（小さいl2ほど強くshrink）
    lambda_i <- alpha / (l2 + 1e-6)
    
    # 有意割合 r を shrinkage により推定
    r <- (l3 + lambda_i * p) / (l2 + lambda_i)
    
    # shrinkageされた有意割合 r を、パスウェイ全体（l1）に適用
    a_star <- round(r * l1)  # ← 検出＋未検出を含めた pathway 内の有意数（推定）
    b_star <- l1 - a_star
    
    # パスウェイ外の期待値（ALL - l1）に対して p を適用
    c_star <- max(0, round(length(ALL) * p - a_star))
    d_star <- max(0, round(length(ALL) * (1 - p) - b_star))
    
    mat <- matrix(c(a_star, c_star, b_star, d_star), nrow = 2)
    
    print(mat)
    
    # フィッシャーの正確確率検定（右側）
    pval <- tryCatch(fisher.test(mat, alternative = "greater")$p.value, error = function(e) NA)
    
    TAB <- rbind(TAB, data.frame(
      pathway = names(M)[i],
      total = l1,
      detected = l2,
      significant = l3,
      r = r,
      lambda_i = lambda_i,
      p.value = pval
    ))
  }
  
  list(table = TAB, all_ratio = p, method = "adaptive_fisher")
}


B7 <- ora_est_adaptive_fisher(SIG, DET, M, alpha = 5)
b7 <- B7[[1]]$p.value

# ----------------------------------------------------

B <- cbind(b1,b2,b3,b4,b5,b6,b7)

colnames(B) <- c("det","est_naive","bino_naive","est_weighted","est_shrink","est_ebayes", "adaptive")

write.csv(B, file="C:/R/B.csv")

