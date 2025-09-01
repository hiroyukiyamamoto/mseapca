## 必要ならパッケージとデータを読み込み
# library(mseapca)  # fasting_mseapca がこのパッケージなら
data(fasting_mseapca)

SIG <- fasting_mseapca$SIG
DET <- fasting_mseapca$DET
M   <- fasting_mseapca$pathway  # named list (各パスウェイの候補物質ベクトル)

## --- 従来 ORA（例：すでにお持ちの関数） ---
B <- ora_est(SIG, DET, M, method="weighted")
ora_res <- B$`Result of MSEA (ORA with adjustment)`

## --- weighted = ルール型(λ_i = N_i - n_i) をベイズで書いた版 ---
bayes_weighted <- function(x, n, N, m){
  lambda <- pmax(N - n, 0)
  a_post <- lambda * m + x
  b_post <- lambda * (1 - m) + (n - x)
  p_hat  <- a_post / (a_post + b_post)
  ci_lo  <- qbeta(0.025, a_post, b_post)
  ci_hi  <- qbeta(0.975, a_post, b_post)
  data.frame(N=N, n=n, x=x, lambda=lambda, p_hat=p_hat, ci_lo=ci_lo, ci_hi=ci_hi)
}

## セットごとの N_i, n_i, x_i を作る
# N_i: そのパスウェイに登録されているユニーク物質数
N_vec <- vapply(M, function(s) length(unique(s)), integer(1))
# n_i: 検出された物質（DET）のうち、そのパスウェイに属する数
n_vec <- vapply(M, function(s) sum(DET %in% s), integer(1))
# x_i: 有意物質（SIG）のうち、そのパスウェイに属する数
x_vec <- vapply(M, function(s) sum(SIG %in% s), integer(1))

## 事前平均 m は全体率（シンプルに） m = |SIG| / |DET|
m_global <- length(unique(SIG)) / max(length(unique(DET)), 1L)

bw_res <- bayes_weighted(x_vec, n_vec, N_vec, m = m_global)
bw_res$pathway <- names(M)

## 見やすく整形（パスウェイ名・サイズ・検出・有意・p_hat と 95%CrI）
bw_show <- bw_res[, c("pathway","N","n","x","lambda","p_hat","ci_lo","ci_hi")]
bw_show <- bw_show[order(-bw_show$p_hat), ]

## 出力例
cat("=== ORA (既存) の先頭 ===\n"); print(head(ora_res, 10))
cat("\n=== weighted=ベイズ表記 の先頭 ===\n"); print(head(bw_show, 10))

# ------------------------------------------

## --- ora_est(method="weighted") と p値を完全一致させる再現コード ---

# 1) baseテーブルを ora_det から取得（※ ora_est ではない）
base_out <- ora_det(SIG, DET, M)  # ココがポイント
base_tabs <- base_out$`Contingency tables`

# 2) ALL（ora_est と同じ定義）
lenALL <- length(unique(c(unique(unlist(M)), DET)))

# 3) p, n, r の定義（ora_est と同じ）
p_glob <- length(SIG) / length(DET)    # p = |SIG| / |DET|

# 4) ループで各パスウェイの Fisher を再計算（完全に同じ式）
p_fisher_bayes <- vapply(seq_along(M), function(i) {
  l1 <- N_vec[i]       # = colSums(L1)[i]
  l2 <- n_vec[i]       # = colSums(L2)[i]
  l3 <- x_vec[i]       # = colSums(Lsig)[i]
  n  <- l1 - l2
  r  <- (l3 + n * p_glob) / l1
  
  # ora_est: a,b は base_tab に n*r, n*(1-r) を足して round
  a <- round(base_tabs[[i]][1,1] + n * r)
  b <- round(base_tabs[[i]][1,2] + n * (1 - r))
  
  # ora_est: c,d は ALL を母数に round＋max(0, …)
  c <- max(0, round(lenALL * p_glob     - a))
  d <- max(0, round(lenALL * (1 - p_glob) - b))
  
  tab <- t(matrix(c(a,b,c,d), nrow = 2))
  fisher.test(tab, alternative = "greater")$p.value
}, numeric(1))

# 5) 比較
cmp <- data.frame(
  pathway        = names(M),
  p_weighted_ORA = ora_res[,1],
  p_fisher_bayes = p_fisher_bayes
)

# 一致確認（ほぼ 0 になるはず。丸め誤差で ~1e-15 程度）
print(max(abs(cmp$p_weighted_ORA - cmp$p_fisher_bayes), na.rm = TRUE))
