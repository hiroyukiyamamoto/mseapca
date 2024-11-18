# Globaltestを計算する関数（Zも計算）
globaltest_manual <- function(X, Y, metabolites, pathway) {
  # パスウェイに含まれる代謝物のみを抽出し、行列に変換
  X_pathway <- as.matrix(X[, colnames(X) %in% pathway, drop = FALSE])
  
  # サンプル間共分散行列 R
  m <- ncol(X_pathway) # パスウェイ内の代謝物数
  R <- (1 / m) * (X_pathway %*% t(X_pathway))
  
  # ロジスティックリンク関数を使用して μ を計算
  alpha <- log(mean(Y) / (1 - mean(Y))) # ロジスティックモデルの切片
  mu <- 1 / (1 + exp(-alpha))           # 期待値 μ
  
  # 残差ベクトル
  residual <- Y - mu
  
  # スコアベクトル Z
  Z <- t(X_pathway) %*% residual
  
  # スコア統計量 Q
  Q <- as.numeric(t(residual) %*% R %*% residual) / (mu^2)
  
  # 結果をリスト形式で返す
  return(list(Q = Q, Z = Z))
}

# パスウェイ全体を解析する関数
run_globaltest <- function(X, Y, metabolites, pathways) {
  results <- list()
  
  for (pathway_name in names(pathways)) {
    pathway <- pathways[[pathway_name]]
    result <- globaltest_manual(X, Y, metabolites, pathway)
    results[[pathway_name]] <- result
  }
  
  return(results)
}

# データの設定
X <- fasting$X
metabolites <- colnames(X)
M <- pathway$fasting

# 表現型データ
Y <- fasting$Y[, 1] # バイナリ表現型の1列目を使用

# Globaltestの実行
results <- run_globaltest(X, Y, metabolites, M)

# 結果の表示
for (pathway_name in names(results)) {
  cat("\nPathway:", pathway_name, "\n")
  cat("Q:", results[[pathway_name]]$Q, "\n")
  cat("Z:", results[[pathway_name]]$Z, "\n")
}

# -----------------------------

# 全物質を使ってスコアベクトル Z を計算する関数
calculate_global_Z <- function(X, Y) {
  # ロジスティックリンク関数を使用して μ を計算
  alpha <- log(mean(Y) / (1 - mean(Y))) # ロジスティックモデルの切片
  mu <- 1 / (1 + exp(-alpha))           # 期待値 μ
  
  # 残差ベクトル
  residual <- Y - mu
  
  # 全物質を使ったスコアベクトル Z の計算
  Z <- t(X) %*% residual
  
  # 結果を返す
  return(Z)
}

# データの設定
X <- fasting$X
Y <- fasting$Y[, 1] # バイナリ表現型の1列目を使用

X <- scale(X)

# 全物質を使った Z の計算
Z_all <- calculate_global_Z(X, Y)

# Z の出力
cat("Global Z (all metabolites):\n")
print(Z_all)
