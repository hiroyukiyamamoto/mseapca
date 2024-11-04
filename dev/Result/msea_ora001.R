###  ORA (Over-Representation Analysis) 

# 必要なライブラリの読み込み
library(mseapca)

# ORA計算用の関数
source("C:/Users/hyama/Documents/R/msea/msea_ora.R")  # msea_ora関数が含まれているスクリプト

# データの読み込み
data(fasting)
data(pathway)

# PCAの実行とロードの取得
pca <- prcomp(fasting$X, scale = TRUE)
pca <- pca_loading(pca)

# 検出された代謝物と有意代謝物の情報
metabolites <- colnames(fasting$X)
SIG <- metabolites[pca$loading$R[, 1] < 0 & pca$loading$p.value[, 1] < 0.05]
DET <- metabolites
ALL <- unique(unlist(pathway$fasting))
M <- pathway$fasting

# ORAの計算
B0 <- msea_ora(SIG, ALL, M) # MetaboAnalyst

B1 <- msea_ora(SIG, DET, M)

source("C:/Users/hyama/Documents/R/msea/msea_ora_all_combination.R")
B2 <- msea_ora_all_combination(SIG, DET, ALL, M)
B3 <- msea_ora_all_combination(SIG, DET, ALL, M, "range")

source("C:/Users/hyama/Documents/R/msea/msea_ora_binomial_ci.R")
B4 <- msea_ora_binomial_ci(SIG, DET, ALL, M)

source("C:/Users/hyama/Documents/R/msea/msea_ora_beta_ci.R")
B5 <- msea_ora_beta_ci(SIG, DET, ALL, M)


