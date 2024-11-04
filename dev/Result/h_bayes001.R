###  Hierarchical Bayesian ORA (Over-Representation Analysis) 

# 必要なライブラリの読み込み
library(mseapca)
library(rstan)

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
B <- msea_ora(SIG, DET, M)

# ORAのp値の取得
ora_p_values <- sapply(B$TAB, function(tab) fisher.test(tab, alternative = "greater")$p.value)

# Stanモデルの定義
stan_model <- "
data {
  int<lower=0> P;                 // パスウェイの数
  int<lower=0> N[P];              // 各パスウェイの検出された物質数
  int<lower=0> x[P];              // 各パスウェイの有意物質数
  real<lower=0, upper=1> mu_theta[P];  // 各パスウェイの観測割合
}
parameters {
  real<lower=0.001> sigma_theta;       // 成功確率のばらつき（分散）、全体で共通
  real<lower=0.001> k[P];              // 各パスウェイの成功確率の集中度
  real<lower=0, upper=1> theta[P];     // 各パスウェイの成功確率
}
model {
  sigma_theta ~ normal(0, 1);           // 成功確率のばらつきに対する共通の事前分布

  for (p in 1:P) {
    k[p] ~ gamma(2, sigma_theta);       
    theta[p] ~ beta(mu_theta[p] * k[p] + 0.001, (1 - mu_theta[p]) * k[p] + 0.001);  // 形状パラメータに小さな値を追加
    x[p] ~ binomial(N[p], theta[p]);     // 観測データに基づく尤度
  }
}
"

# 各パスウェイの観測割合（有意物質の割合）を計算
mu_theta_values <- sapply(1:length(M), function(i) sum(SIG %in% M[[i]]) / sum(DET %in% M[[i]]))

# データリストの設定
data_list <- list(
  P = length(M),           # パスウェイの数
  N = sapply(M, length),    # 各パスウェイの検出された物質数
  x = sapply(M, function(pathway) sum(SIG %in% pathway)),  # 各パスウェイの有意物質数
  mu_theta = mu_theta_values  # 各パスウェイの観測割合
)

# Stanでサンプリング
fit <- stan(model_code = stan_model, data = data_list, iter = 1000, chains = 4)

# Stanの結果からサンプルを取得
posterior_samples <- extract(fit)

# 全体の物質数
M_total <- length(ALL)

# パスウェイごとに事後分布からリサンプリングし、p値の信頼区間を計算
p_value_ci_matrix <- matrix(NA, nrow = length(M), ncol = 2)
rownames(p_value_ci_matrix) <- names(M)
colnames(p_value_ci_matrix) <- c("Lower p-value CI", "Upper p-value CI")

for (i in 1:length(M)) {
  # 各パスウェイの検出されていない物質数
  N_detected <- sum(DET %in% M[[i]])           # パスウェイ内の検出された物質数
  N_undetected <- data_list$N[i] - N_detected  # パスウェイ内の検出されていない物質数
  x_detected <- sum(SIG %in% M[[i]])           # パスウェイ内の有意物質数
  
  theta_samples_pathway <- posterior_samples$theta[, i]  # 各パスウェイのthetaサンプル
  
  # Fisherの正確検定を使用してp値を計算
  p_values <- sapply(theta_samples_pathway, function(theta) {
    a_undetected <- round(theta * N_undetected)               # 未検出物質のうち有意なものの数
    b_undetected <- N_undetected - a_undetected               # 未検出物質のうち非有意なものの数
    
    a <- x_detected + a_undetected                            # パスウェイ内の有意物質数
    b <- N_detected - x_detected + b_undetected               # パスウェイ内の非有意物質数
    c <- max(0, round(M_total * (length(SIG) / length(DET))) - a)  # 他パスウェイの有意物質数
    d <- max(0, M_total - a - b - c)                         # 他パスウェイの非有意物質数
    
    tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    fisher.test(tab, alternative = "greater")$p.value
  })
  
  # 各パスウェイのp値の95%信頼区間を取得
  p_value_ci_matrix[i, ] <- quantile(p_values, probs = c(0.025, 0.975), na.rm = TRUE)
}

# 結果の表示
print("Result of ORA (baseline p-values)")
print(ora_p_values)

print("95% CI for ORA p-values based on Bayesian hierarchical model")
print(p_value_ci_matrix)
