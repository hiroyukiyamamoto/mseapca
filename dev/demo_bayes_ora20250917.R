# ----------------
#   λが共通の場合
# ----------------

rm(list=ls(all=TRUE))

library(mseapca)

source("C:/Users/hyama/Documents/R/bayes_ora/utils.R")
source("C:/Users/hyama/Documents/R/bayes_ora/transforms.R")

## データ読み込み
data(fasting_mseapca); SIG <- unique(fasting_mseapca$SIG)
DET <- unique(fasting_mseapca$DET); M <- fasting_mseapca$pathway


cnt <- .get_counts(SIG, DET, M)
N <- cnt$N; n <- cnt$n; k <- cnt$k

m <- length(SIG)/length(DET) # 事前分布の割合p0

# --------------------------------------------------

## 閾値
Kmin <- 5   # 例: k >= 2 の行だけ使う

## フィルタ & 前処理
i <- (k >= Kmin) & (n >= 0) & is.finite(k) & is.finite(n)
stopifnot(any(i))                      # 1行もなければ停止

k0 <- k
n0 <- n

k <- k[i]; n <- n[i]
m0 <- if (length(m)==1) m else m[1]    # 共通 m

## 観測二次モーメント平均
t  <- mean((k*(k-1))/(n*(n-1)), na.rm=TRUE)
#t <- weighted.mean(k*(k-1)/(n*(n-1)), w = n*(n-1), na.rm=TRUE)


lam <- (m0 - t) / (t - m0^2)

## 事後平均 p（凸結合→必ず 0〜1）
p_hat <- (lam/(lam+n))*m0 + (n/(lam+n))*(k/n)

## 出力（最小）
cat("used_rows =", length(n), "  Kmin =", Kmin, "  lambda_hat =", lam, "\n")
print(head(data.frame(n=n, k=k, p_mle=k/n, p_hat=p_hat)))

# --------------------

## 推定された lambda を固定
m0  <- if (length(m)==1) m else m[1]

## 全パスウェイの p_hat を計算
p_hat_all <- (lam/(lam+n0))*m0 + (n0/(lam+n0))*(k0/n0)

## 確認
cat("lambda_hat =", lam, "\n")
cat("p_hat_all range =", range(p_hat_all, na.rm=TRUE), "\n")

head(data.frame(n=n0, k=k0, p_mle=k0/n0, p_hat=p_hat_all), 20)

plot(k0/n0, p_hat_all)
