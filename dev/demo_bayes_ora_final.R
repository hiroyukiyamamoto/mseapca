## 前提: SIG, DET, M が既にあり、.get_counts() と fisher_from_phat() は定義済み

library(mseapca)

cnt <- .get_counts(SIG, DET, M)
N <- cnt$N; n0 <- cnt$n; k0 <- cnt$k
m0 <- length(SIG) / length(DET)               # 事前割合 p0（共通）

## モーメント法で λ 推定（k>=Kmin のみ使用）
Kmin <- 5
i <- (k0 >= Kmin) & is.finite(k0) & is.finite(n0) & (n0 >= 0)
k <- k0[i]; n <- n0[i]
t <- mean((k*(k-1))/(n*(n-1)), na.rm = TRUE)
lam <- (m0 - t) / (t - m0^2)

## 経験ベイズの p_hat（n=0 の NaN を回避）
p_raw <- ifelse(n0 > 0, k0/n0, NA_real_)
p_hat_all <- (lam/(lam+n0))*m0 + (n0/(lam+n0))*p_raw
p_hat_all[!is.finite(p_hat_all)] <- m0
p_hat_all <- pmin(pmax(p_hat_all, 0), 1)

## パスウェイごとの総候補数 N_vec と m の長さ調整
N_vec <- vapply(M, function(s) length(unique(s)), integer(1))
m_vec <- rep(m0, length(M))

## EB-Fisher 片側 p 値と FDR
p_eb  <- fisher_from_phat(p_hat = p_hat_all, k = k0, n = n0,
                          m = m_vec, M = M, DET = DET, N_vec = N_vec)
padj_eb <- p.adjust(p_eb, method = "BH")

## 結果（上位だけ）
res_eb <- data.frame(pathway = names(M), n = n0, k = k0,
                     phat = p_hat_all, pval = p_eb, padj = padj_eb,
                     check.names = FALSE)
res_eb <- res_eb[order(res_eb$padj), ]
head(res_eb, 10)
