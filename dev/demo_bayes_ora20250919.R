# -------------
#   個別のλi
# -------------

## ==== kmin の設定 ====
kmin <- 5                  # n がこれ未満は λ 計算から除外
ok   <- n >= kmin

## ==== 1) y と λ_i ====
y <- (k * (k - 1)) / (n * (n - 1))
y[!ok] <- NA               # n < kmin を除外
#eps <- 1e-8
#y <- pmin(pmax(y, m^2 + eps), 1 - eps)

lam_i_hat <- (m - y) / (y - m^2)
lam_i_hat[!is.finite(lam_i_hat) | lam_i_hat <= 0] <- NA

## ==== 2) ガンマ MoM ====
L <- lam_i_hat[is.finite(lam_i_hat)]
mu  <- mean(L, na.rm=TRUE)
v   <- var(L,  na.rm=TRUE)
alpha_hat <- mu^2 / v
beta_hat  <- mu / v
lam_common <- mu

## ==== 3) p の推定 ====
p_hat_common <- (k + m * lam_common) / (n + lam_common)

lam_i_plug <- lam_i_hat
lam_i_plug[!is.finite(lam_i_plug)] <- lam_common
p_hat_gamma <- (k + m * lam_i_plug) / (n + lam_i_plug)

## ==== 出力確認 ====
out <- data.frame(k, n, ok, lam_i_hat, p_hat_common, p_hat_gamma)
print(out, digits=4)

## λ_i の推定値（lam_i_hat）と、推定したガンマパラメータ
L <- lam_i_hat[is.finite(lam_i_hat)]
mu  <- mean(L, na.rm=TRUE)   # 実測平均（1次モーメント）
v   <- var(L,  na.rm=TRUE)   # 実測分散（2次中心モーメント）

alpha_hat <- mu^2 / v
beta_hat  <- mu / v

## --- 理論モーメント（ガンマ分布） ---
mean_theory <- alpha_hat / beta_hat
var_theory  <- alpha_hat / (beta_hat^2)
mom2_theory <- var_theory + mean_theory^2  # 2次モーメント

## --- 実測モーメント ---
mom2_emp <- mean(L^2, na.rm=TRUE)

## --- 比較テーブル ---
out <- data.frame(
  moment = c("mean (1st)", "var (2nd central)", "2nd raw"),
  empirical = c(mu, v, mom2_emp),
  theoretical = c(mean_theory, var_theory, mom2_theory)
)
print(out, digits=5)

